#!/usr/bin/env Rscript

# Estimate joint sample scaling from corrected scaledTPM counts and BS syngroups.
# Use TMM or median-of-ratios with all-BS and strict 1:1 reference sets.
# Write factors, scaled matrices and QC to results/02_cross_species_scaling/.
# Estimate factors across species and timepoints; share each Cm factor across B and S.

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(edgeR)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

project_root <- RELEASE_ROOT
sample_sheet <- "0.salmon/sample_info.tsv"
syngroup_file <- "1.normalization/All_pairs_syngroup_w_notes_gene_corrected.csv"
gene_correction_file <- "1.normalization/gene_correction_with_edited_cds_length.tsv"
out_dir <- "1.normalization/results/02_cross_species_scaling"

scaled_count_files <- c(
  Cm = "1.normalization/results/01_tximport/Cm/Cm_scaledTPM_counts.tsv",
  Ls = "1.normalization/results/01_tximport/Ls/Ls_scaledTPM_counts.tsv",
  Ta = "1.normalization/results/01_tximport/Ta/Ta_scaledTPM_counts.tsv"
)

cm_B_source_label <- "BENI"
cm_S_source_label <- "SICY"

min_norm_count <- 5
min_norm_samples_total <- 3
require_detected_in_each_species <- TRUE

tmm_logratio_trim <- 0.30
tmm_sum_trim <- 0.05

reference_sets <- c("all_BS", "strict_1to1")
normalization_methods <- c("TMM", "median_ratio")

stop_on_missing_expression_gene <- TRUE

resolve_path <- function(path, root = project_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(root, path)
}

stop_if_missing <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(label, " not found:\n", paste0("  - ", missing, collapse = "\n"))
  }
}

read_count_matrix <- function(path) {
  x <- read_tsv(
    path,
    col_types = cols(.default = col_double(), GID = col_character()),
    show_col_types = FALSE
  )
  if (!"GID" %in% names(x)) stop("First identifier column must be named GID: ", path)
  if (anyDuplicated(x$GID)) stop("Duplicated GID in expression matrix: ", path)

  mat <- as.matrix(column_to_rownames(as.data.frame(x), "GID"))
  storage.mode(mat) <- "double"
  if (any(!is.finite(mat))) stop("Non-finite expression value in: ", path)
  if (any(mat < 0)) stop("Negative expression value in: ", path)
  mat
}

write_matrix_tsv <- function(mat, path, id_col = "GID") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_tsv(rownames_to_column(as.data.frame(mat), id_col), path)
}

geometric_mean <- function(x) {
  exp(mean(log(x)))
}

center_positive_factors <- function(x) {
  if (any(!is.finite(x)) || any(x <= 0)) stop("Normalization factors must be finite and positive.")
  x / geometric_mean(x)
}

validate_gene_membership <- function(map, gene_col, label) {
  membership <- map %>%
    distinct(comp_id, gene = .data[[gene_col]]) %>%
    dplyr::count(gene, name = "n_comp_id")
  bad <- membership %>% filter(n_comp_id > 1)
  if (nrow(bad) > 0) {
    stop(
      label, " gene(s) occur in more than one comp_id. First IDs: ",
      paste(head(bad$gene, 10), collapse = ", ")
    )
  }
}

check_expression_ids <- function(map_ids, mat, label) {
  missing <- setdiff(unique(map_ids), rownames(mat))
  if (length(missing) > 0) {
    msg <- paste0(
      label, ": ", length(missing),
      " syngroup gene IDs are absent from the expression matrix. First IDs: ",
      paste(head(missing, 10), collapse = ", ")
    )
    if (stop_on_missing_expression_gene) stop(msg) else warning(msg)
  }
  invisible(missing)
}

# Remap obsolete fragments to retained genes using the upstream correction table.
build_gene_correction_map <- function(correction_tbl) {
  required <- c("species", "keep_gene", "merge_genes")
  missing <- setdiff(required, names(correction_tbl))
  if (length(missing) > 0) {
    stop("Gene correction table is missing columns: ", paste(missing, collapse = ", "))
  }

  x <- correction_tbl

  if ("status" %in% names(x)) {
    bad <- x %>% filter(is.na(status) | status != "OK")
    if (nrow(bad) > 0) {
      stop(
        "Gene correction table contains non-OK rows. First keep_gene IDs: ",
        paste(head(bad$keep_gene, 10), collapse = ", ")
      )
    }
  }

  map <- x %>%
    select(species, keep_gene, merge_genes) %>%
    mutate(
      species = trimws(as.character(species)),
      keep_gene = trimws(as.character(keep_gene)),
      merge_genes = as.character(merge_genes)
    ) %>%
    separate_rows(merge_genes, sep = ";") %>%
    mutate(old_gene = trimws(merge_genes)) %>%
    filter(
      !is.na(species), species != "",
      !is.na(keep_gene), keep_gene != "",
      !is.na(old_gene), old_gene != "",
      old_gene != keep_gene
    ) %>%
    select(species, old_gene, new_gene = keep_gene) %>%
    distinct()

  conflict <- map %>%
    distinct(species, old_gene, new_gene) %>%
    count(species, old_gene, name = "n_targets") %>%
    filter(n_targets > 1)

  if (nrow(conflict) > 0) {
    stop(
      "A correction fragment maps to more than one keep_gene. First IDs: ",
      paste(head(conflict$old_gene, 10), collapse = ", ")
    )
  }

  map
}

replace_corrected_gene_ids <- function(x, species, correction_map) {
  x <- as.character(x)
  map_sp <- correction_map %>% filter(.data$species == species)
  if (nrow(map_sp) == 0) return(x)

  idx <- match(x, map_sp$old_gene)
  replace <- !is.na(idx)
  x[replace] <- map_sp$new_gene[idx[replace]]
  x
}

# Count each unique gene once, regardless of repeated synteny edges.
aggregate_unique_gene_mean <- function(mat, map, group_col = "comp_id", gene_col = "gene") {
  map <- map %>%
    transmute(comp_id = as.character(.data[[group_col]]), gene = .data[[gene_col]]) %>%
    filter(!is.na(comp_id), !is.na(gene), gene != "") %>%
    distinct(comp_id, gene)

  groups <- split(map$gene, map$comp_id)
  ans <- t(vapply(groups, function(ids) {
    ids_present <- ids[ids %in% rownames(mat)]
    if (length(ids_present) == 0) return(rep(NA_real_, ncol(mat)))
    colMeans(mat[ids_present, , drop = FALSE])
  }, numeric(ncol(mat))))

  colnames(ans) <- colnames(mat)
  ans
}

median_ratio_factors <- function(counts) {
  # Use geometric means of positive values without pseudocounts.
  row_gm <- apply(counts, 1, function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (length(positive) == 0) NA_real_ else exp(mean(log(positive)))
  })

  ratios <- sweep(counts, 1, row_gm, "/")
  sf <- apply(ratios, 2, function(x) {
    usable <- x[is.finite(x) & x > 0]
    if (length(usable) == 0) NA_real_ else median(usable)
  })
  center_positive_factors(sf)
}

select_normalization_groups <- function(group_mat, sample_species) {
  detected <- group_mat >= min_norm_count
  keep <- rowSums(detected) >= min_norm_samples_total

  if (require_detected_in_each_species) {
    for (sp in unique(sample_species)) {
      sp_cols <- names(sample_species)[sample_species == sp]
      keep <- keep & rowSums(detected[, sp_cols, drop = FALSE]) >= 1
    }
  }
  keep
}

estimate_scaling <- function(group_mat, full_library_sizes, method) {
  if (!identical(colnames(group_mat), names(full_library_sizes))) {
    stop("Column order differs between group matrix and full library sizes.")
  }

  if (method == "TMM") {
    y <- DGEList(counts = group_mat, lib.size = full_library_sizes)
    y <- calcNormFactors(
      y,
      method = "TMM",
      logratioTrim = tmm_logratio_trim,
      sumTrim = tmm_sum_trim
    )

    norm_factor <- y$samples$norm.factors
    effective_library <- y$samples$lib.size * norm_factor
    target_library <- geometric_mean(effective_library)
    multiplier <- target_library / effective_library

    tibble(
      sample = colnames(group_mat),
      raw_full_library_size = as.numeric(y$samples$lib.size),
      norm_factor = norm_factor,
      effective_library_size = effective_library,
      target_library_size = target_library,
      multiplier = multiplier
    )

  } else if (method == "median_ratio") {
    size_factor <- median_ratio_factors(group_mat)
    multiplier <- 1 / size_factor

    tibble(
      sample = colnames(group_mat),
      raw_full_library_size = as.numeric(full_library_sizes),
      norm_factor = size_factor,
      # Median-of-ratios factors already supply the denominator; omit library-size scaling.
      effective_library_size = NA_real_,
      target_library_size = NA_real_,
      multiplier = multiplier
    )

  } else {
    stop("Unknown normalization method: ", method)
  }
}

apply_multipliers <- function(mat, factor_tbl) {
  missing <- setdiff(colnames(mat), factor_tbl$sample)
  if (length(missing) > 0) stop("Factors absent for sample(s): ", paste(missing, collapse = ", "))
  mult <- setNames(factor_tbl$multiplier, factor_tbl$sample)[colnames(mat)]
  sweep(mat, 2, mult, "*")
}

make_output_prefix <- function(reference_set, method) {
  paste(reference_set, method, sep = "__")
}

sample_sheet_path <- resolve_path(sample_sheet)
syngroup_path <- resolve_path(syngroup_file)
gene_correction_path <- resolve_path(gene_correction_file)
scaled_count_paths <- vapply(scaled_count_files, resolve_path, character(1))

stop_if_missing(sample_sheet_path, "Sample sheet")
stop_if_missing(syngroup_path, "Syngroup file")
stop_if_missing(gene_correction_path, "Gene correction file")
stop_if_missing(scaled_count_paths, "scaledTPM matrix/matrices")

samples <- read_tsv(
  sample_sheet_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)
if (!all(c("sample", "species") %in% names(samples))) {
  stop("Sample sheet must contain sample and species columns.")
}
if (anyDuplicated(samples$sample)) stop("Duplicated sample IDs in sample sheet.")

gene_corrections <- read_tsv(
  gene_correction_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)
correction_map <- build_gene_correction_map(gene_corrections)

cm <- read_count_matrix(scaled_count_paths[["Cm"]])
ls <- read_count_matrix(scaled_count_paths[["Ls"]])
ta <- read_count_matrix(scaled_count_paths[["Ta"]])

matrices <- list(Cm = cm, Ls = ls, Ta = ta)
for (sp in names(matrices)) {
  expected <- samples %>% filter(species == sp) %>% pull(sample)
  missing <- setdiff(expected, colnames(matrices[[sp]]))
  extra <- setdiff(colnames(matrices[[sp]]), expected)
  if (length(missing) > 0 || length(extra) > 0) {
    stop(
      "Sample mismatch for ", sp,
      ". Missing from matrix: ", paste(missing, collapse = ", "),
      "; extra in matrix: ", paste(extra, collapse = ", ")
    )
  }
  matrices[[sp]] <- matrices[[sp]][, expected, drop = FALSE]
}
cm <- matrices$Cm
ls <- matrices$Ls
ta <- matrices$Ta

dir.create(resolve_path(out_dir), recursive = TRUE, showWarnings = FALSE)

syngroups_all <- read_csv(
  syngroup_path,
  na = c("", "NA"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

required_group_cols <- c(
  "Group", "Cm1", "Cm2", "Lsi", "Tan",
  "Source1", "Source2", "Type", "comp_id", "complexity"
)
missing_group_cols <- setdiff(required_group_cols, names(syngroups_all))
if (length(missing_group_cols) > 0) {
  stop("Syngroup table is missing columns: ", paste(missing_group_cols, collapse = ", "))
}

bs <- syngroups_all %>%
  filter(complexity == "simple", Type == "BS") %>%
  mutate(comp_id = as.character(comp_id))
if (nrow(bs) == 0) stop("No rows pass complexity == simple & Type == BS.")

correction_usage <- correction_map %>%
  rowwise() %>%
  mutate(
    n_occurrences_before = case_when(
      species == "Cm" ~ sum(bs$Cm1 == old_gene, na.rm = TRUE) + sum(bs$Cm2 == old_gene, na.rm = TRUE),
      species == "Ls" ~ sum(bs$Lsi == old_gene, na.rm = TRUE),
      species == "Ta" ~ sum(bs$Tan == old_gene, na.rm = TRUE),
      TRUE ~ 0L
    )
  ) %>%
  ungroup()

bs <- bs %>%
  mutate(
    Cm1 = replace_corrected_gene_ids(Cm1, "Cm", correction_map),
    Cm2 = replace_corrected_gene_ids(Cm2, "Cm", correction_map),
    Lsi = replace_corrected_gene_ids(Lsi, "Ls", correction_map),
    Tan = replace_corrected_gene_ids(Tan, "Ta", correction_map)
  ) %>%
  distinct()

write_tsv(
  correction_usage,
  file.path(resolve_path(out_dir), "syngroup_gene_correction_audit.tsv")
)

message(
  "Syngroup gene-ID correction: ",
  sum(correction_usage$n_occurrences_before),
  " obsolete fragment occurrence(s) replaced across selected simple BS rows."
)

# Require exactly one Lsi and one Tan gene per retained comp_id.
parent_check <- bs %>%
  group_by(comp_id) %>%
  summarise(
    n_Lsi = n_distinct(Lsi),
    n_Tan = n_distinct(Tan),
    .groups = "drop"
  )
bad_parent <- parent_check %>% filter(n_Lsi != 1 | n_Tan != 1)
if (nrow(bad_parent) > 0) {
  stop("Some retained BS comp_ids do not contain exactly one Lsi and one Tan gene.")
}

cm_long <- bind_rows(
  bs %>% transmute(comp_id, gene = Cm1, source = Source1),
  bs %>% transmute(comp_id, gene = Cm2, source = Source2)
) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(comp_id, gene, source)

bad_source <- cm_long %>%
  filter(is.na(source) | !source %in% c(cm_B_source_label, cm_S_source_label))
if (nrow(bad_source) > 0) {
  stop("Missing or unexpected Source label in retained BS rows.")
}

cm_gene_source_conflict <- cm_long %>%
  distinct(gene, source) %>%
  dplyr::count(gene, name = "n_source") %>%
  filter(n_source > 1)
if (nrow(cm_gene_source_conflict) > 0) {
  stop("A Cm gene has conflicting B/S source assignments.")
}

cm_B_map <- cm_long %>% filter(source == cm_B_source_label) %>% select(comp_id, gene)
cm_S_map <- cm_long %>% filter(source == cm_S_source_label) %>% select(comp_id, gene)
ls_map <- bs %>% distinct(comp_id, gene = Lsi)
ta_map <- bs %>% distinct(comp_id, gene = Tan)

validate_gene_membership(cm_B_map, "gene", "CmB")
validate_gene_membership(cm_S_map, "gene", "CmS")
validate_gene_membership(ls_map, "gene", "Ls")
validate_gene_membership(ta_map, "gene", "Ta")

check_expression_ids(c(cm_B_map$gene, cm_S_map$gene), cm, "Cm")
check_expression_ids(ls_map$gene, ls, "Ls")
check_expression_ids(ta_map$gene, ta, "Ta")

group_structure <- parent_check %>%
  left_join(cm_B_map %>% dplyr::count(comp_id, name = "n_CmB"), by = "comp_id") %>%
  left_join(cm_S_map %>% dplyr::count(comp_id, name = "n_CmS"), by = "comp_id") %>%
  left_join(bs %>% dplyr::count(comp_id, name = "n_synteny_rows"), by = "comp_id") %>%
  mutate(
    strict_1to1 = n_CmB == 1 & n_CmS == 1 & n_Lsi == 1 & n_Tan == 1
  ) %>%
  arrange(as.numeric(comp_id))

write_tsv(group_structure, file.path(resolve_path(out_dir), "BS_comp_id_structure.tsv"))

cm_B_group <- aggregate_unique_gene_mean(cm, cm_B_map)
cm_S_group <- aggregate_unique_gene_mean(cm, cm_S_map)
ls_group   <- aggregate_unique_gene_mean(ls, ls_map)
ta_group   <- aggregate_unique_gene_mean(ta, ta_map)

common_groups <- Reduce(
  intersect,
  list(rownames(cm_B_group), rownames(cm_S_group), rownames(ls_group), rownames(ta_group))
)
if (length(common_groups) == 0) stop("No comp_id remains after group aggregation.")

cm_B_group <- cm_B_group[common_groups, , drop = FALSE]
cm_S_group <- cm_S_group[common_groups, , drop = FALSE]
ls_group   <- ls_group[common_groups, , drop = FALSE]
ta_group   <- ta_group[common_groups, , drop = FALSE]
cm_total_group <- cm_B_group + cm_S_group

# Combine B and S within each Cm library before factor estimation.
group_matrix_all <- cbind(cm_total_group, ls_group, ta_group)

sample_species <- setNames(
  samples$species[match(colnames(group_matrix_all), samples$sample)],
  colnames(group_matrix_all)
)
if (any(is.na(sample_species))) stop("Cannot assign species to every group-matrix sample.")

full_library_sizes <- c(colSums(cm), colSums(ls), colSums(ta))
full_library_sizes <- full_library_sizes[colnames(group_matrix_all)]
if (any(!is.finite(full_library_sizes)) || any(full_library_sizes <= 0)) {
  stop("Invalid full scaledTPM library size.")
}

write_matrix_tsv(group_matrix_all, file.path(resolve_path(out_dir), "BS_all_groups_unscaled_total_matrix.tsv"), "comp_id")
write_matrix_tsv(cm_B_group, file.path(resolve_path(out_dir), "BS_all_groups_unscaled_CmB_mean.tsv"), "comp_id")
write_matrix_tsv(cm_S_group, file.path(resolve_path(out_dir), "BS_all_groups_unscaled_CmS_mean.tsv"), "comp_id")

all_factor_tables <- list()
all_group_filter_qc <- list()

for (reference_set in reference_sets) {
  if (reference_set == "all_BS") {
    candidate_groups <- common_groups
  } else if (reference_set == "strict_1to1") {
    candidate_groups <- group_structure %>%
      filter(strict_1to1, comp_id %in% common_groups) %>%
      pull(comp_id)
  } else {
    stop("Unknown reference set: ", reference_set)
  }

  candidate_mat <- group_matrix_all[candidate_groups, , drop = FALSE]
  keep_for_norm <- select_normalization_groups(candidate_mat, sample_species)
  norm_mat <- candidate_mat[keep_for_norm, , drop = FALSE]

  if (nrow(norm_mat) < 100) {
    stop(
      "Fewer than 100 groups remain for normalization in reference set ",
      reference_set, ". Reconsider filtering thresholds."
    )
  }

  group_filter_qc <- tibble(
    reference_set = reference_set,
    comp_id = rownames(candidate_mat),
    used_for_factor = keep_for_norm,
    n_detected_samples = rowSums(candidate_mat >= min_norm_count),
    detected_Cm = rowSums(candidate_mat[, names(sample_species)[sample_species == "Cm"], drop = FALSE] >= min_norm_count),
    detected_Ls = rowSums(candidate_mat[, names(sample_species)[sample_species == "Ls"], drop = FALSE] >= min_norm_count),
    detected_Ta = rowSums(candidate_mat[, names(sample_species)[sample_species == "Ta"], drop = FALSE] >= min_norm_count)
  )
  write_tsv(
    group_filter_qc,
    file.path(resolve_path(out_dir), paste0(reference_set, "__normalization_group_filter.tsv"))
  )
  all_group_filter_qc[[reference_set]] <- group_filter_qc

  for (method in normalization_methods) {
    prefix <- make_output_prefix(reference_set, method)
    message("Processing ", prefix, " with ", nrow(norm_mat), " normalization groups.")

    factor_tbl <- estimate_scaling(norm_mat, full_library_sizes, method) %>%
      mutate(
        species = unname(sample_species[sample]),
        reference_set = reference_set,
        method = method,
        n_candidate_groups = nrow(candidate_mat),
        n_groups_used = nrow(norm_mat)
      ) %>%
      select(
        reference_set, method, sample, species,
        n_candidate_groups, n_groups_used,
        raw_full_library_size, norm_factor,
        effective_library_size, target_library_size, multiplier
      )

    write_tsv(
      factor_tbl,
      file.path(resolve_path(out_dir), paste0(prefix, "__sample_factors.tsv"))
    )
    all_factor_tables[[prefix]] <- factor_tbl

    cm_scaled <- apply_multipliers(cm, factor_tbl)
    ls_scaled <- apply_multipliers(ls, factor_tbl)
    ta_scaled <- apply_multipliers(ta, factor_tbl)

    write_matrix_tsv(cm_scaled, file.path(resolve_path(out_dir), paste0(prefix, "__Cm_all_genes_scaled.tsv")))
    write_matrix_tsv(ls_scaled, file.path(resolve_path(out_dir), paste0(prefix, "__Ls_all_genes_scaled.tsv")))
    write_matrix_tsv(ta_scaled, file.path(resolve_path(out_dir), paste0(prefix, "__Ta_all_genes_scaled.tsv")))

    # Use the same Cm multiplier for both subgenomes to preserve their balance.
    cm_B_scaled <- apply_multipliers(cm_B_group, factor_tbl)
    cm_S_scaled <- apply_multipliers(cm_S_group, factor_tbl)
    ls_group_scaled <- apply_multipliers(ls_group, factor_tbl)
    ta_group_scaled <- apply_multipliers(ta_group, factor_tbl)
    cm_total_scaled <- cm_B_scaled + cm_S_scaled

    total_group_scaled <- cbind(cm_total_scaled, ls_group_scaled, ta_group_scaled)
    write_matrix_tsv(
      total_group_scaled,
      file.path(resolve_path(out_dir), paste0(prefix, "__BS_group_total_scaled.tsv")),
      "comp_id"
    )

    colnames(cm_B_scaled) <- paste0("CmB__", colnames(cm_B_scaled))
    colnames(cm_S_scaled) <- paste0("CmS__", colnames(cm_S_scaled))
    colnames(ls_group_scaled) <- paste0("Ls__", colnames(ls_group_scaled))
    colnames(ta_group_scaled) <- paste0("Ta__", colnames(ta_group_scaled))
    component_scaled <- cbind(cm_B_scaled, cm_S_scaled, ls_group_scaled, ta_group_scaled)

    write_matrix_tsv(
      component_scaled,
      file.path(resolve_path(out_dir), paste0(prefix, "__BS_group_components_scaled.tsv")),
      "comp_id"
    )

    post_qc <- tibble(
      sample = colnames(total_group_scaled),
      species = unname(sample_species[colnames(total_group_scaled)]),
      column_sum_before = colSums(group_matrix_all[, colnames(total_group_scaled), drop = FALSE]),
      column_sum_after = colSums(total_group_scaled),
      median_before = apply(group_matrix_all[, colnames(total_group_scaled), drop = FALSE], 2, median),
      median_after = apply(total_group_scaled, 2, median),
      reference_set = reference_set,
      method = method
    )
    write_tsv(
      post_qc,
      file.path(resolve_path(out_dir), paste0(prefix, "__post_scaling_QC.tsv"))
    )
  }
}

factor_comparison <- bind_rows(all_factor_tables)
write_tsv(factor_comparison, file.path(resolve_path(out_dir), "all_candidate_sample_factors.tsv"))

factor_wide <- factor_comparison %>%
  mutate(candidate = paste(reference_set, method, sep = "__")) %>%
  select(sample, candidate, multiplier) %>%
  pivot_wider(names_from = candidate, values_from = multiplier)
write_tsv(factor_wide, file.path(resolve_path(out_dir), "all_candidate_multipliers_wide.tsv"))

factor_cor <- factor_wide %>%
  column_to_rownames("sample") %>%
  as.matrix() %>%
  cor(use = "pairwise.complete.obs", method = "pearson")
write_matrix_tsv(
  factor_cor,
  file.path(resolve_path(out_dir), "candidate_multiplier_correlations.tsv"),
  "candidate"
)

summary_qc <- tibble(
  metric = c(
    "input_rows",
    "simple_BS_rows",
    "simple_BS_comp_ids",
    "strict_1to1_comp_ids",
    "multi_copy_comp_ids",
    "unique_Cm_genes",
    "unique_Ls_genes",
    "unique_Ta_genes"
  ),
  value = c(
    nrow(syngroups_all),
    nrow(bs),
    n_distinct(bs$comp_id),
    sum(group_structure$strict_1to1),
    sum(!group_structure$strict_1to1),
    n_distinct(cm_long$gene),
    n_distinct(ls_map$gene),
    n_distinct(ta_map$gene)
  )
)
write_tsv(summary_qc, file.path(resolve_path(out_dir), "BS_syngroup_summary_QC.tsv"))

writeLines(capture.output(sessionInfo()), file.path(resolve_path(out_dir), "sessionInfo.txt"))

message("\nDone. Output directory: ", resolve_path(out_dir))
message("Primary result: all_BS__TMM__BS_group_components_scaled.tsv")
message("Sensitivity results: all_BS__median_ratio and strict_1to1 variants.")
