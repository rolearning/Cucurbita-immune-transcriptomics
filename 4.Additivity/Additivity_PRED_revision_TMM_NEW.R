#!/usr/bin/env Rscript

# Analyze additivity/PRED and produce Fig. 5 and Supplementary Fig. 12 from
# raw scaledTPM count-like matrices, cross-species sample factors, and the OG map.
# Normalize parental samples before calculating weighted MPV and DESeq2 contrasts.
# Use TMM by default or median_ratio via --method, retaining method-specific outputs.
# Export analysis/classification tables, figure panels, and per-OG candidate plots.

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(openxlsx)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

script_args <- commandArgs(trailingOnly = TRUE)
method_arg <- "TMM"
method_arg_index <- match("--method", script_args)
method_equals_index <- grep("^--method=", script_args)
if (length(method_equals_index)) {
  method_arg <- sub("^--method=", "", script_args[method_equals_index[1]])
} else if (!is.na(method_arg_index) && method_arg_index < length(script_args)) {
  method_arg <- script_args[method_arg_index + 1]
}
if (!method_arg %in% c("TMM", "median_ratio")) {
  stop("--method must be TMM or median_ratio; received: ", method_arg)
}

analysis_version <- paste0("Corrected_order_NEW_all_BS_", method_arg)
lfc_threshold <- 1
alpha_threshold <- 0.05
timepoints_expected <- c("00", "24", "48")

# Use factors estimated from all simple BS groups; the selected method is shared across species/timepoints.
normalization_reference_set <- "all_BS"
normalization_method <- method_arg

project_root <- file.path(RELEASE_ROOT, "4.Additivity")

input_files <- c(
  Cm_counts = file.path(RELEASE_ROOT, "1.normalization/results/01_tximport/Cm/Cm_scaledTPM_counts.tsv"),
  Ls_counts = file.path(RELEASE_ROOT, "1.normalization/results/01_tximport/Ls/Ls_scaledTPM_counts.tsv"),
  Ta_counts = file.path(RELEASE_ROOT, "1.normalization/results/01_tximport/Ta/Ta_scaledTPM_counts.tsv"),
  scaling = file.path(RELEASE_ROOT, "1.normalization/results/02_cross_species_scaling/all_candidate_sample_factors.tsv"),
  OG_info = file.path(RELEASE_ROOT, "OG_info/OG_info_260906_OK.xlsx")
)

resolve_path <- function(path, root = project_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(root, path)
}

input_paths <- vapply(input_files, resolve_path, character(1))
output_dir <- resolve_path(file.path("results", analysis_version))
table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "figures")
per_og_dir <- file.path(plot_dir, "5.perOG_figs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(per_og_dir, recursive = TRUE, showWarnings = FALSE)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"))
}

read_count_matrix <- function(path) {
  dat <- read_tsv(
    path,
    col_types = cols(.default = col_double(), GID = col_character()),
    show_col_types = FALSE
  )
  if (!"GID" %in% names(dat)) {
    stop("Expression matrix identifier column must be named GID: ", path)
  }
  if (anyDuplicated(dat$GID)) stop("Duplicated GID in expression matrix: ", path)
  mat <- as.matrix(column_to_rownames(as.data.frame(dat), "GID"))
  storage.mode(mat) <- "double"
  if (any(!is.finite(mat))) stop("Non-finite expression value in: ", path)
  if (any(mat < 0)) stop("Negative expression value in: ", path)
  mat
}

parse_sample_meta <- function(sample_names) {
  match_mat <- str_match(
    sample_names,
    "^(Cm|Ls|Ta)_(\\d{1,2})h_D(\\d+)$"
  )
  if (anyNA(match_mat[, 1])) {
    stop("Unexpected sample names: ",
         paste(sample_names[is.na(match_mat[, 1])], collapse = ", "))
  }
  tibble(
    sample = sample_names,
    species = match_mat[, 2],
    timepoint = sprintf("%02d", as.integer(match_mat[, 3])),
    replicate = paste0("D", match_mat[, 4])
  )
}

aggregate_to_og <- function(count_mat, gene_to_og) {
  mapped_og <- unname(gene_to_og[rownames(count_mat)])
  keep <- !is.na(mapped_og)
  if (!any(keep)) stop("No genes could be mapped to OGs.")

  dat <- as.data.frame(count_mat[keep, , drop = FALSE])
  dat$OG <- mapped_og[keep]
  aggregated <- dat %>%
    group_by(OG) %>%
    summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)),
              .groups = "drop")

  mat <- as.matrix(aggregated[, -1, drop = FALSE])
  rownames(mat) <- aggregated$OG
  storage.mode(mat) <- "double"
  mat
}

align_og_matrix <- function(mat, ogs) {
  out <- matrix(
    0,
    nrow = length(ogs),
    ncol = ncol(mat),
    dimnames = list(ogs, colnames(mat))
  )
  overlap <- intersect(rownames(mat), ogs)
  out[overlap, ] <- mat[overlap, , drop = FALSE]
  storage.mode(out) <- "double"
  out
}

build_parent_map <- function(og_info, subgenome_name) {
  og_info %>%
    filter(
      subgenome == subgenome_name,
      !is.na(GeneName_base), GeneName_base != "",
      !is.na(parent), parent != ""
    ) %>%
    select(OG = GeneName_base, parent) %>%
    mutate(parent = as.character(parent)) %>%
    distinct() %>%
    separate_rows(parent, sep = ",\\s*") %>%
    mutate(parent = trimws(parent)) %>%
    distinct(OG, parent) %>%
    transmute(id = parent, OG)
}

validate_unique_map <- function(map_table, label) {
  id_frequency <- table(map_table$id, useNA = "no")
  duplicated_ids <- names(id_frequency)[id_frequency > 1]
  if (length(duplicated_ids)) {
    stop(label, " IDs map to more than one OG: ",
         paste(head(duplicated_ids, 20), collapse = ", "))
  }
}

validate_map_ids_in_expression <- function(map_table, expression_mat, label) {
  missing_ids <- setdiff(unique(map_table$id), rownames(expression_mat))
  if (length(missing_ids)) {
    stop(
      label, " mapped gene IDs absent from its expression matrix (first 20): ",
      paste(head(missing_ids, 20), collapse = ", ")
    )
  }
}

get_weights <- function(ogs, wB, wS) {
  wb <- unname(wB[ogs])
  ws <- unname(wS[ogs])
  missing <- is.na(wb) | is.na(ws)
  if (any(missing)) {
    stop("Missing Cm copy-number weights for OG(s): ",
         paste(head(ogs[missing], 20), collapse = ", "))
  }
  if (any(abs((wb + ws) - 1) > 1e-12)) {
    stop("Invalid copy-number weights: wB + wS must equal 1.")
  }
  list(wB = wb, wS = ws)
}

align_parent_replicates <- function(B, S, factor_B, factor_S) {
  rep_B <- str_extract(colnames(B), "D\\d+$")
  rep_S <- str_extract(colnames(S), "D\\d+$")
  common <- intersect(rep_B, rep_S)
  if (!length(common)) stop("No matched parental replicates were found.")

  B2 <- B[, match(common, rep_B), drop = FALSE]
  S2 <- S[, match(common, rep_S), drop = FALSE]
  factor_B2 <- unname(factor_B[match(common, rep_B)])
  factor_S2 <- unname(factor_S[match(common, rep_S)])
  if (any(!is.finite(c(factor_B2, factor_S2))) ||
      any(c(factor_B2, factor_S2) <= 0)) {
    stop("Invalid parental scaling factors after replicate matching.")
  }

  list(
    B = B2,
    S = S2,
    factor_B = factor_B2,
    factor_S = factor_S2,
    replicate = common
  )
}

# Normalize each parent before weighted MPV calculation, then apply a synthetic size factor.
# Dividing by that same factor in DESeq2 recovers the normalized MPV.
make_mpv_weighted <- function(B, S, factor_B, factor_S, wB_row, wS_row) {
  aligned <- align_parent_replicates(B, S, factor_B, factor_S)

  scaled_B <- sweep(aligned$B, 2, aligned$factor_B, "/")
  scaled_S <- sweep(aligned$S, 2, aligned$factor_S, "/")
  weighted_B <- sweep(scaled_B, 1, wB_row, `*`)
  weighted_S <- sweep(scaled_S, 1, wS_row, `*`)
  mpv_normalized <- weighted_B + weighted_S

  mpv_factor <- (aligned$factor_B + aligned$factor_S) / 2
  mpv_counts <- sweep(mpv_normalized, 2, mpv_factor, "*")
  colnames(mpv_counts) <- paste0("MPVw_", aligned$replicate)
  names(mpv_factor) <- colnames(mpv_counts)
  list(counts = mpv_counts, factor = mpv_factor)
}

run_deseq2_threshold <- function(count_matrix, size_factors, groups, contrast,
                                 lfc_threshold = 1) {
  numerator <- contrast[1]
  denominator <- contrast[2]
  keep_samples <- groups %in% contrast
  counts2 <- count_matrix[, keep_samples, drop = FALSE]
  if (any(!is.finite(counts2)) || any(counts2 < 0)) {
    stop("Non-finite or negative count-like value for contrast ",
         numerator, " vs ", denominator, ".")
  }
  counts2 <- round(counts2)
  storage.mode(counts2) <- "integer"
  groups2 <- factor(groups[keep_samples], levels = c(denominator, numerator))

  coldata <- data.frame(
    group = groups2,
    row.names = colnames(counts2)
  )
  dds <- DESeqDataSetFromMatrix(counts2, coldata, design = ~ group)
  assigned_factors <- unname(size_factors[colnames(counts2)])
  if (any(!is.finite(assigned_factors)) || any(assigned_factors <= 0)) {
    stop("Missing or invalid size factors for contrast ",
         numerator, " vs ", denominator, ".")
  }
  sizeFactors(dds) <- assigned_factors
  dds <- DESeq(dds, quiet = TRUE)

  if (is.null(lfc_threshold) || lfc_threshold <= 0) {
    res_mle <- results(dds, contrast = c("group", numerator, denominator))
  } else {
    res_mle <- results(
      dds,
      contrast = c("group", numerator, denominator),
      lfcThreshold = lfc_threshold,
      altHypothesis = "greaterAbs"
    )
  }

  coef_name <- grep(
    paste0("^group_", numerator, "_vs_", denominator, "$"),
    resultsNames(dds),
    value = TRUE
  )
  if (length(coef_name) != 1) {
    stop("DESeq2 coefficient not found for ", numerator, " vs ", denominator)
  }
  res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res_mle)
  mle <- as.data.frame(res_mle)
  shr <- as.data.frame(res_shr)
  svalue <- if ("svalue" %in% names(shr)) shr$svalue else rep(NA_real_, nrow(shr))

  normalized <- counts(dds, normalized = TRUE)
  samples_A <- rownames(coldata)[coldata$group == numerator]
  samples_B <- rownames(coldata)[coldata$group == denominator]

  out <- tibble(
    OG = rownames(mle),
    baseMean = mle$baseMean,
    log2FC_MLE = mle$log2FoldChange,
    lfcSE = mle$lfcSE,
    CI95_low = mle$log2FoldChange - 1.96 * mle$lfcSE,
    CI95_high = mle$log2FoldChange + 1.96 * mle$lfcSE,
    pvalue = mle$pvalue,
    padj = mle$padj,
    lfc_shr = shr$log2FoldChange,
    svalue = svalue,
    minusLog10Padj = -log10(mle$padj + 1e-300),
    mean_A = rowMeans(normalized[, samples_A, drop = FALSE]),
    mean_B = rowMeans(normalized[, samples_B, drop = FALSE])
  )
  names(out)[names(out) == "mean_A"] <- paste0("mean_norm_", numerator)
  names(out)[names(out) == "mean_B"] <- paste0("mean_norm_", denominator)
  out
}

cm_raw <- read_count_matrix(input_paths[["Cm_counts"]])
ls_raw <- read_count_matrix(input_paths[["Ls_counts"]])
ta_raw <- read_count_matrix(input_paths[["Ta_counts"]])
og_info <- openxlsx::read.xlsx(input_paths[["OG_info"]], sheet = 1, detectDates = FALSE) %>%
  as.data.frame()
scaling_table <- read_tsv(input_paths[["scaling"]], show_col_types = FALSE)

required_factor_columns <- c(
  "reference_set", "method", "sample", "species", "multiplier"
)
if (!all(required_factor_columns %in% names(scaling_table))) {
  stop(
    "all_candidate_sample_factors.tsv is missing required columns: ",
    paste(setdiff(required_factor_columns, names(scaling_table)), collapse = ", ")
  )
}

required_og_columns <- c(
  "Cmo", "GeneName_base", "Module", "subgenome", "Type", "parent"
)
if (!all(required_og_columns %in% names(og_info))) {
  stop("OG_info is missing required columns: ",
       paste(setdiff(required_og_columns, names(og_info)), collapse = ", "))
}

cmB_map <- og_info %>%
  filter(
    subgenome == "CmB",
    !is.na(Cmo), Cmo != "",
    !is.na(GeneName_base), GeneName_base != ""
  ) %>%
  transmute(id = Cmo, OG = GeneName_base) %>%
  distinct()
cmS_map <- og_info %>%
  filter(
    subgenome == "CmS",
    !is.na(Cmo), Cmo != "",
    !is.na(GeneName_base), GeneName_base != ""
  ) %>%
  transmute(id = Cmo, OG = GeneName_base) %>%
  distinct()
ls_parent_map <- build_parent_map(og_info, "CmB")
ta_parent_map <- build_parent_map(og_info, "CmS")
validate_unique_map(cmB_map, "CmB")
validate_unique_map(cmS_map, "CmS")
validate_unique_map(ls_parent_map, "Ls")
validate_unique_map(ta_parent_map, "Ta")
validate_map_ids_in_expression(cmB_map, cm_raw, "CmB")
validate_map_ids_in_expression(cmS_map, cm_raw, "CmS")
validate_map_ids_in_expression(ls_parent_map, ls_raw, "Ls")
validate_map_ids_in_expression(ta_parent_map, ta_raw, "Ta")

cmB_gene_to_og <- setNames(cmB_map$OG, cmB_map$id)
cmS_gene_to_og <- setNames(cmS_map$OG, cmS_map$id)
ls_gene_to_og <- setNames(ls_parent_map$OG, ls_parent_map$id)
ta_gene_to_og <- setNames(ta_parent_map$OG, ta_parent_map$id)

copy_number <- og_info %>%
  filter(
    subgenome %in% c("CmB", "CmS"),
    !is.na(Cmo), Cmo != "",
    !is.na(GeneName_base), GeneName_base != ""
  ) %>%
  distinct(GeneName_base, subgenome, Cmo) %>%
  dplyr::count(GeneName_base, subgenome, name = "n") %>%
  pivot_wider(names_from = subgenome, values_from = n, values_fill = 0) %>%
  transmute(
    OG = GeneName_base,
    n_B = as.integer(CmB),
    n_S = as.integer(CmS),
    denominator = pmax(n_B + n_S, 1L),
    wB = n_B / denominator,
    wS = n_S / denominator
  )
wB <- setNames(copy_number$wB, copy_number$OG)
wS <- setNames(copy_number$wS, copy_number$OG)

cmB_og <- aggregate_to_og(cm_raw, cmB_gene_to_og)
cmS_og <- aggregate_to_og(cm_raw, cmS_gene_to_og)
ls_og <- aggregate_to_og(ls_raw, ls_gene_to_og)
ta_og <- aggregate_to_og(ta_raw, ta_gene_to_og)

all_ogs <- Reduce(
  union,
  list(rownames(cmB_og), rownames(cmS_og), rownames(ls_og), rownames(ta_og))
)
cmB_og <- align_og_matrix(cmB_og, all_ogs)
cmS_og <- align_og_matrix(cmS_og, all_ogs)
ls_og <- align_og_matrix(ls_og, all_ogs)
ta_og <- align_og_matrix(ta_og, all_ogs)
cm_total_og <- cmB_og + cmS_og
storage.mode(cm_total_og) <- "double"

# Each parental proxy is the mean of its unique mapped genes; retain decimals until DESeq2 input.
n_Ls <- table(ls_parent_map$OG)
n_Ta <- table(ta_parent_map$OG)
n_Ls <- as.numeric(n_Ls[all_ogs])
n_Ta <- as.numeric(n_Ta[all_ogs])
n_Ls[is.na(n_Ls)] <- 1
n_Ta[is.na(n_Ta)] <- 1
ls_og <- sweep(ls_og, 1, n_Ls, "/")
ta_og <- sweep(ta_og, 1, n_Ta, "/")
storage.mode(ls_og) <- "double"
storage.mode(ta_og) <- "double"

cm_meta <- parse_sample_meta(colnames(cm_total_og))
ls_meta <- parse_sample_meta(colnames(ls_og))
ta_meta <- parse_sample_meta(colnames(ta_og))

all_sample_meta <- bind_rows(cm_meta, ls_meta, ta_meta) %>%
  mutate(species = sub("_[BS]$", "", species))
if (anyDuplicated(all_sample_meta$sample)) {
  stop("A biological sample ID occurs in more than one input matrix.")
}

factor_selected <- scaling_table %>%
  filter(
    reference_set == normalization_reference_set,
    method == normalization_method
  ) %>%
  select(reference_set, method, sample, species, multiplier)

if (nrow(factor_selected) == 0) {
  stop(
    "No rows found for reference_set = ", normalization_reference_set,
    " and method = ", normalization_method, "."
  )
}
if (anyDuplicated(factor_selected$sample)) {
  stop("Selected factor table contains duplicated sample IDs.")
}
if (any(!is.finite(factor_selected$multiplier)) ||
    any(factor_selected$multiplier <= 0)) {
  stop("Selected multiplier values must be finite and positive.")
}

expected_samples <- all_sample_meta$sample
missing_factors <- setdiff(expected_samples, factor_selected$sample)
extra_factors <- setdiff(factor_selected$sample, expected_samples)
if (length(missing_factors) || length(extra_factors)) {
  stop(
    "Sample mismatch between expression matrices and selected factors. ",
    "Missing factors: ", paste(missing_factors, collapse = ", "),
    "; extra factor rows: ", paste(extra_factors, collapse = ", ")
  )
}

factor_selected <- factor_selected %>%
  left_join(
    all_sample_meta %>% select(sample, matrix_species = species),
    by = "sample"
  )
if (any(is.na(factor_selected$species)) ||
    any(factor_selected$species != factor_selected$matrix_species)) {
  bad <- factor_selected$sample[
    is.na(factor_selected$species) |
      factor_selected$species != factor_selected$matrix_species
  ]
  stop("Species mismatch in factor table for sample(s): ", paste(bad, collapse = ", "))
}

# Center factors jointly across species and timepoints to retain cross-species scale differences.
factor_center <- exp(mean(log(1 / factor_selected$multiplier)))
factor_selected <- factor_selected %>%
  mutate(
    raw_deseq_size_factor = 1 / multiplier,
    joint_geometric_center = .env$factor_center,
    deseq_size_factor = raw_deseq_size_factor / joint_geometric_center
  )

sf_all <- setNames(factor_selected$deseq_size_factor, factor_selected$sample)
sf_cm_new <- sf_all[cm_meta$sample]
sf_ls_new <- sf_all[ls_meta$sample]
sf_ta_new <- sf_all[ta_meta$sample]

write_tsv(
  factor_selected %>%
    select(
      reference_set, method, sample, species, multiplier,
      raw_deseq_size_factor, joint_geometric_center, deseq_size_factor
    ),
  file.path(table_dir, "selected_cross_species_sample_factors.tsv")
)

timepoints <- Reduce(
  intersect,
  list(unique(cm_meta$timepoint), unique(ls_meta$timepoint), unique(ta_meta$timepoint))
)
timepoints <- timepoints_expected[timepoints_expected %in% timepoints]
if (!identical(timepoints, timepoints_expected)) {
  stop("Expected timepoints 00, 24, and 48 in all count matrices.")
}

bs_ogs <- unique(og_info$GeneName_base[og_info$Type == "BS"])
res_additivity <- list()
res_pred_B <- list()
res_pred_S <- list()

for (tp in timepoints) {
  tp_h <- paste0(as.integer(tp), "h")
  cm_cols <- cm_meta %>% filter(timepoint == tp) %>% pull(sample)
  ls_cols <- ls_meta %>% filter(timepoint == tp) %>% pull(sample)
  ta_cols <- ta_meta %>% filter(timepoint == tp) %>% pull(sample)

  cm_reps <- str_extract(cm_cols, "D\\d+$")
  ls_reps <- str_extract(ls_cols, "D\\d+$")
  ta_reps <- str_extract(ta_cols, "D\\d+$")
  common_reps <- Reduce(intersect, list(cm_reps, ls_reps, ta_reps))
  if (!length(common_reps)) stop("No common replicates at timepoint ", tp)
  cm_cols <- cm_cols[match(common_reps, cm_reps)]
  ls_cols <- ls_cols[match(common_reps, ls_reps)]
  ta_cols <- ta_cols[match(common_reps, ta_reps)]

  Cm <- cm_total_og[, cm_cols, drop = FALSE]
  ELs <- ls_og[, ls_cols, drop = FALSE]
  ETa <- ta_og[, ta_cols, drop = FALSE]
  weights <- get_weights(rownames(Cm), wB, wS)
  mpv <- make_mpv_weighted(
    ELs,
    ETa,
    sf_ls_new[ls_cols],
    sf_ta_new[ta_cols],
    weights$wB,
    weights$wS
  )

  colnames(Cm) <- paste0("Cm_", tp_h, "_", common_reps)
  colnames(ELs) <- paste0("ELs_", tp_h, "_", common_reps)
  colnames(ETa) <- paste0("ETa_", tp_h, "_", common_reps)
  colnames(mpv$counts) <- paste0("MPVw_", tp_h, "_", common_reps)
  names(mpv$factor) <- colnames(mpv$counts)

  add_counts <- round(cbind(Cm, mpv$counts, ELs, ETa))
  add_counts <- add_counts[rownames(add_counts) %in% bs_ogs, , drop = FALSE]
  add_counts <- add_counts[rowSums(add_counts > 0) >= 2, , drop = FALSE]
  add_factors <- c(
    sf_cm_new[colnames(Cm)],
    mpv$factor,
    sf_ls_new[sub("^ELs_", "Ls_", colnames(ELs))],
    sf_ta_new[sub("^ETa_", "Ta_", colnames(ETa))]
  )
  names(add_factors) <- colnames(add_counts)
  add_groups <- case_when(
    str_starts(colnames(add_counts), "Cm_") ~ "Cm",
    str_starts(colnames(add_counts), "MPVw_") ~ "MPVw",
    str_starts(colnames(add_counts), "ELs_") ~ "ELs",
    str_starts(colnames(add_counts), "ETa_") ~ "ETa",
    TRUE ~ NA_character_
  )

  res_additivity[[tp]] <- list(
    Cm_vs_MPVw = run_deseq2_threshold(
      add_counts, add_factors, add_groups, c("Cm", "MPVw"), lfc_threshold
    ),
    Cm_vs_ELs = run_deseq2_threshold(
      add_counts, add_factors, add_groups, c("Cm", "ELs"), lfc_threshold
    ),
    Cm_vs_ETa = run_deseq2_threshold(
      add_counts, add_factors, add_groups, c("Cm", "ETa"), lfc_threshold
    )
  )

  CmB <- cmB_og[, cm_cols, drop = FALSE]
  CmS <- cmS_og[, cm_cols, drop = FALSE]
  colnames(CmB) <- paste0("CmB_", tp_h, "_", common_reps)
  colnames(CmS) <- paste0("CmS_", tp_h, "_", common_reps)
  colnames(ELs) <- paste0("Ls_", tp_h, "_", common_reps)
  colnames(ETa) <- paste0("Ta_", tp_h, "_", common_reps)

  pred_weights <- get_weights(rownames(ELs), wB, wS)
  weighted_ELs <- sweep(ELs, 1, pred_weights$wB, `*`)
  weighted_ETa <- sweep(ETa, 1, pred_weights$wS, `*`)
  pred_B_counts <- round(cbind(CmB, weighted_ELs))
  pred_S_counts <- round(cbind(CmS, weighted_ETa))

  pred_B_factors <- c(
    sf_cm_new[sub("^CmB_", "Cm_", colnames(CmB))],
    sf_ls_new[colnames(weighted_ELs)]
  )
  pred_S_factors <- c(
    sf_cm_new[sub("^CmS_", "Cm_", colnames(CmS))],
    sf_ta_new[colnames(weighted_ETa)]
  )
  names(pred_B_factors) <- colnames(pred_B_counts)
  names(pred_S_factors) <- colnames(pred_S_counts)

  pred_B_counts <- pred_B_counts[rowSums(pred_B_counts > 0) >= 2, , drop = FALSE]
  pred_S_counts <- pred_S_counts[rowSums(pred_S_counts > 0) >= 2, , drop = FALSE]
  pred_B_groups <- ifelse(str_starts(colnames(pred_B_counts), "CmB_"), "CmB", "Ls")
  pred_S_groups <- ifelse(str_starts(colnames(pred_S_counts), "CmS_"), "CmS", "Ta")

  res_pred_B[[tp]] <- run_deseq2_threshold(
    pred_B_counts,
    pred_B_factors,
    pred_B_groups,
    c("CmB", "Ls"),
    lfc_threshold
  )
  res_pred_S[[tp]] <- run_deseq2_threshold(
    pred_S_counts,
    pred_S_factors,
    pred_S_groups,
    c("CmS", "Ta"),
    lfc_threshold
  )
}

results <- list(
  additivity = res_additivity,
  pred = list(CmB_vs_Ls = res_pred_B, CmS_vs_Ta = res_pred_S)
)

standardize_result <- function(df, cm_component, proxy, timepoint,
                               cm_mean_column, proxy_mean_column) {
  df %>%
    transmute(
      OG,
      baseMean,
      log2FC_MLE,
      lfcSE,
      CI95_low,
      CI95_high,
      pvalue,
      padj,
      lfc_shr,
      svalue,
      minusLog10Padj,
      mean_Cm_range = .data[[cm_mean_column]],
      mean_Compare = .data[[proxy_mean_column]],
      Cm_range = cm_component,
      Compare = proxy,
      timepoint = timepoint
    )
}

result_parts <- list()
for (tp in timepoints) {
  result_parts[[paste(tp, "MPVw")]] <- standardize_result(
    results$additivity[[tp]]$Cm_vs_MPVw,
    "Cm_total", "MPVw", tp, "mean_norm_Cm", "mean_norm_MPVw"
  )
  result_parts[[paste(tp, "ELs")]] <- standardize_result(
    results$additivity[[tp]]$Cm_vs_ELs,
    "Cm_total", "ELs", tp, "mean_norm_Cm", "mean_norm_ELs"
  )
  result_parts[[paste(tp, "ETa")]] <- standardize_result(
    results$additivity[[tp]]$Cm_vs_ETa,
    "Cm_total", "ETa", tp, "mean_norm_Cm", "mean_norm_ETa"
  )
  result_parts[[paste(tp, "PRED_B")]] <- standardize_result(
    results$pred$CmB_vs_Ls[[tp]],
    "CmB", "ELs", tp, "mean_norm_CmB", "mean_norm_Ls"
  )
  result_parts[[paste(tp, "PRED_S")]] <- standardize_result(
    results$pred$CmS_vs_Ta[[tp]],
    "CmS", "ETa", tp, "mean_norm_CmS", "mean_norm_Ta"
  )
}

og_annotation <- og_info %>%
  transmute(OG = GeneName_base, Module, Type) %>%
  distinct()
result_internal <- bind_rows(result_parts) %>%
  left_join(og_annotation, by = "OG") %>%
  left_join(copy_number %>% select(OG, wB, wS), by = "OG")

is_sig <- function(p, alpha) !is.na(p) && p <= alpha
is_nonsig <- function(p, alpha) !is.na(p) && p > alpha
is_greater <- function(a, b, tolerance) {
  !is.na(a) && !is.na(b) && a > b + tolerance
}
is_less <- function(a, b, tolerance) {
  !is.na(a) && !is.na(b) && a < b - tolerance
}

classify_one <- function(row, alpha = 0.05, tolerance = 1e-9) {
  lfc_mpv <- row$lfc_shr_MPVw
  p_mpv <- row$padj_MPVw
  p_els <- row$padj_ELs
  p_eta <- row$padj_ETa
  Cm <- row$mean_Cm_total
  ELs <- row$mean_ELs
  ETa <- row$mean_ETa

  if (is.na(p_mpv)) return("ambiguous")

  # "Additive" denotes no detected deviation beyond the tested threshold.
  if (is_nonsig(p_mpv, alpha)) return("additive")

  greater_Ls <- is_greater(Cm, ELs, tolerance)
  less_Ls <- is_less(Cm, ELs, tolerance)
  greater_Ta <- is_greater(Cm, ETa, tolerance)
  less_Ta <- is_less(Cm, ETa, tolerance)
  between_parents <- (greater_Ls && less_Ta) || (less_Ls && greater_Ta)

  if (greater_Ls && greater_Ta && is_sig(p_els, alpha) &&
      is_sig(p_eta, alpha) && lfc_mpv > 1) return("transgressive_up")
  if (less_Ls && less_Ta && is_sig(p_els, alpha) &&
      is_sig(p_eta, alpha) && lfc_mpv < -1) return("transgressive_down")

  if (is_nonsig(p_els, alpha) && is_sig(p_eta, alpha)) {
    if (!is.na(lfc_mpv) && lfc_mpv > 1) return("bias_toward_Ls_up")
    if (!is.na(lfc_mpv) && lfc_mpv < -1) return("bias_toward_Ls_down")
    return("bias_toward_Ls")
  }
  if (is_nonsig(p_eta, alpha) && is_sig(p_els, alpha)) {
    if (!is.na(lfc_mpv) && lfc_mpv > 1) return("bias_toward_Ta_up")
    if (!is.na(lfc_mpv) && lfc_mpv < -1) return("bias_toward_Ta_down")
    return("bias_toward_Ta")
  }
  if (between_parents && is_sig(p_els, alpha) && is_sig(p_eta, alpha)) {
    return("incomplete_dominance")
  }
  if (is_sig(p_mpv, alpha) && (is.na(p_els) || is.na(p_eta))) {
    return("nonadditive_unclear")
  }
  "ambiguous"
}

classification <- result_internal %>%
  filter(Cm_range == "Cm_total", Compare %in% c("MPVw", "ELs", "ETa")) %>%
  select(OG, timepoint, Compare, lfc_shr, padj, mean_Cm_range, mean_Compare) %>%
  mutate(Compare = factor(Compare, levels = c("MPVw", "ELs", "ETa"))) %>%
  pivot_wider(
    names_from = Compare,
    values_from = c(lfc_shr, padj, mean_Compare),
    names_sep = "_"
  ) %>%
  rename(
    mean_Cm_total = mean_Cm_range,
    mean_ELs = mean_Compare_ELs,
    mean_ETa = mean_Compare_ETa
  ) %>%
  rowwise() %>%
  mutate(category = classify_one(pick(everything()), alpha_threshold)) %>%
  ungroup() %>%
  left_join(og_annotation, by = "OG") %>%
  left_join(copy_number %>% select(OG, wB, wS), by = "OG")

supp_table_s9 <- result_internal %>%
  transmute(
    OG,
    Cm_expression_component = Cm_range,
    Progenitor_proxy = Compare,
    Timepoint = as.integer(timepoint),
    Module,
    Type,
    wB,
    wS,
    mean_Cm_group = mean_Cm_range,
    mean_Compare,
    log2FC_MLE,
    lfcSE,
    pvalue,
    padj
  )
write.csv(
  supp_table_s9,
  file.path(table_dir, "result_tab_full.csv"),
  row.names = FALSE
)

supp_table_s10 <- classification %>%
  select(OG, Module, timepoint, category) %>%
  mutate(timepoint = paste0("category_", as.integer(timepoint))) %>%
  pivot_wider(names_from = timepoint, values_from = category) %>%
  select(OG, Module, category_0, category_24, category_48)
write.csv(
  supp_table_s10,
  file.path(table_dir, "additivity_catogory_all_timepoints_padj.csv"),
  row.names = FALSE
)

category_colors <- c(
  additive = "#C5C5C5",
  bias_toward_Ls = "#5A9BD4",
  bias_toward_Ls_down = "#5A9BD4",
  bias_toward_Ls_up = "#5A9BD4",
  bias_toward_Ta = "#C0392B",
  bias_toward_Ta_down = "#C0392B",
  bias_toward_Ta_up = "#C0392B",
  incomplete_dominance = "#58B793",
  transgressive_down = "#D98841",
  transgressive_up = "#D98841",
  nonadditive_unclear = "grey40",
  ambiguous = "grey70"
)

heatmap_category_colors <- c(
  additive = "#C5C5C5",
  bias_toward_Ls = "#5A9BD4",
  bias_toward_Ls_down = "#9ee5f0",
  bias_toward_Ls_up = "#5A9BD4",
  bias_toward_Ta = "#C0392B",
  bias_toward_Ta_down = "#E6A5A5",
  bias_toward_Ta_up = "#C0392B",
  incomplete_dominance = "#58B793",
  transgressive_down = "#f9d55d",
  transgressive_up = "#D98841",
  nonadditive_unclear = "grey40",
  ambiguous = "grey70"
)
bind_additivity_timepoint <- function(result_list) {
  imap_dfr(
    list(
      "Cm vs MPVw" = result_list$Cm_vs_MPVw,
      "Cm vs ELs" = result_list$Cm_vs_ELs,
      "Cm vs ETa" = result_list$Cm_vs_ETa
    ),
    ~ .x %>% select(OG, log2FC_MLE, padj, minusLog10Padj) %>%
      mutate(contrast = .y)
  )
}

plot_additivity_density <- function(result_list, x_limits) {
  dat <- bind_additivity_timepoint(result_list) %>%
    mutate(contrast = factor(
      contrast,
      levels = c("Cm vs MPVw", "Cm vs ELs", "Cm vs ETa")
    )) %>%
    arrange(contrast)
  colors <- c(
    "Cm vs MPVw" = "#66b2a0",
    "Cm vs ELs" = "#1e88e5",
    "Cm vs ETa" = "#ee6055"
  )
  alpha_values <- c("Cm vs MPVw" = 0.9, "Cm vs ELs" = 0.4, "Cm vs ETa" = 0.4)

  ggplot(dat, aes(x = log2FC_MLE, fill = contrast, alpha = contrast)) +
    geom_density(color = NA, adjust = 1, position = "identity") +
    scale_alpha_manual(values = alpha_values, guide = "none") +
    scale_fill_manual(values = alpha(colors, 0.25), guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(
      limits = x_limits,
      position = "top",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) +
    labs(x = "log2FC", y = NULL)
}

plot_pred_dual <- function(B_result, S_result, category_table, tp, x_limits) {
  dat <- inner_join(
    B_result %>% select(OG, LFC_B = log2FC_MLE),
    S_result %>% select(OG, LFC_S = log2FC_MLE),
    by = "OG"
  ) %>%
    left_join(
      category_table %>% filter(timepoint == tp) %>% select(OG, category),
      by = "OG"
    )

  ggplot(dat, aes(LFC_B, LFC_S, color = category)) +
    geom_point(alpha = 0.9, size = 3) +
    scale_color_manual(values = category_colors, drop = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_cartesian(xlim = x_limits) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15)
    ) +
    labs(
      x = "log2FC (CmB vs B proxy)",
      y = "log2FC (CmS vs S proxy)",
      title = ifelse(tp == "00", "0 hpi", paste0(tp, " hpi"))
    )
}

common_panel_xlim <- function(results, tp) {
  values <- c(
    bind_additivity_timepoint(results$additivity[[tp]])$log2FC_MLE,
    results$pred$CmB_vs_Ls[[tp]]$log2FC_MLE
  )
  values <- values[is.finite(values)]
  limits <- quantile(values, c(0.005, 0.995), na.rm = TRUE)
  retained <- values[values >= limits[1] & values <= limits[2]]
  span <- diff(range(retained))
  if (!is.finite(span) || span == 0) span <- 1
  range(retained) + c(-0.02, 0.02) * span
}

extract_additivity_per_og <- function(results_additivity, og) {
  bind_rows(lapply(timepoints_expected, function(tp) {
    if (!tp %in% names(results_additivity)) return(NULL)
    bind_rows(
      results_additivity[[tp]]$Cm_vs_MPVw %>%
        filter(OG == og) %>% mutate(contrast = "Cm vs MPVw"),
      results_additivity[[tp]]$Cm_vs_ELs %>%
        filter(OG == og) %>% mutate(contrast = "Cm vs ELs"),
      results_additivity[[tp]]$Cm_vs_ETa %>%
        filter(OG == og) %>% mutate(contrast = "Cm vs ETa")
    ) %>% mutate(timepoint = tp)
  }))
}

extract_pred_per_og <- function(results_pred, og) {
  bind_rows(lapply(timepoints_expected, function(tp) {
    if (!tp %in% names(results_pred$CmB_vs_Ls) ||
        !tp %in% names(results_pred$CmS_vs_Ta)) return(NULL)
    inner_join(
      results_pred$CmB_vs_Ls[[tp]] %>%
        filter(OG == og) %>% select(OG, LFC_B = log2FC_MLE),
      results_pred$CmS_vs_Ta[[tp]] %>%
        filter(OG == og) %>% select(OG, LFC_S = log2FC_MLE),
      by = "OG"
    ) %>% mutate(timepoint = tp)
  }))
}

plot_per_og <- function(results, og, output_path) {
  add_dat <- extract_additivity_per_og(results$additivity, og)
  pred_dat <- extract_pred_per_og(results$pred, og)
  all_values <- c(add_dat$log2FC_MLE, pred_dat$LFC_B, pred_dat$LFC_S)
  all_values <- all_values[is.finite(all_values)]
  if (!length(all_values)) return(invisible(NULL))
  axis_max <- max(2.5, max(abs(all_values)) * 1.04)
  common_limits <- c(-axis_max, axis_max)

  add_dat <- add_dat %>%
    mutate(
      timepoint_label = factor(
        ifelse(timepoint == "00", "0 hpi", paste0(timepoint, " hpi")),
        levels = c("48 hpi", "24 hpi", "0 hpi")
      )
    )
  add_colors <- c(
    "Cm vs MPVw" = "#f1a93b",
    "Cm vs ELs" = "#1e88e5",
    "Cm vs ETa" = "#ee6055"
  )
  top <- ggplot(add_dat, aes(log2FC_MLE, timepoint_label, color = contrast)) +
    geom_path(aes(group = contrast), alpha = 0.5, linewidth = 0.6) +
    geom_point(
      size = 3,
      alpha = 0.8,
      position = position_jitter(height = 0.05, width = 0, seed = 1)
    ) +
    scale_color_manual(values = add_colors, name = NULL) +
    scale_x_continuous(position = "top") +
    coord_cartesian(xlim = common_limits) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    labs(x = "log2FC (additivity)", y = NULL, tag = og) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      plot.tag = element_text(face = "bold", size = 11, hjust = 0),
      plot.tag.position = "top"
    )

  pred_dat <- pred_dat %>%
    mutate(timepoint = factor(timepoint, levels = timepoints_expected))
  time_colors <- c("00" = "#8ef0b4", "24" = "#4abf8a", "48" = "#2c6d4f")
  bottom <- ggplot(pred_dat, aes(LFC_B, LFC_S)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_path(aes(group = 1), linewidth = 0.6, alpha = 0.6, color = "grey60") +
    geom_point(aes(color = timepoint), size = 4, alpha = 0.95) +
    scale_color_manual(values = time_colors, name = "Timepoint") +
    coord_cartesian(xlim = common_limits, ylim = common_limits) +
    labs(
      x = "log2FC (CmB vs B proxy)",
      y = "log2FC (CmS vs S proxy)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14)
    )

  panel <- top / bottom + plot_layout(heights = c(1, 2))
  ggsave(output_path, panel, width = 5, height = 4.5)
  invisible(panel)
}

plot_data <- result_internal %>%
  mutate(
    Analysis = ifelse(
      Cm_range == "Cm_total",
      "additive",
      ifelse(Type == "BS", "PRED_BS", "PRED_Single")
    ),
    Source = ifelse(Type == "BS", "BS", "Single"),
    timepoint = sprintf("%02d", as.integer(timepoint))
  )

timepoint_colors <- c("00" = "#7dd8c3", "24" = "#ffb84d", "48" = "#9b79b9")
main_5a_left <- ggplot(
  filter(plot_data, Compare == "MPVw"),
  aes(log2FC_MLE, fill = timepoint)
) +
  geom_density(color = NA, adjust = 1, alpha = 0.7, position = "identity") +
  scale_fill_manual(
    values = timepoint_colors,
    breaks = c("00", "24", "48"),
    labels = c("0 hpi", "24 hpi", "48 hpi"),
    name = NULL
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-3, 5)) +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13)) +
  labs(x = "log2FC", y = NULL)
ggsave(
  file.path(plot_dir, "13.ND_1D_compare_BS_Single_FINAL-1.pdf"),
  main_5a_left,
  width = 4,
  height = 3.5
)

source_colors <- c("BS" = "#9b79b9", "Single" = "#7dd8c3")
plot_data2 <- plot_data %>%
  filter(Cm_range != "Cm_total") %>%
  mutate(Group = "All") %>%
  rbind(plot_data %>% filter(Cm_range == "CmB") %>% mutate(Group = "CmB")) %>%
  rbind(plot_data %>% filter(Cm_range == "CmS") %>% mutate(Group = "CmS")) %>%
  filter(Group != "All")

main_5a_right <- ggplot(
  filter(plot_data2, Cm_range != "Cm_total"),
  aes(log2FC_MLE, fill = Source)
) +
  geom_density(color = NA, adjust = 1, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = source_colors, name = NULL) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(Group ~ timepoint, axes = "all_x") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  labs(x = "log2FC", y = NULL, title = NULL)
ggsave(
  file.path(plot_dir, "13.ND_1D_compare_BS_Single_FINAL-2.pdf"),
  main_5a_right,
  width = 7.5,
  height = 4
)

panel_columns <- lapply(timepoints, function(tp) {
  x_limits <- common_panel_xlim(results, tp)
  top <- plot_additivity_density(results$additivity[[tp]], x_limits)
  bottom <- plot_pred_dual(
    results$pred$CmB_vs_Ls[[tp]],
    results$pred$CmS_vs_Ta[[tp]],
    classification,
    tp,
    x_limits
  )
  top / bottom + plot_layout(heights = c(1, 2))
})
main_5b_upper <- wrap_plots(panel_columns, nrow = 1) +
  plot_annotation(title = "Additivity LFC (top) and dual PRED (bottom)")
ggsave(
  file.path(plot_dir, "4.AddLFC_vs_NDdual_2x3_category_FINAL.pdf"),
  main_5b_upper,
  width = 10,
  height = 5.5
)

quadrant_data <- bind_rows(lapply(timepoints, function(tp) {
  inner_join(
    results$pred$CmB_vs_Ls[[tp]] %>%
      filter(OG %in% bs_ogs) %>%
      select(
        OG,
        LFC_B = log2FC_MLE,
        padj_B = padj,
        CI_low_B = CI95_low,
        CI_high_B = CI95_high
      ),
    results$pred$CmS_vs_Ta[[tp]] %>%
      filter(OG %in% bs_ogs) %>%
      select(
        OG,
        LFC_S = log2FC_MLE,
        padj_S = padj,
        CI_low_S = CI95_low,
        CI_high_S = CI95_high
      ),
    by = "OG"
  ) %>% mutate(timepoint = tp)
})) %>%
  mutate(
    Quadrant = case_when(
      LFC_B > 0 & LFC_S > 0 ~ "Q1",
      LFC_B > 0 & LFC_S < 0 ~ "Q2",
      LFC_B < 0 & LFC_S < 0 ~ "Q3",
      TRUE ~ "Q4"
    ),
    Quadrant = factor(Quadrant, levels = c("Q1", "Q2", "Q3", "Q4")),
    timepoint = factor(paste0(as.integer(timepoint), " hpi")),
    sig_padj = if_else(padj_B < 0.05 | padj_S < 0.05, "sig", "ns"),
    sig_ci = if_else(
      CI_low_B > 1 | CI_high_B < -1 | CI_low_S > 1 | CI_high_S < -1,
      "sig",
      "ns"
    )
  )

polar_all <- quadrant_data %>%
  dplyr::count(timepoint, Quadrant) %>%
  ggplot(aes(Quadrant, n, fill = Quadrant)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  coord_polar() +
  facet_wrap(~ timepoint) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "All genes: Quadrant distribution", y = "Count", x = NULL) +
  theme_minimal(base_size = 14)

make_significance_polar <- function(data, significance_column, title) {
  data %>%
    dplyr::count(timepoint, Quadrant, significance = .data[[significance_column]]) %>%
    ggplot(aes(Quadrant, n, fill = significance)) +
    geom_col(position = "stack", width = 1, color = "black", linewidth = 0.3) +
    coord_polar() +
    facet_wrap(~ timepoint) +
    scale_fill_manual(values = c(sig = "#7dd8c3", ns = "grey85")) +
    labs(title = title, y = "Count", x = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      panel.spacing = grid::unit(0.8, "cm")
    )
}

polar_ci <- make_significance_polar(
  quadrant_data, "sig_ci", "Quadrant distribution by CI ±1 significance"
)
polar_padj <- make_significance_polar(
  quadrant_data, "sig_padj", "Quadrant distribution by padj significance"
)
main_5b_lower <- polar_all / polar_ci / polar_padj
ggsave(
  file.path(plot_dir, "Quadrant_distribution_PolarBar_FINAL.pdf"),
  main_5b_lower,
  width = 9,
  height = 9
)

nonadditive_ogs <- classification %>%
  filter(!category %in% c("additive", "ambiguous"), !is.na(category)) %>%
  distinct(OG)
heat_data <- result_internal %>%
  filter(
    Cm_range == "Cm_total",
    Compare == "MPVw",
    OG %in% nonadditive_ogs$OG
  ) %>%
  left_join(classification %>% select(OG, timepoint, category),
            by = c("OG", "timepoint")) %>%
  mutate(
    timepoint = factor(as.integer(timepoint), levels = c(0, 24, 48)),
    OG = as.character(OG)
  )
og_levels <- heat_data %>% distinct(Module, OG) %>% arrange(Module, OG) %>% pull(OG)
heat_data <- heat_data %>% mutate(OG = factor(OG, levels = unique(og_levels)))

heat_long <- heat_data %>%
  select(OG, Module, timepoint, mean_Cm_range, mean_Compare) %>%
  pivot_longer(
    c(mean_Cm_range, mean_Compare),
    names_to = "Group",
    values_to = "value"
  ) %>%
  mutate(
    Group = recode(Group, mean_Cm_range = "Cm", mean_Compare = "MPVw"),
    Group = factor(Group, levels = c("Cm", "MPVw"))
  ) %>%
  group_by(OG) %>%
  mutate(value_scaled = as.numeric(scale(value))) %>%
  ungroup()

heatmap_plot <- ggplot(heat_long, aes(Group, OG, fill = value_scaled)) +
  scale_fill_gradient2(
    low = "#fee5d9", mid = "red", high = "#b01404", midpoint = 0,
    name = "z(mean)"
  ) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(
    rows = vars(Module), cols = vars(timepoint),
    scales = "free_y", space = "free_y", switch = "y"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text.y.left = element_text(angle = 0, size = 10),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

dot_plot <- ggplot(heat_data, aes(timepoint, OG, color = category)) +
  geom_point(size = 4.5, alpha = 0.95) +
  scale_color_manual(values = heatmap_category_colors) +
  scale_x_discrete(position = "top") +
  facet_grid(rows = vars(Module), scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y.left = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 11),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )
main_5c <- (heatmap_plot | dot_plot) + plot_layout(widths = c(3, 1.5))
ggsave(
  file.path(plot_dir, "18.Heatmap_for_non_add.pdf"),
  main_5c,
  width = 5.5,
  height = 3.5
)

walk(bs_ogs, function(og) {
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", og)
  plot_per_og(
    results,
    og,
    file.path(per_og_dir, paste0("OG_", safe_name, "_add_vs_nd.pdf"))
  )
})

pred_time_data <- result_internal %>%
  filter(Cm_range != "Cm_total", Type == "BS") %>%
  select(OG, timepoint, Module, Type, Cm_range, log2FC_MLE) %>%
  mutate(
    component = ifelse(Cm_range == "CmB", "B", "S"),
    timepoint = as.integer(timepoint)
  ) %>%
  select(-Cm_range) %>%
  pivot_wider(names_from = component, values_from = log2FC_MLE,
              names_prefix = "log2FC_") %>%
  pivot_wider(
    names_from = timepoint,
    values_from = c(log2FC_B, log2FC_S),
    names_glue = "{.value}_{timepoint}"
  ) %>%
  mutate(
    diff_log2FC_B_24_0 = log2FC_B_24 - log2FC_B_0,
    diff_log2FC_B_48_24 = log2FC_B_48 - log2FC_B_24,
    diff_log2FC_B_48_0 = log2FC_B_48 - log2FC_B_0,
    diff_log2FC_S_24_0 = log2FC_S_24 - log2FC_S_0,
    diff_log2FC_S_48_24 = log2FC_S_48 - log2FC_S_24,
    diff_log2FC_S_48_0 = log2FC_S_48 - log2FC_S_0
  ) %>%
  filter(Module %in% c("MAPKs", "RLCKs", "WRKYs")) %>%
  rowwise() %>%
  mutate(
    B_max_abs = max(abs(c(
      diff_log2FC_B_24_0, diff_log2FC_B_48_24, diff_log2FC_B_48_0
    )), na.rm = TRUE),
    S_max_abs = max(abs(c(
      diff_log2FC_S_24_0, diff_log2FC_S_48_24, diff_log2FC_S_48_0
    )), na.rm = TRUE),
    B_flag = B_max_abs > 1,
    S_flag = S_max_abs > 1,
    state = case_when(
      !B_flag & !S_flag ~ "FALSE/FALSE",
      B_flag & !S_flag ~ "TRUE/FALSE",
      !B_flag & S_flag ~ "FALSE/TRUE",
      TRUE ~ "TRUE/TRUE"
    )
  ) %>%
  ungroup()

state_levels <- c("FALSE/FALSE", "TRUE/FALSE", "FALSE/TRUE", "TRUE/TRUE")
supp_figure_12_data <- pred_time_data %>%
  dplyr::count(Module, state) %>%
  group_by(Module) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    Module = factor(Module, levels = c("MAPKs", "RLCKs", "WRKYs")),
    state = factor(state, levels = state_levels)
  )

supp_figure_12 <- ggplot(
  supp_figure_12_data,
  aes(Module, proportion, fill = state)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_discrete(labels = c(
    "FALSE/FALSE" = "None",
    "TRUE/FALSE" = "B",
    "FALSE/TRUE" = "S",
    "TRUE/TRUE" = "Both"
  )) +
  labs(
    x = NULL,
    y = "Proportion of BS pairs",
    fill = "B/S deviation state",
    title = "Homoeolog deviation states"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )
ggsave(
  file.path(plot_dir, "stacked_bar_100percent.pdf"),
  supp_figure_12,
  width = 5,
  height = 5
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message(
  "Normalization: ", normalization_reference_set, " + ", normalization_method,
  "; DESeq2 size factor = centered(1 / multiplier)."
)
message("Completed ", analysis_version, " analysis: ", output_dir)
