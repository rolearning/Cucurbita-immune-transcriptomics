#!/usr/bin/env Rscript

# Extract syngroup expression and generate the Fig. 2a cluster heatmaps.
# Read corrected syngroup maps and normalized scaledTPM matrices for 25 libraries.
# Sum Cm copy pairs before time-point averaging; retain TMM and median branches.
# Write syngroup matrices, k-means assignments, heatmaps, and sensitivity audits.

.run_release_stage <- function(code) {
  prior_wd <- getwd()
  on.exit(setwd(prior_wd), add = TRUE)
  eval(code, envir = new.env(parent = globalenv()))
}

.run_release_stage(quote({
RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)
options(stringsAsFactors = FALSE)

project_dir <- RELEASE_ROOT
output_dir <- file.path(project_dir, "2.global_expression/Fig2A_GO/kmeans")
syngroup_dir <- file.path(project_dir, "1.normalization")
normalization_dir <- file.path(
  project_dir,
  "1.normalization/results/02_cross_species_scaling"
)
tximport_dir <- file.path(project_dir, "1.normalization/results/01_tximport")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "TMM"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "median"), recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# Exclude Cm_48h_D2 and Ls_24h_D1 from the final 25-library set.
retained_samples <- list(
  Cm = c(
    "Cm_0h_D1", "Cm_0h_D2", "Cm_0h_D3",
    "Cm_24h_D1", "Cm_24h_D2", "Cm_24h_D3",
    "Cm_48h_D1", "Cm_48h_D3"
  ),
  Ls = c(
    "Ls_0h_D1", "Ls_0h_D2", "Ls_0h_D3",
    "Ls_24h_D2", "Ls_24h_D3",
    "Ls_48h_D1", "Ls_48h_D2", "Ls_48h_D3"
  ),
  Ta = c(
    "Ta_0h_D1", "Ta_0h_D2", "Ta_0h_D3",
    "Ta_24h_D1", "Ta_24h_D2", "Ta_24h_D3",
    "Ta_48h_D1", "Ta_48h_D2", "Ta_48h_D3"
  )
)
all_samples <- unlist(retained_samples, use.names = FALSE)
stopifnot(length(all_samples) == 25L, !anyDuplicated(all_samples))

library(dplyr)

pairs_file <- file.path(syngroup_dir, "All_pairs_syngroup_w_notes_gene_corrected.csv")
singles_file <- file.path(syngroup_dir, "All_single_syngroup_w_notes_gene_corrected.csv")
if (!file.exists(pairs_file) || !file.exists(singles_file)) {
  stop("Final syngroup definition files were not found.")
}

pairs_df <- read.csv(pairs_file, check.names = FALSE, stringsAsFactors = FALSE)
singles_df <- read.csv(singles_file, check.names = FALSE, stringsAsFactors = FALSE)
simple_pairs <- pairs_df %>% dplyr::filter(trimws(complexity) == "simple")
simple_singles <- singles_df %>% dplyr::filter(trimws(complexity) == "simple")

if (!all(c("Cm1", "Cm2", "Lsi", "Tan", "comp_id", "Type") %in% colnames(simple_pairs))) {
  stop("The pair syngroup definitions lack required columns.")
}
if (!all(c("Cm", "Lsi", "Tan", "comp_id") %in% colnames(simple_singles))) {
  stop("The single syngroup definitions lack required columns.")
}

# Reject fragment IDs superseded by the curated split-gene corrections.
legacy_gene_ids <- c(
  "CmoCh13G000470", "Lsi09G017820", "Lsi09G017830",
  "Lsi02G019240", "Lsi05G000620", "CmoCh19G004650",
  "Lsi06G003480", "CmoCh04G011740", "CmoCh19G008720"
)
assert_no_legacy_ids <- function(df, id_columns, label) {
  observed <- unique(unlist(df[id_columns], use.names = FALSE))
  observed <- trimws(observed[!is.na(observed) & observed != ""])
  residual <- intersect(legacy_gene_ids, observed)
  if (length(residual) > 0L) {
    stop(label, " still contains uncorrected gene IDs: ", paste(residual, collapse = ", "))
  }
}
assert_no_legacy_ids(pairs_df, c("Cm1", "Cm2", "Lsi", "Tan"), "pair syngroup table")
assert_no_legacy_ids(singles_df, c("Cm", "Lsi", "Tan"), "single syngroup table")

write.table(
  data.frame(
    file = c(basename(pairs_file), basename(singles_file)),
    total_rows = c(nrow(pairs_df), nrow(singles_df)),
    simple_rows = c(nrow(simple_pairs), nrow(simple_singles)),
    complex_rows = c(sum(trimws(pairs_df$complexity) == "complex"),
                     sum(trimws(singles_df$complexity) == "complex")),
    legacy_gene_id_check = "PASS",
    stringsAsFactors = FALSE
  ),
  "Fig2A_syngroup_input_audit.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

read_gene_matrix <- function(file, id_name, sample_names) {
  if (!file.exists(file)) stop("Input file does not exist: ", file)
  x <- read.delim(file, header = TRUE, check.names = FALSE,
                  quote = "", comment.char = "")
  if (ncol(x) < 2L) stop("Unexpected number of expression matrix columns: ", file)
  colnames(x)[1] <- id_name
  if (!all(sample_names %in% colnames(x))) {
    stop("Expression matrix lacks retained samples: ", file)
  }
  x <- x[, c(id_name, sample_names), drop = FALSE]
  if (anyDuplicated(x[[id_name]])) {
    stop("Expression matrix contains duplicate gene IDs: ", file)
  }
  value_mat <- as.matrix(x[, sample_names, drop = FALSE])
  storage.mode(value_mat) <- "numeric"
  if (any(!is.finite(value_mat)) || any(value_mat < 0)) {
    stop("Expression matrix contains NA/NaN/Inf or negative values: ", file)
  }
  x
}

read_scaled_branch <- function(branch) {
  branch_name <- if (branch == "TMM") "all_BS__TMM__" else "all_BS__median_ratio__"
  list(
    Cm = read_gene_matrix(
      file.path(normalization_dir, paste0(branch_name, "Cm_all_genes_scaled.tsv")),
      "Cm", retained_samples$Cm
    ),
    Ls = read_gene_matrix(
      file.path(normalization_dir, paste0(branch_name, "Ls_all_genes_scaled.tsv")),
      "Ls", retained_samples$Ls
    ),
    Ta = read_gene_matrix(
      file.path(normalization_dir, paste0(branch_name, "Ta_all_genes_scaled.tsv")),
      "Ta", retained_samples$Ta
    )
  )
}

# Unscaled TPM supplies the shared expression and MAD filters for both branches.
read_txi_abundance <- function(species) {
  txi_file <- file.path(tximport_dir, species, paste0(species, "_txi_scaledTPM.rds"))
  if (!file.exists(txi_file)) stop("Cannot find tximport abundance RDS: ", txi_file)
  txi <- readRDS(txi_file)
  if (!is.list(txi) || is.null(txi$abundance)) {
    stop("tximport RDS lacks abundance: ", txi_file)
  }
  abundance <- as.data.frame(txi$abundance, check.names = FALSE)
  colnames(abundance)[1:ncol(abundance)] <- colnames(txi$abundance)
  sample_names <- retained_samples[[species]]
  if (!all(sample_names %in% colnames(abundance))) {
    stop("Abundance RDS lacks retained samples: ", txi_file)
  }
  abundance$GID <- rownames(abundance)
  abundance <- abundance[, c("GID", sample_names), drop = FALSE]
  colnames(abundance)[1] <- species
  if (anyDuplicated(abundance[[species]])) stop("Abundance matrix contains duplicate gene IDs: ", txi_file)
  value_mat <- as.matrix(abundance[, sample_names, drop = FALSE])
  storage.mode(value_mat) <- "numeric"
  if (any(!is.finite(value_mat)) || any(value_mat < 0)) {
    stop("Abundance matrix contains NA/NaN/Inf or negative values: ", txi_file)
  }
  abundance
}

check_mapping_coverage <- function(mapping_ids, expr_df, id_col, label) {
  ids <- unique(mapping_ids[!is.na(mapping_ids) & mapping_ids != ""])
  missing_ids <- setdiff(ids, expr_df[[id_col]])
  if (length(missing_ids) > 0L) {
    stop(label, " contains gene IDs absent from the expression matrix. Examples: ",
         paste(head(missing_ids, 10), collapse = ", "))
  }
}

mean_numeric <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

aggregate_one_branch <- function(expr_list, branch_label) {
  cm <- expr_list$Cm
  ls <- expr_list$Ls
  ta <- expr_list$Ta

  check_mapping_coverage(c(simple_pairs$Cm1, simple_pairs$Cm2, simple_singles$Cm), cm, "Cm", "Cm")
  check_mapping_coverage(c(simple_pairs$Lsi, simple_singles$Lsi), ls, "Ls", "Ls")
  check_mapping_coverage(c(simple_pairs$Tan, simple_singles$Tan), ta, "Ta", "Ta")

  if (!"Cm" %in% names(cm)) {
    stop("The Cm expression matrix lacks the renamed Cm gene ID column; present columns: ",
         paste(names(cm), collapse = ", "))
  }
  if (!"Ls" %in% names(ls) || !"Ta" %in% names(ta)) {
    stop("The Ls/Ta expression matrix lacks the expected gene ID column.")
  }

  cm1 <- cm %>%
    dplyr::rename(Cm1 = Cm) %>%
    dplyr::rename_with(~ paste0("Cm1__", .x), dplyr::all_of(retained_samples$Cm))
  cm2 <- cm %>%
    dplyr::rename(Cm2 = Cm) %>%
    dplyr::rename_with(~ paste0("Cm2__", .x), dplyr::all_of(retained_samples$Cm))

  pair_joined <- simple_pairs %>%
    dplyr::mutate(unique_id = paste0("pair_", comp_id)) %>%
    dplyr::left_join(cm1, by = "Cm1") %>%
    dplyr::left_join(cm2, by = "Cm2") %>%
    dplyr::left_join(ta, by = c("Tan" = "Ta")) %>%
    dplyr::left_join(ls, by = c("Lsi" = "Ls"))

  pair_numeric <- c(
    paste0("Cm1__", retained_samples$Cm),
    paste0("Cm2__", retained_samples$Cm),
    retained_samples$Ls,
    retained_samples$Ta
  )
  # Average duplicate mapping rows before summing the two Cm copies.
  pair_agg <- pair_joined %>%
    dplyr::select(unique_id, Lsi, Tan, dplyr::all_of(pair_numeric)) %>%
    dplyr::group_by(unique_id, Lsi, Tan) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(pair_numeric), mean_numeric), .groups = "drop")

  pairs_final <- pair_agg %>% dplyr::select(unique_id)
  for (sample_name in retained_samples$Cm) {
    pairs_final[[sample_name]] <-
      pair_agg[[paste0("Cm1__", sample_name)]] +
      pair_agg[[paste0("Cm2__", sample_name)]]
  }
  pairs_final <- dplyr::bind_cols(
    pairs_final,
    pair_agg %>% dplyr::select(dplyr::all_of(retained_samples$Ls), dplyr::all_of(retained_samples$Ta))
  )

  single_joined <- simple_singles %>%
    dplyr::mutate(unique_id = paste0("single_", comp_id)) %>%
    dplyr::left_join(cm, by = "Cm") %>%
    dplyr::left_join(ta, by = c("Tan" = "Ta")) %>%
    dplyr::left_join(ls, by = c("Lsi" = "Ls"))
  single_numeric <- c(retained_samples$Cm, retained_samples$Ls, retained_samples$Ta)
  singles_final <- single_joined %>%
    dplyr::select(unique_id, Lsi, Tan, dplyr::all_of(single_numeric)) %>%
    dplyr::group_by(unique_id, Lsi, Tan) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(single_numeric), mean_numeric), .groups = "drop") %>%
    dplyr::select(unique_id, dplyr::all_of(single_numeric))

  total <- dplyr::bind_rows(
    pairs_final %>% dplyr::select(unique_id, dplyr::all_of(all_samples)),
    singles_final %>% dplyr::select(unique_id, dplyr::all_of(all_samples))
  )
  if (anyDuplicated(total$unique_id)) stop(branch_label, " contains duplicate unique_id values after aggregation.")
  value_mat <- as.matrix(total[, all_samples, drop = FALSE])
  storage.mode(value_mat) <- "numeric"
  if (any(!is.finite(value_mat)) || any(value_mat < 0)) {
    stop(branch_label, " syngroup expression matrix contains NA/NaN/Inf or negative values.")
  }
  total
}

branch_results <- list(
  TMM = aggregate_one_branch(read_scaled_branch("TMM"), "TMM"),
  median = aggregate_one_branch(read_scaled_branch("median"), "median")
)

tpm_results <- list(
  Cm = read_txi_abundance("Cm"),
  Ls = read_txi_abundance("Ls"),
  Ta = read_txi_abundance("Ta")
)
tpm_aggregated <- aggregate_one_branch(tpm_results, "TPM reference")

id_map <- dplyr::bind_rows(
  simple_pairs %>% dplyr::mutate(unique_id = paste0("pair_", comp_id)) %>%
    dplyr::select(unique_id, Cm1, Cm2, Lsi, Tan),
  simple_singles %>% dplyr::mutate(unique_id = paste0("single_", comp_id), Cm1 = Cm, Cm2 = NA_character_) %>%
    dplyr::select(unique_id, Cm1, Cm2, Lsi, Tan)
)
write.csv(id_map, "All_syngroup_ID_list.csv", row.names = FALSE, na = "")

time_mean <- function(expr_df, hour) {
  pattern <- if (hour == "00") "_0h_" else paste0("_", hour, "h_")
  out <- lapply(c("Cm", "Ls", "Ta"), function(species) {
    cols <- grep(paste0("^", species, pattern), all_samples, value = TRUE)
    if (length(cols) == 0L) stop("No samples found: ", species, " ", hour)
    rowMeans(as.matrix(expr_df[, cols, drop = FALSE]), na.rm = TRUE)
  })
  out <- as.data.frame(out, check.names = FALSE)
  colnames(out) <- paste0(c("Cm", "Ls", "Ta"), "_", hour)
  out$unique_id <- expr_df$unique_id
  out[, c("unique_id", paste0(c("Cm", "Ls", "Ta"), "_", hour)), drop = FALSE]
}

for (branch in names(branch_results)) {
  branch_output_dir <- file.path(output_dir, branch)
  write.table(
    branch_results[[branch]],
    file.path(branch_output_dir, paste0("All_syngroup_expr_", branch, "_sample_level.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  for (hour in c("00", "24", "48")) {
    write.table(
      time_mean(branch_results[[branch]], hour),
      file.path(branch_output_dir, paste0("All_syngroup_expr_", branch, "_", hour, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

write.table(
  tpm_aggregated,
  "All_syngroup_filter_TPM_sample_level.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
for (hour in c("00", "24", "48")) {
  write.table(
    time_mean(tpm_aggregated, hour),
    paste0("All_syngroup_filter_TPM_", hour, ".tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

write.table(
  data.frame(
    branch = c("TMM", "median"),
    input_prefix = c("all_BS__TMM__", "all_BS__median_ratio__"),
    role = c("Main_Fig2A_primary", "Fig2A_sensitivity"),
    stringsAsFactors = FALSE
  ),
  "cross_species_normalization_branch_status.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(sample = all_samples,
             species = sub("_.*", "", all_samples),
             time = sub("^[^_]+_([^_]+)_.*$", "\\1", all_samples)),
  "Fig2A_final_sample_audit.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "Fig2A_expression_extraction_sessionInfo.txt")
)

message("Expression extraction completed for the median and TMM branches.")
}))

.run_release_stage(quote({
RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)
options(stringsAsFactors = FALSE)

project_dir <- RELEASE_ROOT
work_dir <- file.path(project_dir, "2.global_expression/Fig2A_GO/kmeans")
dir.create(file.path(work_dir, "TMM"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "median"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work_dir, "sensitivity_comparison"), recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

library(ggplot2)
library(reshape2)
library(ggh4x)

METHODS <- c("TMM", "median")
HOURS <- c("00", "24", "48")
K_CLUSTERS <- 6L
KMEANS_NSTART <- 1L
MAD_QUANTILE <- 0.25
TPM_CUTOFF <- 1
MIN_SPECIES_DETECTED <- 2L

read_time_matrix <- function(method, hour) {
  method_dir <- file.path(work_dir, method)
  file <- file.path(method_dir, paste0("All_syngroup_expr_", method, "_", hour, ".tsv"))
  if (!file.exists(file)) stop("Expression matrix is missing: ", file)
  x <- read.delim(file, header = TRUE, row.names = 1, check.names = FALSE,
                  quote = "", comment.char = "")
  expected <- paste0(c("Cm", "Ls", "Ta"), "_", hour)
  if (!identical(colnames(x), expected)) {
    stop("Unexpected time-point expression matrix column names: ", file)
  }
  value_mat <- as.matrix(x)
  storage.mode(value_mat) <- "numeric"
  if (any(!is.finite(value_mat)) || any(value_mat < 0)) {
    stop("Expression matrix contains NA/NaN/Inf or negative values: ", file)
  }
  x
}

read_tpm_filter_matrix <- function(hour) {
  file <- file.path(work_dir, paste0("All_syngroup_filter_TPM_", hour, ".tsv"))
  if (!file.exists(file)) stop("TPM filtering matrix is missing: ", file)
  x <- read.delim(file, header = TRUE, row.names = 1, check.names = FALSE,
                  quote = "", comment.char = "")
  expected <- paste0(c("Cm", "Ls", "Ta"), "_", hour)
  if (!identical(colnames(x), expected)) stop("Unexpected TPM filtering matrix column names: ", file)
  value_mat <- as.matrix(x)
  storage.mode(value_mat) <- "numeric"
  if (any(!is.finite(value_mat)) || any(value_mat < 0)) {
    stop("TPM filtering matrix contains NA/NaN/Inf or negative values: ", file)
  }
  x
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  total_pairs <- choose2(sum(tab))
  if (total_pairs == 0) return(NA_real_)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / total_pairs
  denominator <- 0.5 * (sum_rows + sum_cols) - expected
  if (denominator == 0) return(NA_real_)
  (sum_cells - expected) / denominator
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_permutations(x[-i]))
  }))
}

# Match cluster labels by maximizing agreement across all label permutations.
best_cluster_mapping <- function(x, y) {
  tab <- table(x, y)
  x_levels <- sort(unique(x))
  y_levels <- sort(unique(y))
  permutations <- all_permutations(y_levels)
  scores <- apply(permutations, 1, function(candidate) {
    sum(mapply(function(a, b) tab[as.character(a), as.character(b)],
               x_levels, candidate))
  })
  best <- which.max(scores)
  list(
    matched_n = scores[best],
    mapping = stats::setNames(permutations[best, ], x_levels),
    contingency = tab
  )
}

run_syngroup_cluster <- function(method, hour) {
  expr <- read_time_matrix(method, hour)
  filter_tpm <- read_tpm_filter_matrix(hour)
  if (!identical(rownames(expr), rownames(filter_tpm))) {
    stop(method, " and the TPM filtering matrix have different syngroup orders: ", hour)
  }

  # Require mean TPM >= 1 in at least two species at each time point.
  keep_expr <- rowSums(filter_tpm >= TPM_CUTOFF, na.rm = TRUE) >= MIN_SPECIES_DETECTED
  filtered_tpm <- filter_tpm[keep_expr, , drop = FALSE]

  # Use log2 unscaled TPM for MAD; retain values strictly above its 25th percentile.
  log_tpm_reference <- log2(filtered_tpm + 1)
  mad_values <- apply(log_tpm_reference, 1, stats::mad, na.rm = TRUE)
  mad_cutoff <- stats::quantile(mad_values, MAD_QUANTILE, na.rm = TRUE)
  keep_mad <- names(mad_values)[is.finite(mad_values) & mad_values > mad_cutoff]

  log_expr <- log2(expr + 1)
  log_expr <- log_expr[keep_mad, , drop = FALSE]
  row_sd <- apply(log_expr, 1, stats::sd, na.rm = TRUE)
  keep_sd <- is.finite(row_sd) & row_sd > 0
  log_expr <- log_expr[keep_sd, , drop = FALSE]
  if (nrow(log_expr) < K_CLUSTERS) {
    stop(method, " ", hour, " hpi has fewer clusterable syngroups than k = ", K_CLUSTERS)
  }

  # Cluster row Z-scores after removing non-finite or zero row SD.
  expr_z <- t(scale(t(log_expr), center = TRUE, scale = TRUE))
  set.seed(1)
  km <- kmeans(expr_z, centers = K_CLUSTERS, nstart = KMEANS_NSTART)

  cluster_df <- as.data.frame(expr_z, check.names = FALSE)
  cluster_df$cluster <- as.integer(km$cluster)
  cluster_df <- cluster_df[order(cluster_df$cluster), , drop = FALSE]
  method_dir <- file.path(work_dir, method)
  cluster_file <- file.path(method_dir,
                             paste0("all_syngroup_expr_", method, "_", hour, "_cluster.txt"))
  write.table(cluster_df, cluster_file, sep = "\t", quote = FALSE)

  audit <- data.frame(
    method = method,
    hour = hour,
    n_input = nrow(expr),
    n_after_expression_filter = sum(keep_expr),
    n_after_mad_filter = length(keep_mad),
    n_after_row_sd_filter = nrow(expr_z),
    mad_cutoff = as.numeric(mad_cutoff),
    tpm_cutoff = TPM_CUTOFF,
    min_species_detected = MIN_SPECIES_DETECTED,
    mad_quantile = MAD_QUANTILE,
    k = K_CLUSTERS,
    seed = 1L,
    nstart = KMEANS_NSTART,
    stringsAsFactors = FALSE
  )
  write.table(audit,
              file.path(method_dir, paste0("kmeans_filter_audit_", method, "_", hour, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  species_names <- c("Cm", "Ls", "Ta")
  cluster_pattern <- do.call(rbind, lapply(sort(unique(km$cluster)), function(cluster_id) {
    members <- which(km$cluster == cluster_id)
    centroid <- colMeans(expr_z[members, , drop = FALSE])
    data.frame(
      method = method,
      hour = hour,
      cluster = cluster_id,
      n_syngroup = length(members),
      expression_order = paste(species_names[order(centroid, decreasing = TRUE)], collapse = ">"),
      Cm_centroid = unname(centroid[paste0("Cm_", hour)]),
      Ls_centroid = unname(centroid[paste0("Ls_", hour)]),
      Ta_centroid = unname(centroid[paste0("Ta_", hour)]),
      stringsAsFactors = FALSE
    )
  }))
  write.table(cluster_pattern,
              file.path(method_dir, paste0("cluster_pattern_", method, "_", hour, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  plot_df <- cluster_df
  plot_df$Row <- rownames(plot_df)
  plot_long <- reshape2::melt(plot_df, id.vars = c("Row", "cluster"))
  colnames(plot_long) <- c("Row", "Cluster", "Column", "Value")
  column_order <- paste0(c("Cm", "Ls", "Ta"), "_", hour)
  plot_long$Column <- factor(plot_long$Column, levels = column_order)
  plot_long$Cluster <- factor(plot_long$Cluster, levels = sort(unique(plot_long$Cluster)))

  axis_labels <- stats::setNames(c("Cm", "Ls", "Ta"), column_order)
  cluster_labels <- stats::setNames(paste0("C", sort(unique(km$cluster))),
                                    sort(unique(km$cluster)))

  p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = Column, y = Row, fill = Value)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::scale_fill_gradientn(
      colors = c("#1B66AE", "#4194C4", "#93C6DF", "#D2E6F0",
                 "#FDDCC8", "#D94021", "#B31127")
    ) +
    ggplot2::scale_x_discrete(labels = axis_labels) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 24, angle = 0, vjust = 0.5, hjust = 0.5,
        family = "Arial", face = "italic", color = "black"
      ),
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(size = 22, family = "Arial", angle = 0),
      plot.title = ggplot2::element_text(size = 24, family = "Arial", hjust = 0),
      panel.spacing = grid::unit(0.1, "lines")
    ) +
    ggplot2::labs(title = paste0(ifelse(hour == "00", "0", hour), " hpi")) +
    ggh4x::facet_grid2(
      Cluster ~ ., scales = "free", space = "free",
      labeller = ggplot2::labeller(Cluster = cluster_labels)
    )

  pdf_file <- file.path(
    method_dir,
    paste0("sumCm_all_Cluster_", method, "_", ifelse(hour == "00", "0", hour),
           "hpi_heatmap_anti2.0.pdf")
  )
  ggplot2::ggsave(pdf_file, p, width = 4, height = 4,
                  device = grDevices::cairo_pdf)

  invisible(list(cluster = cluster_df, pattern = cluster_pattern, plot = p))
}

for (method in METHODS) {
  for (hour in HOURS) {
    run_syngroup_cluster(method, hour)
  }
}

comparison_rows <- list()
for (hour in HOURS) {
  tmm_file <- file.path(work_dir, "TMM", paste0("all_syngroup_expr_TMM_", hour, "_cluster.txt"))
  median_file <- file.path(work_dir, "median", paste0("all_syngroup_expr_median_", hour, "_cluster.txt"))
  tmm <- read.delim(tmm_file, header = TRUE, row.names = 1, check.names = FALSE)
  median <- read.delim(median_file, header = TRUE, row.names = 1, check.names = FALSE)
  common_ids <- intersect(rownames(tmm), rownames(median))
  union_ids <- union(rownames(tmm), rownames(median))
  tmm_cluster <- as.integer(tmm[common_ids, "cluster"])
  median_cluster <- as.integer(median[common_ids, "cluster"])
  mapping <- best_cluster_mapping(tmm_cluster, median_cluster)

  comparison <- data.frame(
    hour = hour,
    n_TMM = nrow(tmm),
    n_median = nrow(median),
    n_common = length(common_ids),
    gene_set_Jaccard = length(common_ids) / length(union_ids),
    adjusted_Rand_index = adjusted_rand_index(tmm_cluster, median_cluster),
    best_label_matched_n = mapping$matched_n,
    best_label_match_rate = mapping$matched_n / length(common_ids),
    role = "TMM_primary_median_sensitivity",
    stringsAsFactors = FALSE
  )
  comparison_rows[[hour]] <- comparison
  write.table(comparison,
              file.path(work_dir, "sensitivity_comparison",
                        paste0("TMM_vs_median_cluster_comparison_", hour, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(mapping$contingency,
              file.path(work_dir, "sensitivity_comparison",
                        paste0("TMM_vs_median_cluster_contingency_", hour, ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
}
write.table(do.call(rbind, comparison_rows),
            file.path(work_dir, "sensitivity_comparison",
                      "TMM_vs_median_cluster_comparison_all_hours.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(work_dir, "Fig2A_kmeans_sessionInfo.txt")
)

message("K-means clustering and heatmaps completed for the median and TMM branches.")
}))
