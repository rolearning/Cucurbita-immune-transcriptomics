#!/usr/bin/env Rscript

# Summarize genome-wide CmB/CmS expression bias for Table S3 from the final Cm
# sample manifest, scaledTPM count-like matrix, and corrected simple BS syngroup map.
# Fit paired DESeq2 models at each timepoint using all-gene poscounts size factors.
# Export pair-level results, summary/QC tables, and TableS3_NEW.xlsx.

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(openxlsx)
})

project_root <- RELEASE_ROOT
out_dir <- file.path(project_root, "3.SC_index/FigS6_TableS3S6/TableS3_NEW")
sample_info_path <- file.path(project_root, "0.salmon/sample_info.tsv")
counts_path <- file.path(
  project_root,
  "1.normalization/results/01_tximport/Cm/Cm_scaledTPM_counts.tsv"
)
pairs_path <- file.path(
  project_root,
  "1.normalization/All_pairs_syngroup_w_notes_gene_corrected.csv"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(sample_info_path, counts_path, pairs_path)
if (any(!file.exists(required))) {
  stop("Missing required project input(s):\n", paste(required[!file.exists(required)], collapse = "\n"))
}

samples <- readr::read_tsv(sample_info_path, show_col_types = FALSE,
                           col_types = cols(.default = col_character()))
required_sample_cols <- c("sample", "species", "timepoint", "replicate")
if (!all(required_sample_cols %in% names(samples))) {
  stop("sample_info.tsv is missing: ", paste(setdiff(required_sample_cols, names(samples)), collapse = ", "))
}
cm_samples <- samples %>% filter(species == "Cm")
if (nrow(cm_samples) != 8L || anyDuplicated(cm_samples$sample)) {
  stop("Current sample_info.tsv must contain exactly 8 unique Cm libraries.")
}
if (any(cm_samples$sample %in% c("Cm_48h_D2"))) {
  stop("Cm_48h_D2 must not be present in the final Cm sample set.")
}
if (any(is.na(cm_samples$timepoint) | is.na(cm_samples$replicate) |
        cm_samples$timepoint == "" | cm_samples$replicate == "")) {
  stop("Cm sample_info contains missing timepoint or replicate values.")
}

counts <- readr::read_tsv(counts_path, show_col_types = FALSE,
                          col_types = cols(.default = col_double(), GID = col_character())) %>%
  as.data.frame(check.names = FALSE)
if (!"GID" %in% names(counts)) stop("Counts matrix must contain GID.")
if (anyDuplicated(counts$GID)) stop("Counts matrix contains duplicated GID.")
missing_samples <- setdiff(cm_samples$sample, names(counts))
if (length(missing_samples)) stop("Counts matrix is missing Cm sample(s): ", paste(missing_samples, collapse = ", "))
count_mat_full <- as.matrix(counts[, cm_samples$sample, drop = FALSE])
storage.mode(count_mat_full) <- "double"
if (any(!is.finite(count_mat_full)) || any(count_mat_full < 0)) stop("Counts matrix contains invalid values.")

# Estimate factors across all Cm genes before restricting to simple BS pairs.
global_sf <- DESeq2::estimateSizeFactorsForMatrix(count_mat_full, type = "poscounts")
names(global_sf) <- cm_samples$sample
if (any(!is.finite(global_sf)) || any(global_sf <= 0)) stop("Invalid global poscounts size factors.")
readr::write_tsv(tibble(sample = names(global_sf), poscounts_size_factor = as.numeric(global_sf)),
                 file.path(out_dir, "Cm_global_poscounts_size_factors.tsv"))

pairs <- read.csv(pairs_path, check.names = FALSE, stringsAsFactors = FALSE)
required_pair_cols <- c("Group", "Cm1", "Cm2", "Source1", "Source2", "Type", "complexity")
if (!all(required_pair_cols %in% names(pairs))) {
  stop("Current syngroup file is missing: ", paste(setdiff(required_pair_cols, names(pairs)), collapse = ", "))
}
simple_bs <- pairs %>% filter(Type == "BS", complexity == "simple")
n_simple_bs <- nrow(simple_bs)
if (n_simple_bs == 0L) stop("No simple BS rows found in current syngroup file.")

bs_map <- simple_bs %>%
  mutate(
    CmB = case_when(Source1 == "BENI" & Source2 == "SICY" ~ Cm1,
                    Source2 == "BENI" & Source1 == "SICY" ~ Cm2,
                    TRUE ~ NA_character_),
    CmS = case_when(Source1 == "BENI" & Source2 == "SICY" ~ Cm2,
                    Source2 == "BENI" & Source1 == "SICY" ~ Cm1,
                    TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(CmB), !is.na(CmS)) %>%
  transmute(Group, CmB, CmS, Source1, Source2, Type, complexity)

duplicated_group_n <- sum(duplicated(bs_map$Group))
if (duplicated_group_n > 0L) {
  stop("Duplicated Group values among valid simple BS pairs: ", duplicated_group_n)
}
missing_expression_genes <- setdiff(unique(c(bs_map$CmB, bs_map$CmS)), counts$GID)
if (length(missing_expression_genes) > 0L) {
  stop("BS pair gene(s) absent from Cm counts matrix: ", paste(head(missing_expression_genes, 20), collapse = ", "))
}
readr::write_tsv(bs_map, file.path(out_dir, "pairs_CmBCmS_list_NEW.tsv"))

timepoints <- unique(as.character(cm_samples$timepoint))
timepoints <- timepoints[order(suppressWarnings(as.numeric(timepoints)))]
if (!all(timepoints %in% c("00", "24", "48"))) stop("Unexpected Cm timepoint in sample_info: ", paste(timepoints, collapse = ", "))

status_levels <- c("No bias", "B dominant", "S dominant", "Undetermined")
result_list <- list()
qc_list <- list()

for (tp in timepoints) {
  smp <- cm_samples %>% filter(timepoint == tp)
  if (nrow(smp) < 2L) stop("Fewer than two Cm samples at timepoint ", tp)
  if (anyDuplicated(smp$replicate)) stop("Replicate labels are duplicated within timepoint ", tp)

  b <- count_mat_full[match(bs_map$CmB, counts$GID), smp$sample, drop = FALSE]
  s <- count_mat_full[match(bs_map$CmS, counts$GID), smp$sample, drop = FALSE]
  combined <- cbind(b, s)
  colnames(combined) <- c(paste0("B_", smp$sample), paste0("S_", smp$sample))
  rownames(combined) <- bs_map$Group
  if (any(!is.finite(combined))) stop("Non-finite B/S count value at timepoint ", tp)

  cd <- data.frame(
    rep = factor(rep(smp$replicate, 2)),
    subgenome = factor(c(rep("BENI", nrow(smp)), rep("SICY", nrow(smp))),
                       levels = c("SICY", "BENI")),
    row.names = colnames(combined)
  )
  sf <- c(global_sf[smp$sample], global_sf[smp$sample])
  names(sf) <- colnames(combined)

  # Retain fractional counts until model construction; block on matched biological replicates.
  dds <- DESeqDataSetFromMatrix(
    countData = round(combined),
    colData = cd,
    design = ~ rep + subgenome
  )
  sizeFactors(dds) <- sf
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("subgenome", "BENI", "SICY"))) %>%
    rownames_to_column("Group") %>%
    transmute(
      Group,
      timepoint = tp,
      log2FC = log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      # Bias requires both |log2FC| >= 1 and padj < 0.05; missing tests remain undetermined.
      Status = case_when(
        is.na(log2FC) | is.na(padj) ~ "Undetermined",
        log2FC >= 1 & padj < 0.05 ~ "B dominant",
        log2FC <= -1 & padj < 0.05 ~ "S dominant",
        TRUE ~ "No bias"
      )
    )
  result_list[[tp]] <- res
  summary_tp <- res %>% count(Status, name = "n") %>%
    right_join(tibble(Status = status_levels), by = "Status") %>%
    mutate(n = coalesce(n, 0L), proportion = n / sum(n), timepoint = tp) %>%
    select(timepoint, Status, n, proportion)
  qc_list[[tp]] <- tibble(
    timepoint = tp,
    n_timepoint_samples = nrow(smp),
    timepoint_samples = paste(smp$sample, collapse = ";"),
    n_valid_CmB_CmS_pairs = nrow(bs_map),
    n_entered_DESeq2 = nrow(res),
    n_NA_log2FC = sum(is.na(res$log2FC)),
    n_NA_padj = sum(is.na(res$padj)),
    n_no_bias = sum(res$Status == "No bias"),
    n_B_dominant = sum(res$Status == "B dominant"),
    n_S_dominant = sum(res$Status == "S dominant"),
    n_undetermined = sum(res$Status == "Undetermined")
  )
  readr::write_tsv(res, file.path(out_dir, paste0("TableS3_", tp, "hpi_NEW.tsv")))
  readr::write_tsv(summary_tp, file.path(out_dir, paste0("TableS3_summary_", tp, "hpi_NEW.tsv")))
}

summary_all <- bind_rows(lapply(result_list, function(x) {
  x %>% count(Status, name = "n") %>% right_join(tibble(Status = status_levels), by = "Status") %>%
    mutate(n = coalesce(n, 0L), proportion = n / sum(n), timepoint = unique(x$timepoint)) %>%
    select(timepoint, Status, n, proportion)
}))
readr::write_tsv(summary_all, file.path(out_dir, "TableS3_summary_NEW.tsv"))
readr::write_tsv(bind_rows(qc_list) %>% mutate(
  n_simple_BS_rows_before_mapping = n_simple_bs,
  n_valid_pairs_after_mapping = nrow(bs_map),
  duplicated_Group_count = duplicated_group_n,
  n_missing_expression_genes = length(missing_expression_genes)
), file.path(out_dir, "TableS3_QC_NEW.tsv"))

summary_for_excel <- bind_rows(lapply(result_list, function(x) {
  x %>%
    count(Status, name = "Count") %>%
    right_join(tibble(Status = status_levels), by = "Status") %>%
    mutate(
      Count = coalesce(Count, 0L),
      Proportion = Count / sum(Count),
      timepoint = unique(x$timepoint)
    ) %>%
    select(Status, timepoint, Count, Proportion)
}))

excel_table <- tibble(Status = status_levels)
for (tp in c("00", "24", "48")) {
  tp_data <- summary_for_excel %>% filter(timepoint == tp)
  excel_table[[paste0(tp, " hpi Count")]] <-
    tp_data$Count[match(status_levels, tp_data$Status)]
  excel_table[[paste0(tp, " hpi Proportion")]] <-
    tp_data$Proportion[match(status_levels, tp_data$Status)]
}

excel_path <- file.path(out_dir, "TableS3_NEW.xlsx")
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "S3-genome-wide-bias")

title_text <- paste0(
  "Supplementary Table S3 Genome-wide distribution for subgenome expression bias of ",
  nrow(bs_map), " BS pairs in Cucurbita moschata (Cm)."
)
openxlsx::mergeCells(wb, "S3-genome-wide-bias", cols = 1:7, rows = 1)
openxlsx::writeData(wb, "S3-genome-wide-bias", title_text, startCol = 1, startRow = 1)

openxlsx::mergeCells(wb, "S3-genome-wide-bias", cols = 2:3, rows = 2)
openxlsx::mergeCells(wb, "S3-genome-wide-bias", cols = 4:5, rows = 2)
openxlsx::mergeCells(wb, "S3-genome-wide-bias", cols = 6:7, rows = 2)
openxlsx::writeData(
  wb, "S3-genome-wide-bias",
  matrix(c("", "0 hpi", "", "24 hpi", "", "48 hpi", ""), nrow = 1),
  startCol = 1, startRow = 2, colNames = FALSE
)
openxlsx::writeData(
  wb, "S3-genome-wide-bias",
  matrix(c("Status", "Count", "Proportion", "Count", "Proportion", "Count", "Proportion"), nrow = 1),
  startCol = 1, startRow = 3, colNames = FALSE
)
openxlsx::writeData(
  wb, "S3-genome-wide-bias",
  excel_table,
  startCol = 1, startRow = 4, colNames = FALSE, rowNames = FALSE
)

title_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 11, textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
header_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 11, textDecoration = "bold",
  halign = "center", valign = "center", border = "Bottom"
)
body_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 11)
number_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 11, numFmt = "0.0000"
)

openxlsx::addStyle(wb, "S3-genome-wide-bias", title_style,
                   rows = 1, cols = 1:7, gridExpand = TRUE)
openxlsx::addStyle(wb, "S3-genome-wide-bias", header_style,
                   rows = 2:3, cols = 1:7, gridExpand = TRUE)
openxlsx::addStyle(wb, "S3-genome-wide-bias", body_style,
                   rows = 4:(3 + nrow(excel_table)), cols = 1:7, gridExpand = TRUE)
openxlsx::addStyle(wb, "S3-genome-wide-bias", number_style,
                   rows = 4:(3 + nrow(excel_table)), cols = c(3, 5, 7), gridExpand = TRUE)
openxlsx::setColWidths(wb, "S3-genome-wide-bias", cols = 1, widths = 24)
openxlsx::setColWidths(wb, "S3-genome-wide-bias", cols = 2:7, widths = 15)
openxlsx::setRowHeights(wb, "S3-genome-wide-bias", rows = 1, heights = 34)
openxlsx::freezePane(wb, "S3-genome-wide-bias", firstActiveRow = 4)
openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)

message("Table S3 TSV and Excel outputs written: ", excel_path)
