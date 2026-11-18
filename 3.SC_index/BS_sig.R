#!/usr/bin/env Rscript

# Calculate copy-level CmB/CmS bias and produce Table S6 and Supplementary Fig. 6.
# Read the final sample manifest, OG map, normalization QC, and scaledTPM counts;
# use corrected tximport TPM abundance for the heatmap colors.
# Fit paired DESeq2 models per timepoint with global ratio size factors.
# Export Table S6 Excel/TSV files, per-module heatmaps, and a shared legend.

.run_release_stage <- function(code) {
  prior_wd <- getwd()
  on.exit(setwd(prior_wd), add = TRUE)
  eval(code, envir = new.env(parent = globalenv()))
}

.run_release_stage(quote({

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(readxl)
  library(tidyr)
  library(tibble)
  library(openxlsx)
})

project_root <- RELEASE_ROOT
out_dir <- file.path(project_root, "3.SC_index/FigS6_TableS3S6/TableS6_NEW")
sample_info_path <- file.path(project_root, "0.salmon/sample_info.tsv")
counts_path <- file.path(project_root, "1.normalization/results/01_tximport/Cm/Cm_scaledTPM_counts.tsv")
qc_path <- file.path(project_root, "1.normalization/results/01_tximport/all_species_sample_normalization_QC.tsv")
og_path <- file.path(project_root, "OG_info/OG_info_260906_OK.xlsx")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(sample_info_path, counts_path, qc_path, og_path)
if (any(!file.exists(required))) stop("Missing required input(s):\n", paste(required[!file.exists(required)], collapse = "\n"))

samples <- read_tsv(sample_info_path, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
  filter(species == "Cm")
if (nrow(samples) != 8L || anyDuplicated(samples$sample) || any(samples$sample == "Cm_48h_D2")) {
  stop("sample_info.tsv does not contain the required final 8 Cm libraries.")
}
counts <- read_tsv(counts_path, show_col_types = FALSE,
                   col_types = cols(.default = col_double(), GID = col_character())) %>% as.data.frame(check.names = FALSE)
if (!"GID" %in% names(counts) || anyDuplicated(counts$GID)) stop("Invalid Cm counts matrix GID column.")
if (length(setdiff(samples$sample, names(counts)))) stop("Counts matrix is missing final Cm sample(s).")
full_mat <- as.matrix(counts[, samples$sample, drop = FALSE]); storage.mode(full_mat) <- "double"
if (any(!is.finite(full_mat)) || any(full_mat < 0)) stop("Counts matrix contains invalid values.")

og <- readxl::read_excel(og_path, sheet = 1, .name_repair = "minimal") %>% as.data.frame(check.names = FALSE)
required_og <- c("Cmo", "GeneName", "GeneName_base", "Module", "subgenome", "Type")
if (!all(required_og %in% names(og))) stop("OG_info is missing: ", paste(setdiff(required_og, names(og)), collapse = ", "))
og <- og %>% mutate(across(everything(), as.character))
if (any(!og$subgenome %in% c("CmB", "CmS"))) stop("OG_info contains unexpected subgenome values.")
if (any(!og$Type %in% c("BS", "B", "S", "BB", "SS"))) stop("OG_info contains unexpected Type values.")

# Test true BS pairs only; BB/SS records are audited and later displayed without significance tests.
pair_audit <- og %>%
  group_by(GeneName_base) %>%
  summarise(
    n_CmB = sum(subgenome == "CmB"), n_CmS = sum(subgenome == "CmS"),
    n_BS_type_CmB = sum(subgenome == "CmB" & Type == "BS"),
    n_BS_type_CmS = sum(subgenome == "CmS" & Type == "BS"),
    .groups = "drop"
  ) %>%
  mutate(
    eligible_true_BS = n_BS_type_CmB > 0 & n_BS_type_CmS > 0,
    excluded_reason = case_when(
      eligible_true_BS ~ "eligible",
      n_CmB == 0 ~ "no CmB counterpart",
      n_CmS == 0 ~ "no CmS counterpart",
      TRUE ~ "not both sides marked Type=BS"
    )
  )
readr::write_tsv(pair_audit, file.path(out_dir, "TableS6_BS_BB_SS_inclusion_audit_NEW.tsv"))

bs <- og %>% filter(subgenome == "CmB", Type == "BS") %>%
  select(GeneName_base, Cmo_B = Cmo, GeneName_B = GeneName, Module_B = Module)
ss <- og %>% filter(subgenome == "CmS", Type == "BS") %>%
  select(GeneName_base, Cmo_S = Cmo, GeneName_S = GeneName, Module_S = Module)
pairs <- inner_join(bs, ss, by = "GeneName_base") %>%
  mutate(
    PairID = paste(GeneName_base, Cmo_B, Cmo_S, sep = "|"),
    Module = ifelse(Module_B == Module_S, Module_B, NA_character_)
  ) %>%
  relocate(PairID)
if (anyDuplicated(pairs$PairID)) stop("Duplicated S6 PairID.")
if (any(is.na(pairs$Module))) stop("B/S copy pair has inconsistent Module assignment.")
missing_genes <- setdiff(unique(c(pairs$Cmo_B, pairs$Cmo_S)), counts$GID)
if (length(missing_genes)) stop("S6 B/S copy gene(s) absent from counts matrix: ", paste(missing_genes, collapse = ", "))
write_tsv(pairs, file.path(out_dir, "TableS6_BS_copy_pairs_NEW.tsv"))

# Estimate ratio factors over all final Cm genes and verify them against upstream QC.
sf_ratio_recomputed <- estimateSizeFactorsForMatrix(full_mat, type = "ratio")
names(sf_ratio_recomputed) <- samples$sample
qc <- read_tsv(qc_path, show_col_types = FALSE)
if (!all(c("sample", "species", "scaledTPM_size_factor") %in% names(qc))) {
  stop("Formal sample normalization QC lacks scaledTPM_size_factor.")
}
qc_cm <- qc %>% filter(species == "Cm") %>% select(sample, scaledTPM_size_factor)
if (setequal(qc_cm$sample, samples$sample) == FALSE) stop("QC sample set does not match final Cm sample set.")
qc_sf <- setNames(qc_cm$scaledTPM_size_factor, qc_cm$sample)[samples$sample]
if (any(!is.finite(qc_sf)) || any(qc_sf <= 0)) stop("Invalid ratio factors in formal QC.")
if (max(abs(as.numeric(sf_ratio_recomputed) - as.numeric(qc_sf))) > 1e-8) {
  stop("Recomputed ratio factors do not match all_species_sample_normalization_QC.tsv beyond numerical tolerance.")
}
write_tsv(tibble(sample = samples$sample, ratio_size_factor_recomputed = as.numeric(sf_ratio_recomputed),
                 ratio_size_factor_QC = as.numeric(qc_sf)),
          file.path(out_dir, "Cm_global_ratio_size_factors_NEW.tsv"))

timepoints <- unique(as.character(samples$timepoint))
timepoints <- timepoints[order(suppressWarnings(as.numeric(timepoints)))]
long <- list(); qc_out <- list()
for (tp in timepoints) {
  smp <- samples %>% filter(timepoint == tp)
  if (nrow(smp) < 2L || anyDuplicated(smp$replicate)) stop("Invalid replicate set at timepoint ", tp)
  b <- full_mat[match(pairs$Cmo_B, counts$GID), smp$sample, drop = FALSE]
  s <- full_mat[match(pairs$Cmo_S, counts$GID), smp$sample, drop = FALSE]
  mat <- cbind(b, s)
  colnames(mat) <- c(paste0("B_", smp$sample), paste0("S_", smp$sample))
  rownames(mat) <- pairs$PairID
  cd <- data.frame(
    rep = factor(rep(smp$replicate, 2)),
    condition = factor(c(rep("B", nrow(smp)), rep("S", nrow(smp))), levels = c("S", "B")),
    row.names = colnames(mat)
  )
  # Each B/S observation from the same library shares its externally estimated size factor.
  sf <- c(qc_sf[smp$sample], qc_sf[smp$sample]); names(sf) <- colnames(mat)
  # Round only at model construction; the blocking factor pairs B and S within each replicate.
  dds <- DESeqDataSetFromMatrix(round(mat), cd, design = ~ rep + condition)
  sizeFactors(dds) <- sf
  dds <- DESeq(dds, quiet = TRUE)
  # Positive log2FC favors B; stars require |log2FC| >= 1 and padj < 0.05.
  rr <- as.data.frame(results(dds, contrast = c("condition", "B", "S"))) %>%
    rownames_to_column("PairID") %>%
    transmute(
      PairID, timepoint = tp, log2FC = log2FoldChange, lfcSE, stat, pvalue, padj,
      direction = case_when(is.na(log2FC) ~ "Undetermined", log2FC > 0 ~ "B>S", log2FC < 0 ~ "S>B", TRUE ~ "Equal"),
      star = ifelse(!is.na(log2FC) & !is.na(padj) & abs(log2FC) >= 1 & padj < 0.05, "*", "")
    )
  long[[tp]] <- rr
  qc_out[[tp]] <- tibble(timepoint = tp, n_samples = nrow(smp), samples = paste(smp$sample, collapse = ";"),
                          n_pairs = nrow(rr), n_NA_log2FC = sum(is.na(rr$log2FC)), n_NA_padj = sum(is.na(rr$padj)))
}

long_df <- bind_rows(long)
wide <- long_df %>%
  pivot_wider(id_cols = PairID, names_from = timepoint,
              values_from = c(log2FC, lfcSE, stat, pvalue, padj, direction, star),
              names_glue = "{.value}_{timepoint}") %>%
  left_join(pairs, by = "PairID") %>%
  relocate(PairID, GeneName_base, Cmo_B, GeneName_B, Cmo_S, GeneName_S, Module)
write_tsv(wide, file.path(out_dir, "TableS6_copy_level_BS_bias_NEW.tsv"))
write_tsv(long_df, file.path(out_dir, "TableS6_copy_level_BS_bias_long_NEW.tsv"))
write_tsv(bind_rows(qc_out) %>% mutate(n_BS_pairs = nrow(pairs), n_BB_or_SS_rows_not_tested = sum(!pair_audit$eligible_true_BS)),
          file.path(out_dir, "TableS6_QC_NEW.tsv"))

table_s6_formatted <- wide %>%
  transmute(
    OG = GeneName_base,
    GeneName_CmB = GeneName_B,
    GeneName_CmS = GeneName_S,
    Module,
    Type = "BS",
    log2FC_0hpi = log2FC_00,
    log2FC_24hpi = log2FC_24,
    log2FC_48hpi = log2FC_48,
    padj_0hpi = padj_00,
    padj_24hpi = padj_24,
    padj_48hpi = padj_48,
    sig_0hpi = ifelse(!is.na(star_00), star_00, ""),
    sig_24hpi = ifelse(!is.na(star_24), star_24, ""),
    sig_48hpi = ifelse(!is.na(star_48), star_48, ""),
    Bias_0hpi = ifelse(!is.na(direction_00), direction_00, ""),
    Bias_24hpi = ifelse(!is.na(direction_24), direction_24, ""),
    Bias_48hpi = ifelse(!is.na(direction_48), direction_48, "")
  )

formatted_tsv <- file.path(out_dir, "TableS6_copy_level_BS_bias_formatted_NEW.tsv")
write_tsv(table_s6_formatted, formatted_tsv, na = "")

formatted_xlsx <- file.path(out_dir, "TableS6_NEW.xlsx")
wb_s6 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_s6, "Table S6")
openxlsx::writeData(wb_s6, "Table S6", table_s6_formatted, startRow = 1, startCol = 1)
header_style_s6 <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 10, textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE,
  border = "Bottom"
)
body_style_s6 <- openxlsx::createStyle(fontName = "Arial", fontSize = 10)
padj_style_s6 <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, numFmt = "0.0000E+00")
openxlsx::addStyle(wb_s6, "Table S6", header_style_s6,
                   rows = 1, cols = seq_len(ncol(table_s6_formatted)), gridExpand = TRUE)
openxlsx::addStyle(wb_s6, "Table S6", body_style_s6,
                   rows = 2:(nrow(table_s6_formatted) + 1),
                   cols = seq_len(ncol(table_s6_formatted)), gridExpand = TRUE)
openxlsx::addStyle(wb_s6, "Table S6", padj_style_s6,
                   rows = 2:(nrow(table_s6_formatted) + 1),
                   cols = c(9, 10, 11), gridExpand = TRUE)
openxlsx::setColWidths(wb_s6, "Table S6", cols = 1, widths = 22)
openxlsx::setColWidths(wb_s6, "Table S6", cols = 2:5, widths = 20)
openxlsx::setColWidths(wb_s6, "Table S6", cols = 6:11, widths = 15)
openxlsx::setColWidths(wb_s6, "Table S6", cols = 12:17, widths = 12)
openxlsx::freezePane(wb_s6, "Table S6", firstActiveRow = 2)
openxlsx::saveWorkbook(wb_s6, formatted_xlsx, overwrite = TRUE)

message("Table S6 analysis and export completed.")
}))

# Plot corrected TPM abundance, using Table S6 only for BS-pair significance stars.
.run_release_stage(quote({

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(grid)
})

project_root <- RELEASE_ROOT
out_dir <- file.path(project_root, "3.SC_index/FigS6_TableS3S6/FigS6_NEW")
sample_info_path <- file.path(project_root, "0.salmon/sample_info.tsv")
expression_rds_path <- file.path(project_root, "1.normalization/results/01_tximport/Cm/Cm_txi_original_counts_offset.rds")
counts_path <- file.path(project_root, "1.normalization/results/01_tximport/Cm/Cm_scaledTPM_counts.tsv")
og_path <- file.path(project_root, "OG_info/OG_info_260906_OK.xlsx")
table_s6_path <- file.path(project_root, "3.SC_index/FigS6_TableS3S6/TableS6_NEW/TableS6_copy_level_BS_bias_NEW.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(sample_info_path, expression_rds_path, counts_path, og_path, table_s6_path)
if (any(!file.exists(required))) stop("Missing required Fig.S6 input(s):\n", paste(required[!file.exists(required)], collapse = "\n"))

samples <- read_tsv(sample_info_path, show_col_types = FALSE, col_types = cols(.default = col_character())) %>% filter(species == "Cm")
if (nrow(samples) != 8L || anyDuplicated(samples$sample) || any(samples$sample == "Cm_48h_D2")) stop("Final Cm sample audit failed.")
txi <- readRDS(expression_rds_path)
if (!is.list(txi) || !all(c("counts", "abundance", "length") %in% names(txi))) {
  stop("Corrected tximport RDS must contain counts, abundance, and length.")
}
if (is.null(rownames(txi$abundance)) || anyDuplicated(rownames(txi$abundance))) {
  stop("Corrected tximport abundance must have unique gene IDs as row names.")
}
if (!identical(rownames(txi$abundance), rownames(txi$counts)) ||
    !identical(colnames(txi$abundance), colnames(txi$counts))) {
  stop("Corrected tximport counts and abundance dimensions/IDs do not match.")
}
if (length(setdiff(samples$sample, colnames(txi$abundance)))) {
  stop("Corrected tximport abundance is missing final Cm sample(s).")
}
expr <- as.data.frame(txi$abundance[, samples$sample, drop = FALSE], check.names = FALSE)
expr <- cbind(GID = rownames(expr), expr)
expr_mat <- as.matrix(expr[, samples$sample, drop = FALSE]); storage.mode(expr_mat) <- "double"
if (any(!is.finite(expr_mat)) || any(expr_mat < 0)) stop("Corrected tximport TPM abundance contains invalid values.")
if (any(colSums(expr_mat) <= 0)) stop("Corrected tximport TPM abundance has a non-positive sample total.")

# Missing TPM rows may be restored only if the count-like reference is zero in every final library.
counts_ref <- read_tsv(counts_path, show_col_types = FALSE,
                       col_types = cols(.default = col_double(), GID = col_character())) %>%
  as.data.frame(check.names = FALSE)
if (!"GID" %in% names(counts_ref) || anyDuplicated(counts_ref$GID)) stop("Count reference must contain unique GID.")
if (length(setdiff(samples$sample, names(counts_ref)))) stop("Count reference is missing final Cm sample(s).")
counts_ref_mat <- as.matrix(counts_ref[, samples$sample, drop = FALSE]); storage.mode(counts_ref_mat) <- "double"
if (any(!is.finite(counts_ref_mat)) || any(counts_ref_mat < 0)) stop("Count reference contains invalid values.")

og <- readxl::read_excel(og_path, sheet = 1, .name_repair = "minimal") %>% as.data.frame(check.names = FALSE) %>% mutate(across(everything(), as.character))
required_og <- c("Cmo", "GeneName", "GeneName_base", "Module", "subgenome", "Type")
if (!all(required_og %in% names(og))) stop("OG_info is missing: ", paste(setdiff(required_og, names(og)), collapse = ", "))
if (any(!og$subgenome %in% c("CmB", "CmS")) || any(!og$Type %in% c("BS", "B", "S", "BB", "SS"))) stop("OG_info contains unexpected subgenome/Type values.")

b <- og %>% filter(subgenome == "CmB") %>% select(GeneName_base, Cmo_B = Cmo, GeneName_B = GeneName, Module_B = Module, Type_B = Type)
s <- og %>% filter(subgenome == "CmS") %>% select(GeneName_base, Cmo_S = Cmo, GeneName_S = GeneName, Module_S = Module, Type_S = Type)
# BB/SS display pairs place the second same-subgenome copy on the right, preserving its B/S color role.
make_same_side_pairs <- function(x, role) {
  if (nrow(x) < 2L) return(NULL)
  out <- lapply(split(x, x$GeneName_base), function(z) {
    if (nrow(z) < 2L) return(NULL)
    ij <- utils::combn(seq_len(nrow(z)), 2L)
    do.call(rbind, lapply(seq_len(ncol(ij)), function(k) {
      i <- ij[1L, k]; j <- ij[2L, k]
      if (role == "B") {
        data.frame(
          GeneName_base = z$GeneName_base[i],
          Cmo_B = z$Cmo_B[i], GeneName_B = z$GeneName_B[i], Module_B = z$Module_B[i], Type_B = z$Type_B[i],
          Cmo_S = z$Cmo_B[j], GeneName_S = z$GeneName_B[j], Module_S = z$Module_B[j], Type_S = z$Type_B[j],
          ColorRole_B = "B", ColorRole_S = "B", PairDisplay = "BB-paired",
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          GeneName_base = z$GeneName_base[i],
          Cmo_B = z$Cmo_S[i], GeneName_B = z$GeneName_S[i], Module_B = z$Module_S[i], Type_B = z$Type_S[i],
          Cmo_S = z$Cmo_S[j], GeneName_S = z$GeneName_S[j], Module_S = z$Module_S[j], Type_S = z$Type_S[j],
          ColorRole_B = "S", ColorRole_S = "S", PairDisplay = "SS-paired",
          stringsAsFactors = FALSE
        )
      }
    }))
  })
  bind_rows(out)
}

both <- inner_join(b, s, by = "GeneName_base", relationship = "many-to-many") %>%
  mutate(
    ColorRole_B = "B", ColorRole_S = "S",
    PairDisplay = ifelse(Type_B == "BS" & Type_S == "BS", "BS", "non-BS-paired")
  )
bb_source <- b %>% filter(Type_B == "BB")
ss_source <- s %>% filter(Type_S == "SS")
bb_pairs <- make_same_side_pairs(bb_source, "B")
ss_pairs <- make_same_side_pairs(ss_source, "S")
bb_used <- if (is.null(bb_pairs)) character() else unique(c(bb_pairs$Cmo_B, bb_pairs$Cmo_S))
ss_used <- if (is.null(ss_pairs)) character() else unique(c(ss_pairs$Cmo_B, ss_pairs$Cmo_S))
b_only <- anti_join(b, s, by = "GeneName_base") %>%
  filter(!Cmo_B %in% bb_used) %>%
  mutate(Cmo_S = NA_character_, GeneName_S = NA_character_, Module_S = NA_character_, Type_S = NA_character_,
         ColorRole_B = "B", ColorRole_S = NA_character_, PairDisplay = "B-only")
s_only <- anti_join(s, b, by = "GeneName_base") %>%
  filter(!Cmo_S %in% ss_used) %>%
  mutate(Cmo_B = NA_character_, GeneName_B = NA_character_, Module_B = NA_character_, Type_B = NA_character_,
         ColorRole_B = NA_character_, ColorRole_S = "S", PairDisplay = "S-only")
rows <- bind_rows(both, bb_pairs, ss_pairs, b_only, s_only) %>%
  mutate(
    PairID = paste(GeneName_base, ifelse(is.na(Cmo_B), "NA", Cmo_B), ifelse(is.na(Cmo_S), "NA", Cmo_S), sep = "|"),
    Module = coalesce(Module_B, Module_S),
    Type_display = PairDisplay
  )
if (anyDuplicated(rows$PairID)) stop("Fig.S6 display PairID is not unique.")
missing_expression_genes <- setdiff(unique(c(rows$Cmo_B, rows$Cmo_S)), expr$GID)
missing_expression_genes <- missing_expression_genes[!is.na(missing_expression_genes)]
if (length(missing_expression_genes)) {
  missing_from_counts <- setdiff(missing_expression_genes, counts_ref$GID)
  if (length(missing_from_counts)) {
    stop("Fig.S6 gene(s) absent from both normalized expression and count reference: ",
         paste(head(missing_from_counts, 20), collapse = ", "))
  }
  missing_idx <- match(missing_expression_genes, counts_ref$GID)
  missing_count_values <- counts_ref_mat[missing_idx, , drop = FALSE]
  if (any(missing_count_values != 0)) {
    stop("Non-zero Fig.S6 gene(s) are absent from normalized expression: ",
         paste(head(missing_expression_genes[rowSums(missing_count_values != 0) > 0], 20), collapse = ", "))
  }
  zero_rows <- as.data.frame(matrix(0, nrow = length(missing_expression_genes), ncol = length(samples$sample)),
                             check.names = FALSE)
  names(zero_rows) <- samples$sample
  zero_rows <- cbind(GID = missing_expression_genes, zero_rows)
  expr <- bind_rows(expr, zero_rows)
  expr_mat <- as.matrix(expr[, samples$sample, drop = FALSE])
  storage.mode(expr_mat) <- "double"
}

timepoints <- unique(as.character(samples$timepoint)); timepoints <- timepoints[order(suppressWarnings(as.numeric(timepoints)))]
mean_side <- function(ids, side) {
  out <- matrix(NA_real_, nrow = length(ids), ncol = length(timepoints), dimnames = list(NULL, paste0(side, "_", timepoints)))
  for (j in seq_along(timepoints)) {
    smp <- samples$sample[samples$timepoint == timepoints[j]]
    idx <- match(ids, expr$GID)
    vals <- expr_mat[idx, smp, drop = FALSE]
    out[, j] <- ifelse(is.na(ids), NA_real_, rowMeans(vals, na.rm = TRUE))
  }
  out
}
mean_B <- mean_side(rows$Cmo_B, "B"); mean_S <- mean_side(rows$Cmo_S, "S")
mean_tbl <- bind_cols(rows %>% select(PairID, GeneName_base, Cmo_B, GeneName_B, Cmo_S, GeneName_S, Module, Type_display, ColorRole_B, ColorRole_S),
                      as.data.frame(mean_B), as.data.frame(mean_S))

s6 <- read_tsv(table_s6_path, show_col_types = FALSE)
if (!"PairID" %in% names(s6) || anyDuplicated(s6$PairID)) stop("Table S6 lacks a unique PairID column.")
s6_pair_ids <- s6$PairID
if (any(!s6_pair_ids %in% rows$PairID)) stop("Table S6 contains PairID not present in Fig.S6 rows.")
bs_rows <- rows %>% filter(Type_display == "BS")
if (setequal(s6_pair_ids, bs_rows$PairID) == FALSE) stop("True BS Fig.S6 rows and Table S6 PairID set are not identical.")

star_cols <- intersect(grep("^star_(00|24|48)$", names(s6), value = TRUE), names(s6))
if (!all(paste0("star_", timepoints) %in% star_cols)) stop("Table S6 lacks star columns for all actual timepoints.")
star_tbl <- s6 %>% select(PairID,
                          all_of(paste0("star_", timepoints)),
                          all_of(paste0("direction_", timepoints)))
mapped <- mean_tbl %>% left_join(star_tbl, by = "PairID")
if (nrow(mapped) != nrow(mean_tbl) || any(is.na(mapped$PairID))) stop("Fig.S6/Table S6 row mapping failed.")
write_tsv(mapped, file.path(out_dir, "FigS6_expression_means_NEW.tsv"))
write_tsv(tibble(
  n_heatmap_rows = nrow(rows), n_true_BS_rows = nrow(bs_rows), n_tableS6_rows = nrow(s6),
  unmatched_heatmap_rows = sum(rows$Type_display == "BS" & !rows$PairID %in% s6_pair_ids),
  unmatched_tableS6_rows = sum(!s6_pair_ids %in% rows$PairID), duplicated_tableS6_PairID = anyDuplicated(s6$PairID) > 0
), file.path(out_dir, "FigS6_TableS6_star_mapping_QC.tsv"))

# Preserve annotation order within modules and Type groups; do not cluster the heatmaps.
module_order <- unique(rows$Module[!is.na(rows$Module) & rows$Module != ""])
type_order <- unique(rows$Type_display[!is.na(rows$Type_display)])
rows$module_factor <- factor(rows$Module, levels = module_order, ordered = TRUE)
rows$type_factor <- factor(rows$Type_display, levels = type_order, ordered = TRUE)
ord <- order(rows$module_factor, rows$type_factor, seq_len(nrow(rows)), na.last = TRUE)
rows <- rows[ord, , drop = FALSE]; mapped <- mapped[match(rows$PairID, mapped$PairID), , drop = FALSE]

all_B <- log10(as.matrix(mapped[, paste0("B_", timepoints), drop = FALSE]) + 1)
all_S <- log10(as.matrix(mapped[, paste0("S_", timepoints), drop = FALSE]) + 1)
rownames(all_B) <- rows$PairID; rownames(all_S) <- rows$PairID

# Derive one numerical range from all log10(mean TPM + 1) values for both B/S palettes.
transformed_values <- c(all_B[is.finite(all_B)], all_S[is.finite(all_S)])
if (!length(transformed_values)) stop("No finite transformed corrected TPM values available for Fig.S6 scale.")
global_heatmap_limits <- range(transformed_values)
global_heatmap_breaks <- pretty(global_heatmap_limits, n = 5)

global_heatmap_breaks <- global_heatmap_breaks[
  global_heatmap_breaks >= global_heatmap_limits[1] &
    global_heatmap_breaks <= global_heatmap_limits[2]
]
if (length(global_heatmap_breaks) < 2L) global_heatmap_breaks <- global_heatmap_limits
global_color_breaks <- seq(global_heatmap_limits[1], global_heatmap_limits[2], length.out = 50)
global_blue_palette <- c("white", "#93C6DF", "#1B66AE", "#00008B")
global_red_palette <- c("white", "#FDDCC8", "#D94021", "#B31127")
global_blue_colors <- colorRampPalette(global_blue_palette)(length(global_color_breaks))
global_red_colors <- colorRampPalette(global_red_palette)(length(global_color_breaks))
global_na_color <- "grey"
global_oob <- "squish"

clamp_global <- function(x) {
  pmin(pmax(x, global_heatmap_limits[1]), global_heatmap_limits[2])
}
col_B_fun <- colorRamp2(global_color_breaks, global_blue_colors)
col_S_fun <- colorRamp2(global_color_breaks, global_red_colors)
col_B <- function(x) col_B_fun(clamp_global(x))
col_S <- function(x) col_S_fun(clamp_global(x))

global_legend_pdf <- file.path(out_dir, "FigS6_global_heatmap_legend_NEW.pdf")
grDevices::cairo_pdf(global_legend_pdf, width = 5, height = 2.2,
                     bg = "transparent", family = "Arial")
lg_B <- Legend(col_fun = col_B, title = "CmB",
               direction = "horizontal",
               at = global_heatmap_breaks, labels = as.character(global_heatmap_breaks),
               title_gp = gpar(fontfamily = "Arial", fontsize = 14),
               grid_width = unit(60, "mm"),
               grid_height = unit(8, "mm"),
               labels_gp = gpar(fontfamily = "Arial", fontsize = 12))
lg_S <- Legend(col_fun = col_S, title = "CmS",
               direction = "horizontal",
               at = global_heatmap_breaks, labels = as.character(global_heatmap_breaks),
               title_gp = gpar(fontfamily = "Arial", fontsize = 14),
               grid_width = unit(60, "mm"),
               grid_height = unit(8, "mm"),
               labels_gp = gpar(fontfamily = "Arial", fontsize = 12))
draw(packLegend(lg_B, lg_S, direction = "horizontal", gap = unit(8, "mm")))
dev.off()

safe_name <- function(x) gsub("[^[:alnum:]_.-]+", "_", x)
for (m in module_order) {
  sel <- !is.na(rows$Module) & rows$Module == m
  if (!any(sel)) next
  row_ids <- rows$PairID[sel]
  mat_B <- all_B[row_ids, , drop = FALSE]
  mat_S <- all_S[row_ids, , drop = FALSE]
  lab_B <- ifelse(is.na(rows$GeneName_B[sel]), "", rows$GeneName_B[sel])
  lab_S <- ifelse(is.na(rows$GeneName_S[sel]), "", rows$GeneName_S[sel])
  color_role_B <- as.character(rows$ColorRole_B[sel])
  color_role_S <- as.character(rows$ColorRole_S[sel])
  row_split <- factor(rows$Type_display[sel], levels = type_order, ordered = TRUE)
  star_B_matrix <- matrix("", nrow = sum(sel), ncol = length(timepoints),
                          dimnames = list(row_ids, paste0("star_B_", timepoints)))
  star_S_matrix <- matrix("", nrow = sum(sel), ncol = length(timepoints),
                          dimnames = list(row_ids, paste0("star_S_", timepoints)))
  for (j in seq_along(timepoints)) {
    star_flag <- mapped[[paste0("star_", timepoints[j])]][sel]
    direction <- mapped[[paste0("direction_", timepoints[j])]][sel]
    star_B_matrix[, j] <- ifelse(!is.na(star_flag) & !is.na(direction) &
                                   star_flag == "*" & direction == "B>S", "*", "")
    star_S_matrix[, j] <- ifelse(!is.na(star_flag) & !is.na(direction) &
                                   star_flag == "*" & direction == "S>B", "*", "")
  }
  cell_B <- function(j, i, x, y, width, height, fill) {
    if (!is.na(mat_B[i, j])) {
      fill_col <- if (identical(color_role_B[i], "S")) col_S(mat_B[i, j]) else col_B(mat_B[i, j])
      grid.rect(x=x, y=y, width=width, height=height,
                gp=gpar(fill=fill_col, col=NA))
    }
    if (!is.na(mat_B[i, j]) && mat_B[i, j] < 1e-5) {
      grid.segments(x-width*.4, y-height*.4, x+width*.4, y+height*.4, gp=gpar(col="grey50", lwd=1))
      grid.segments(x-width*.4, y+height*.4, x+width*.4, y-height*.4, gp=gpar(col="grey50", lwd=1))
    }
    if (star_B_matrix[i, j] == "*") grid.text("*", x, y - unit(1, "mm"), gp=gpar(fontsize=25, col="black"))
  }
  cell_S <- function(j, i, x, y, width, height, fill) {
    if (!is.na(mat_S[i, j])) {
      fill_col <- if (identical(color_role_S[i], "B")) col_B(mat_S[i, j]) else col_S(mat_S[i, j])
      grid.rect(x=x, y=y, width=width, height=height,
                gp=gpar(fill=fill_col, col=NA))
    }
    if (!is.na(mat_S[i, j]) && mat_S[i, j] < 1e-5) {
      grid.segments(x-width*.4, y-height*.4, x+width*.4, y+height*.4, gp=gpar(col="grey50", lwd=1))
      grid.segments(x-width*.4, y+height*.4, x+width*.4, y-height*.4, gp=gpar(col="grey50", lwd=1))
    }
    if (star_S_matrix[i, j] == "*") grid.text("*", x, y - unit(1, "mm"), gp=gpar(fontsize=25, col="black"))
  }
  ht_B <- Heatmap(mat_B, name="CmB", col=col_B, width=unit(10*ncol(mat_B), "mm"), height=unit(10*nrow(mat_B), "mm"),
                  cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=TRUE, row_labels=lab_B,
                  row_names_side="left", row_names_gp=gpar(fontsize=22, fontface="italic", fontfamily="Arial"),
                  row_names_max_width=unit(60, "mm"), row_split=row_split, row_title=NULL, show_column_names=FALSE,
                  cell_fun=cell_B, na_col=global_na_color, show_heatmap_legend=FALSE)
  ht_S <- Heatmap(mat_S, name="CmS", col=col_S, width=unit(10*ncol(mat_S), "mm"), height=unit(10*nrow(mat_S), "mm"),
                  cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=TRUE, row_labels=lab_S,
                  row_names_side="right", row_names_gp=gpar(fontsize=22, fontface="italic", fontfamily="Arial"),
                  row_names_max_width=unit(60, "mm"), row_split=row_split, row_title=NULL, show_column_names=FALSE,
                  cell_fun=cell_S, na_col=global_na_color, show_heatmap_legend=FALSE)
  grDevices::cairo_pdf(file.path(out_dir, paste0("Heatmap_Module_", safe_name(m), "_NEW.pdf")),
                       width=9.5, height=max(8, 0.5*nrow(mat_B)), bg="transparent", family="Arial")
  draw(ht_B + ht_S, show_heatmap_legend=FALSE, show_annotation_legend=FALSE)
  dev.off()
}

message("Supplementary Fig. 6 heatmaps written.")
}))
