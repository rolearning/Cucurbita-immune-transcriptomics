#!/usr/bin/env Rscript

# Import Salmon estimates using sample metadata, GFF annotations and split-gene corrections.
# Export DDS, contrasts, VST, scaledTPM and QC under results/01_tximport/.
# Merge split-gene counts and anchor corrected lengths to curated CDS lengths.
# Fit DESeq2 to corrected counts with tximport offsets; export scaledTPM separately.

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(rtracklayer)
  library(S4Vectors)
  library(readr)
  library(dplyr)
  library(tibble)
})

project_root <- RELEASE_ROOT
sample_sheet <- "0.salmon/sample_info.tsv"
out_dir      <- "1.normalization/results/01_tximport"
gene_correction_file <- "1.normalization/gene_correction_with_edited_cds_length.tsv"

# For tx2gene input, the first two columns must be transcript and gene IDs.
species_config <- tribble(
  ~species, ~annotation_file,          ~annotation_format,
  "Cm",     "ref/Cm.gff3",    "gff3",
  "Ls",     "ref/Ls.gff3",    "gff3",
  "Ta",     "ref/Ta.gff3",    "gff3"
)

# Match transcript-ID handling to the Salmon index and inspect mapping QC.
ignore_tx_version <- FALSE
ignore_after_bar  <- FALSE
min_tx_match_fraction <- 0.90

reference_timepoint <- "00"
use_batch_if_present <- TRUE
within_sizefactor_type <- "ratio"

# This prefilter applies only to DESeq2/VST, not scaledTPM exports.
deseq_min_count   <- 10
deseq_min_samples <- 2
export_deg_contrasts <- TRUE
deg_alpha <- 0.05

resolve_path <- function(path, root = project_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(root, path)
}

stop_if_missing <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(label, " not found:\n", paste0("  - ", missing, collapse = "\n"))
  }
}

write_matrix_tsv <- function(mat, path, id_col = "GID") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_tsv(rownames_to_column(as.data.frame(mat), id_col), path)
}

normalize_tx_ids <- function(x) {
  x <- as.character(x)
  if (ignore_after_bar)  x <- sub("\\|.*$", "", x)
  if (ignore_tx_version) x <- sub("\\.[0-9]+$", "", x)
  x
}

first_existing_column <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

to_character_scalar <- function(x) {
  vapply(seq_along(x), function(i) {
    value <- as.character(x[[i]])
    if (length(value) == 0) NA_character_ else value[[1]]
  }, character(1))
}

read_tx2gene <- function(annotation_file, annotation_format) {
  annotation_format <- tolower(annotation_format)

  if (annotation_format == "tx2gene") {
    x <- read_tsv(
      annotation_file,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    )
    if (ncol(x) < 2) stop("tx2gene file must contain at least two columns: ", annotation_file)

    tx_col <- first_existing_column(x, c("tx_id", "transcript_id", "TXNAME", "transcript"))
    gene_col <- first_existing_column(x, c("gene_id", "GENEID", "gene"))
    if (is.null(tx_col)) tx_col <- names(x)[1]
    if (is.null(gene_col)) gene_col <- names(x)[2]

    ans <- transmute(x, tx_id = .data[[tx_col]], gene_id = .data[[gene_col]])

  } else if (annotation_format == "gff3") {
    gff <- rtracklayer::import(annotation_file)
    gff_df <- as.data.frame(gff)

    type_col <- first_existing_column(gff_df, c("type"))
    if (is.null(type_col)) stop("Cannot find feature type in GFF: ", annotation_file)

    is_tx <- as.character(gff_df[[type_col]]) %in% c("mRNA", "transcript")
    tx_df <- gff_df[is_tx, , drop = FALSE]
    if (nrow(tx_df) == 0) {
      stop("No mRNA/transcript features found in GFF: ", annotation_file)
    }

    tx_col <- first_existing_column(tx_df, c("transcript_id", "ID", "Name"))
    gene_col <- first_existing_column(tx_df, c("gene_id", "Parent", "gene", "gene_name"))
    if (is.null(tx_col) || is.null(gene_col)) {
      stop(
        "Cannot identify transcript/gene attributes in GFF: ", annotation_file,
        "\nAvailable columns: ", paste(names(tx_df), collapse = ", ")
      )
    }

    ans <- tibble(
      tx_id   = to_character_scalar(tx_df[[tx_col]]),
      gene_id = to_character_scalar(tx_df[[gene_col]])
    )

  } else {
    stop("annotation_format must be 'gff3' or 'tx2gene'; got: ", annotation_format)
  }

  ans <- ans %>%
    filter(!is.na(tx_id), tx_id != "", !is.na(gene_id), gene_id != "") %>%
    distinct(tx_id, gene_id)

  # Each transcript must map to one gene for tximport.
  ambiguous <- ans %>% dplyr::count(tx_id, name = "n_gene") %>% filter(n_gene > 1)
  if (nrow(ambiguous) > 0) {
    stop(
      "A transcript maps to more than one gene in ", annotation_file,
      ". First IDs: ", paste(head(ambiguous$tx_id, 10), collapse = ", ")
    )
  }

  ans
}

sort_timepoints <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  if (all(is.finite(numeric_x))) x[order(numeric_x)] else sort(x)
}
estimate_scaled_size_factors <- function(count_mat, type = "ratio") {
  keep <- rowSums(count_mat) > 0
  if (!any(keep)) stop("The scaledTPM matrix contains no non-zero genes.")

  sf <- DESeq2::estimateSizeFactorsForMatrix(count_mat[keep, , drop = FALSE], type = type)
  if (any(!is.finite(sf)) || any(sf <= 0)) {
    stop("Invalid size factors estimated from scaledTPM counts.")
  }
  sf
}
estimate_tpm_rate_scale <- function(txi, excluded_genes = character()) {
  usable <- setdiff(rownames(txi$counts), excluded_genes)
  if (length(usable) < 100) {
    stop("Too few unaffected genes available to estimate TPM/rate scale.")
  }

  rate <- txi$counts[usable, , drop = FALSE] /
    txi$length[usable, , drop = FALSE]

  abund <- txi$abundance[usable, , drop = FALSE]

  scale <- vapply(seq_len(ncol(rate)), function(j) {
    ok <- is.finite(rate[, j]) & rate[, j] > 0 &
      is.finite(abund[, j]) & abund[, j] > 0

    if (sum(ok) < 100) {
      stop(
        "Too few positive unaffected genes to estimate TPM/rate scale for sample ",
        colnames(rate)[[j]], "."
      )
    }

    median(abund[ok, j] / rate[ok, j])

  }, numeric(1))

  names(scale) <- colnames(rate)
  scale
}

apply_split_gene_corrections <- function(txi, correction_tbl, species, audit_path) {
  if (nrow(correction_tbl) == 0) {
    message("No split-gene corrections for species ", species, ".")
    return(txi)
  }

  if (!all(c("counts", "abundance", "length") %in% names(txi))) {
    stop("tximport object does not contain counts/abundance/length.")
  }

  if (!identical(rownames(txi$counts), rownames(txi$length)) ||
      !identical(rownames(txi$counts), rownames(txi$abundance))) {
    stop("tximport counts/abundance/length row names are not aligned for ", species, ".")
  }

  merge_sets <- lapply(correction_tbl$merge_genes, function(x) {
    y <- strsplit(x, ";", fixed = TRUE)[[1]]
    trimws(y[nzchar(trimws(y))])
  })
  all_corrected_ids <- unique(unlist(merge_sets))

  tpm_rate_scale <- estimate_tpm_rate_scale(
    txi,
    excluded_genes = all_corrected_ids
  )
  audit_list <- vector("list", nrow(correction_tbl))

  for (i in seq_len(nrow(correction_tbl))) {
    keep_gene <- correction_tbl$keep_gene[[i]]
    merge_genes <- merge_sets[[i]]

    edited_cds_length <- suppressWarnings(
      as.numeric(correction_tbl$edited_cds_length[[i]])
    )

    if (!is.finite(edited_cds_length) || edited_cds_length <= 0) {
      stop(
        "Invalid edited_cds_length for ", species, " / ", keep_gene, ": ",
        correction_tbl$edited_cds_length[[i]]
      )
    }

    if (!keep_gene %in% merge_genes) {
      stop("keep_gene is not included in merge_genes: ", keep_gene)
    }

    missing_genes <- setdiff(merge_genes, rownames(txi$counts))
    if (length(missing_genes) > 0) {
      stop(
        "Split-gene correction genes missing from tximport matrix for ", species,
        " / ", keep_gene, ": ", paste(missing_genes, collapse = ", ")
      )
    }

    counts_before <- txi$counts[merge_genes, , drop = FALSE]
    merged_counts <- colSums(counts_before)
    # Preserve relative sample lengths while anchoring their mean to the edited CDS.
    fragment_length_sum <- colSums(
      txi$length[merge_genes, , drop = FALSE]
    )

    if (any(!is.finite(fragment_length_sum)) ||
        any(fragment_length_sum <= 0)) {
      stop(
        "Invalid fragment length(s) for ",
        species, " / ", keep_gene, "."
      )
    }

    relative_length_pattern <-
      fragment_length_sum / mean(fragment_length_sum)

    corrected_length <-
      edited_cds_length * relative_length_pattern

    corrected_rate <- merged_counts / corrected_length

    corrected_abundance <-
      corrected_rate *
      tpm_rate_scale[colnames(txi$counts)]

    txi$counts[keep_gene, ] <- merged_counts
    txi$length[keep_gene, ] <- corrected_length
    txi$abundance[keep_gene, ] <- corrected_abundance

    remove_genes <- setdiff(merge_genes, keep_gene)
    if (length(remove_genes) > 0) {
      keep_rows <- setdiff(rownames(txi$counts), remove_genes)
      txi$counts <- txi$counts[keep_rows, , drop = FALSE]
      txi$length <- txi$length[keep_rows, , drop = FALSE]
      txi$abundance <- txi$abundance[keep_rows, , drop = FALSE]
    }

    audit_list[[i]] <- tibble(
      species = species,
      keep_gene = keep_gene,
      merge_genes = paste(merge_genes, collapse = ";"),
      sample = colnames(txi$counts),
      merged_raw_count = as.numeric(merged_counts),
      fragment_length_sum_before = as.numeric(fragment_length_sum),
      corrected_length_used = as.numeric(corrected_length),
      edited_cds_length_target = edited_cds_length,
      corrected_abundance_before_renorm =
        as.numeric(corrected_abundance)
    )

    message(
      "Applied split-gene correction: ", species, " / ",
      paste(merge_genes, collapse = " + "), " -> ", keep_gene,
      " ; edited CDS length = ", edited_cds_length, " bp"
    )
  }
  # Restore each TPM column total to 1e6 after split-gene replacement.
  abundance_totals <- colSums(txi$abundance)

  if (any(!is.finite(abundance_totals)) ||
      any(abundance_totals <= 0)) {
    stop(
      "Invalid abundance totals after split-gene correction for ",
      species, "."
    )
  }

  txi$abundance <- sweep(
    txi$abundance,
    2,
    abundance_totals / 1e6,
    "/"
  )
  audit_tbl <- bind_rows(audit_list)
  write_tsv(audit_tbl, audit_path)

  txi
}

make_scaled_txi <- function(txi) {
  out <- txi
  out$counts <- tximport::makeCountsFromAbundance(
    countsMat = txi$counts,
    abundanceMat = txi$abundance,
    lengthMat = txi$length,
    countsFromAbundance = "scaledTPM"
  )
  out$countsFromAbundance <- "scaledTPM"
  out
}

sample_sheet_path <- resolve_path(sample_sheet)
stop_if_missing(sample_sheet_path, "Sample sheet")

samples <- read_tsv(
  sample_sheet_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

required_sample_cols <- c("sample", "species", "quant_file", "timepoint", "replicate")
missing_sample_cols <- setdiff(required_sample_cols, names(samples))
if (length(missing_sample_cols) > 0) {
  stop("Sample sheet is missing columns: ", paste(missing_sample_cols, collapse = ", "))
}
if (anyDuplicated(samples$sample)) stop("Duplicated sample IDs in sample sheet.")
if (any(is.na(samples$sample) | samples$sample == "")) stop("Empty sample ID in sample sheet.")

# Require 25 libraries; exclude Cm_48h_D2 and Ls_24h_D1.
excluded_final_samples <- c("Cm_48h_D2", "Ls_24h_D1")
unexpected_excluded <- intersect(samples$sample, excluded_final_samples)
if (length(unexpected_excluded) > 0) {
  stop(
    "Excluded libraries must not appear in the FINAL sample sheet: ",
    paste(unexpected_excluded, collapse = ", ")
  )
}

expected_final_samples <- c(
  paste0("Cm_", c("0h_D1", "0h_D2", "0h_D3", "24h_D1", "24h_D2", "24h_D3", "48h_D1", "48h_D3")),
  paste0("Ls_", c("0h_D1", "0h_D2", "0h_D3", "24h_D2", "24h_D3", "48h_D1", "48h_D2", "48h_D3")),
  paste0("Ta_", c("0h_D1", "0h_D2", "0h_D3", "24h_D1", "24h_D2", "24h_D3", "48h_D1", "48h_D2", "48h_D3"))
)
observed_final_samples <- samples$sample
missing_final_samples <- setdiff(expected_final_samples, observed_final_samples)
unexpected_samples <- setdiff(observed_final_samples, expected_final_samples)
expected_species_by_sample <- setNames(
  sub("_.*$", "", expected_final_samples),
  expected_final_samples
)
species_mismatch <- samples$sample[
  !is.na(samples$sample) &
    samples$sample %in% expected_final_samples &
    samples$species != unname(expected_species_by_sample[samples$sample])
]
if (length(observed_final_samples) != 25L || length(missing_final_samples) > 0L ||
    length(unexpected_samples) > 0L || length(species_mismatch) > 0L ||
    !identical(as.integer(table(factor(samples$species, levels = c("Cm", "Ls", "Ta")))), c(8L, 8L, 9L))) {
  stop(
    "sample_info.tsv must contain exactly the FINAL 25 libraries (Cm=8, Ls=8, Ta=9).\n",
    "Missing: ", if (length(missing_final_samples)) paste(missing_final_samples, collapse = ", ") else "none", "\n",
    "Unexpected: ", if (length(unexpected_samples)) paste(unexpected_samples, collapse = ", ") else "none", "\n",
    "Species-label mismatch: ", if (length(species_mismatch)) paste(species_mismatch, collapse = ", ") else "none", "\n",
    "Observed species counts: ", paste(names(table(samples$species)), table(samples$species), sep = "=", collapse = "; ")
  )
}

samples <- samples %>%
  mutate(quant_file = vapply(quant_file, resolve_path, character(1)))
stop_if_missing(samples$quant_file, "Salmon quant.sf file(s)")

if (anyDuplicated(species_config$species)) stop("Duplicated species in species_config.")
species_config <- species_config %>%
  mutate(annotation_file = vapply(annotation_file, resolve_path, character(1)))
stop_if_missing(species_config$annotation_file, "Annotation file(s)")

missing_species_config <- setdiff(unique(samples$species), species_config$species)
if (length(missing_species_config) > 0) {
  stop("Species missing from species_config: ", paste(missing_species_config, collapse = ", "))
}

gene_correction_path <- resolve_path(gene_correction_file)
stop_if_missing(gene_correction_path, "Gene-correction table")

gene_corrections <- read_tsv(
  gene_correction_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

required_correction_cols <- c(
  "species", "keep_gene", "merge_genes", "edited_cds_length", "status"
)
missing_correction_cols <- setdiff(required_correction_cols, names(gene_corrections))
if (length(missing_correction_cols) > 0) {
  stop(
    "Gene-correction table is missing columns: ",
    paste(missing_correction_cols, collapse = ", ")
  )
}

bad_corrections <- gene_corrections %>%
  filter(is.na(status) | status != "OK")
if (nrow(bad_corrections) > 0) {
  stop(
    "Gene-correction table contains non-OK rows. Resolve them before running 01. ",
    "keep_gene(s): ", paste(bad_corrections$keep_gene, collapse = ", ")
  )
}

if (anyDuplicated(gene_corrections$keep_gene)) {
  stop("Duplicated keep_gene values in gene-correction table.")
}

all_merge_ids <- unlist(strsplit(gene_corrections$merge_genes, ";", fixed = TRUE))
all_merge_ids <- trimws(all_merge_ids)
if (anyDuplicated(all_merge_ids)) {
  dup_ids <- unique(all_merge_ids[duplicated(all_merge_ids)])
  stop(
    "A gene occurs in more than one split-gene correction group: ",
    paste(dup_ids, collapse = ", ")
  )
}

dir.create(resolve_path(out_dir), recursive = TRUE, showWarnings = FALSE)

all_import_qc <- list()
all_sample_qc <- list()

for (sp in species_config$species) {
  message("\n===== Processing species: ", sp, " =====")

  cfg <- species_config %>% filter(species == sp)
  smp <- samples %>% filter(species == sp)
  if (nrow(smp) == 0) {
    warning("No samples found for species ", sp, "; skipping.")
    next
  }

  sp_out <- file.path(resolve_path(out_dir), sp)
  dir.create(sp_out, recursive = TRUE, showWarnings = FALSE)

  tx2gene <- read_tx2gene(cfg$annotation_file[[1]], cfg$annotation_format[[1]])

  # Check transcript IDs in every sample quantified against the same index.
  mapping_qc <- lapply(seq_len(nrow(smp)), function(i) {
    q <- read_tsv(
      smp$quant_file[[i]],
      col_select = "Name",
      col_types = cols(Name = col_character()),
      show_col_types = FALSE
    )
    salmon_ids <- unique(normalize_tx_ids(q$Name))
    map_ids <- unique(normalize_tx_ids(tx2gene$tx_id))
    n_match <- sum(salmon_ids %in% map_ids)
    tibble(
      sample = smp$sample[[i]],
      species = sp,
      salmon_transcripts = length(salmon_ids),
      matched_transcripts = n_match,
      unmatched_transcripts = length(salmon_ids) - n_match,
      match_fraction = n_match / length(salmon_ids)
    )
  }) %>% bind_rows()

  write_tsv(mapping_qc, file.path(sp_out, paste0(sp, "_transcript_mapping_QC.tsv")))
  if (any(mapping_qc$match_fraction < min_tx_match_fraction)) {
    bad <- mapping_qc %>% filter(match_fraction < min_tx_match_fraction)
    stop(
      "Transcript-ID match fraction below ", min_tx_match_fraction,
      " for: ", paste(bad$sample, collapse = ", "),
      ". Check GFF/FASTA versions and ID prefixes before continuing."
    )
  }

  files <- setNames(smp$quant_file, smp$sample)

  txi_original <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance = "no",
    ignoreTxVersion = ignore_tx_version,
    ignoreAfterBar = ignore_after_bar,
    dropInfReps = TRUE
  )

  # Correct split genes before normalization.
  correction_sp <- gene_corrections %>% filter(species == sp)

  txi_original <- apply_split_gene_corrections(
    txi = txi_original,
    correction_tbl = correction_sp,
    species = sp,
    audit_path = file.path(sp_out, paste0(sp, "_split_gene_correction_audit.tsv"))
  )

  # Recompute scaledTPM after correction; summing the old scaledTPM rows is invalid.

  txi_scaled <- make_scaled_txi(txi_original)

  saveRDS(
    txi_original,
    file.path(sp_out, paste0(sp, "_txi_original_counts_offset.rds"))
  )

  saveRDS(
    txi_scaled,
    file.path(sp_out, paste0(sp, "_txi_scaledTPM.rds"))
  )

  write_matrix_tsv(
    txi_scaled$counts,
    file.path(sp_out, paste0(sp, "_scaledTPM_counts.tsv"))
  )

  coldata <- as.data.frame(smp)
  rownames(coldata) <- coldata$sample
  coldata <- coldata[colnames(txi_original$counts), , drop = FALSE]

  tp_levels <- sort_timepoints(coldata$timepoint)
  if (!reference_timepoint %in% tp_levels) {
    stop(
      "reference_timepoint '", reference_timepoint,
      "' is absent for species ", sp,
      ". Available: ", paste(tp_levels, collapse = ", ")
    )
  }
  coldata$timepoint <- factor(coldata$timepoint, levels = tp_levels)
  coldata$timepoint <- relevel(coldata$timepoint, ref = reference_timepoint)
  coldata$replicate <- factor(coldata$replicate)

  has_batch <-
    use_batch_if_present &&
    "batch" %in% names(coldata) &&
    any(!is.na(coldata$batch) & coldata$batch != "") &&
    length(unique(coldata$batch[!is.na(coldata$batch) & coldata$batch != ""])) > 1

  if (has_batch) {
    if (any(is.na(coldata$batch) | coldata$batch == "")) {
      stop("Batch is missing for some, but not all, samples of species ", sp, ".")
    }
    coldata$batch <- factor(coldata$batch)
    design_formula <- ~ batch + timepoint
  } else {
    design_formula <- ~ timepoint
  }

  dds <- DESeqDataSetFromTximport(
    txi_original,
    colData = coldata,
    design = design_formula
  )

  keep_deseq <- rowSums(counts(dds) >= deseq_min_count) >= deseq_min_samples
  if (sum(keep_deseq) < 100) {
    stop("Fewer than 100 genes pass the DESeq2 prefilter for species ", sp, ".")
  }
  dds <- dds[keep_deseq, ]

  dds <- DESeq(dds, sfType = within_sizefactor_type)
  saveRDS(dds, file.path(sp_out, paste0(sp, "_DESeq2_dds.rds")))

  deseq_norm <- counts(dds, normalized = TRUE)
  write_matrix_tsv(
    deseq_norm,
    file.path(sp_out, paste0(sp, "_DESeq2_normalized_counts_with_offset.tsv"))
  )

  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)
  write_matrix_tsv(vst_mat, file.path(sp_out, paste0(sp, "_VST_blindFALSE.tsv")))
  saveRDS(vsd, file.path(sp_out, paste0(sp, "_VST_blindFALSE.rds")))

  # Retain fractional scaledTPM values after sample-wise normalization.
  ls_counts <- txi_scaled$counts
  ls_counts <- ls_counts[rowSums(ls_counts) > 0, , drop = FALSE]
  ls_sf <- estimate_scaled_size_factors(ls_counts, within_sizefactor_type)
  names(ls_sf) <- colnames(ls_counts)
  ls_norm <- sweep(ls_counts, 2, ls_sf, "/")

  write_matrix_tsv(
    ls_norm,
    file.path(sp_out, paste0(sp, "_scaledTPM_withinSpecies_normalized.tsv"))
  )

if (export_deg_contrasts) {

  deg_contrasts <- list(
    "24_vs_00" = c("24", "00"),
    "48_vs_00" = c("48", "00"),
    "48_vs_24" = c("48", "24")
  )

  for (contrast_name in names(deg_contrasts)) {

    numerator <- deg_contrasts[[contrast_name]][1]
    denominator <- deg_contrasts[[contrast_name]][2]

    res <- results(
      dds,
      contrast = c("timepoint", numerator, denominator),
      alpha = deg_alpha
    )

    res_tbl <- as.data.frame(res) %>%
      rownames_to_column("GID") %>%
      arrange(padj)

    write_tsv(
      res_tbl,
      file.path(sp_out, paste0(sp, "_DEG_", contrast_name, ".tsv"))
    )
  }
}

  dds_sf <- sizeFactors(dds)
  if (is.null(dds_sf)) dds_sf <- rep(NA_real_, ncol(dds))

  dds_nf <- normalizationFactors(dds)
  if (is.null(dds_nf)) {
    dds_nf_col_gmean <- rep(NA_real_, ncol(dds))
  } else {
    dds_nf_col_gmean <- apply(dds_nf, 2, function(x) exp(mean(log(x[x > 0]))))
  }

  sample_qc <- tibble(
    sample = colnames(txi_original$counts),
    species = sp,
    original_estimated_library = colSums(txi_original$counts),
    scaledTPM_library = colSums(txi_scaled$counts),
    DESeq2_size_factor = as.numeric(dds_sf),
    DESeq2_normalization_factor_col_gmean = as.numeric(dds_nf_col_gmean),
    scaledTPM_size_factor = as.numeric(ls_sf[colnames(txi_original$counts)]),
    n_genes_original = nrow(txi_original$counts),
    n_genes_DESeq2_prefiltered = nrow(dds),
    n_genes_scaledTPM_nonzero = nrow(ls_counts)
  )
  write_tsv(sample_qc, file.path(sp_out, paste0(sp, "_sample_normalization_QC.tsv")))

  all_import_qc[[sp]] <- mapping_qc
  all_sample_qc[[sp]] <- sample_qc
}

write_tsv(
  bind_rows(all_import_qc),
  file.path(resolve_path(out_dir), "all_species_transcript_mapping_QC.tsv")
)
write_tsv(
  bind_rows(all_sample_qc),
  file.path(resolve_path(out_dir), "all_species_sample_normalization_QC.tsv")
)

writeLines(capture.output(sessionInfo()), file.path(resolve_path(out_dir), "sessionInfo.txt"))

message("\nDone. Main output directory: ", resolve_path(out_dir))
message("Use *_DESeq2_dds.rds / *_VST_blindFALSE.tsv for within-species analyses.")
message("Use *_scaledTPM_counts.tsv as input to 02_cross_species_scaling.R.")
