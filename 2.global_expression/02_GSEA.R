#!/usr/bin/env Rscript

# Run Fig. 2b GO BP GSEA from per-species DESeq2 DDS and eggNOG annotations.
# Rank genes by Wald statistics for 24 vs 0 h and 48 vs 24 h.
# The DDS uses Salmon counts with length offsets rather than scaledTPM counts.
# Write full fgsea results and separate main and supplementary heatmaps.

.run_release_stage <- function(code) {
  prior_wd <- getwd()
  on.exit(setwd(prior_wd), add = TRUE)
  eval(code, envir = new.env(parent = globalenv()))
}

.run_release_stage(quote({
RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

library(DESeq2)
library(fgsea)
library(tidyverse)
library(stringr)
library(AnnotationDbi)
library(GO.db)
library(ggplot2)
library(patchwork)

project_root <- RELEASE_ROOT

deseq_dir <- file.path(
  project_root,
  "1.normalization/results/01_tximport"
)

gsea_dir <- file.path(
  project_root,
  "2.global_expression/Fig2B_GSEA/results"
)

eggnog_dir <- file.path(
  project_root,
  "eggnog_Cm_Ls_Ta"
)

dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
setwd(gsea_dir)

eggnog_files <- c(
  Cm = file.path(eggnog_dir, "Cmos_eggnog_anno.emapper.annotations"),
  Ls = file.path(eggnog_dir, "Lsic_eggnog_anno.emapper.annotations"),
  Ta = file.path(eggnog_dir, "Tang_eggnog_anno.emapper.annotations")
)

dds_files <- c(
  Cm = file.path(deseq_dir, "Cm/Cm_DESeq2_dds.rds"),
  Ls = file.path(deseq_dir, "Ls/Ls_DESeq2_dds.rds"),
  Ta = file.path(deseq_dir, "Ta/Ta_DESeq2_dds.rds")
)

stopifnot(all(file.exists(eggnog_files)))
stopifnot(all(file.exists(dds_files)))

# Positive statistics indicate higher expression at the later time point.
extract_deseq_result <- function(dds, numerator, denominator) {
  res <- DESeq2::results(
    dds,
    contrast = c("timepoint", numerator, denominator),
    alpha = 0.05
  )

  df <- as.data.frame(res)
  df$id <- rownames(df)

  df %>%
    dplyr::select(
      id,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    )
}

dds_list <- lapply(dds_files, readRDS)
results <- list()
for (sp in c("Cm", "Ls", "Ta")) {

  dds <- dds_list[[sp]]

  results[[paste0(sp, "_24vs00")]] <-
    extract_deseq_result(
      dds,
      numerator = "24",
      denominator = "00"
    )

  results[[paste0(sp, "_48vs24")]] <-
    extract_deseq_result(
      dds,
      numerator = "48",
      denominator = "24"
    )
}
# Remove terminal transcript/protein version suffixes to match gene IDs.
clean_id <- function(x) {
  sub("\\.[0-9]+$", "", x)
}
read_eggnog_gene2go <- function(file) {
  annotation_lines <- readLines(file, warn = FALSE)
  header_line <- grep("^#query\\t", annotation_lines)[1]
  stopifnot(length(header_line) == 1, !is.na(header_line))
  column_names <- strsplit(sub("^#", "", annotation_lines[header_line]), "\\t")[[1]]

  annotation_df <- read.delim(
    file,
    sep = "\t",
    header = FALSE,
    skip = header_line,
    comment.char = "#",
    quote = "",
    fill = TRUE,
    col.names = column_names,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot(all(c("query", "GOs") %in% colnames(annotation_df)))

  annotation_df %>%
    transmute(
      gene_id = clean_id(query),
      GO = GOs
    ) %>%
    filter(!is.na(gene_id), !is.na(GO), GO != "-", GO != "") %>%
    tidyr::separate_rows(GO, sep = ",") %>%
    filter(grepl("^GO:[0-9]{7}$", GO)) %>%
    distinct(gene_id, GO)
}

keep_go_bp <- function(gene2go) {
  go_info <- AnnotationDbi::select(
    GO.db,
    keys = unique(gene2go$GO),
    keytype = "GOID",
    columns = c("ONTOLOGY", "TERM")
  ) %>%
    as_tibble() %>%
    filter(!is.na(GOID), !is.na(ONTOLOGY), ONTOLOGY == "BP") %>%
    distinct(GOID, .keep_all = TRUE)

  gene2go %>%
    inner_join(go_info, by = c("GO" = "GOID")) %>%
    distinct(gene_id, GO, .keep_all = TRUE)
}

add_go_term <- function(fgsea_res) {
  term_map <- AnnotationDbi::select(
    GO.db,
    keys = fgsea_res$pathway,
    keytype = "GOID",
    columns = c("TERM")
  ) %>%
    as_tibble() %>%
    distinct(GOID, TERM)

  fgsea_res %>%
    left_join(term_map, by = c("pathway" = "GOID"))
}

make_go_sets <- function(gene2go_bp) {
  go_list <- gene2go_bp %>%
    group_by(GO) %>%
    summarise(genes = list(unique(gene_id)), .groups = "drop") %>%
    deframe()

  go_list
}

make_ranks_from_deseq <- function(res_df, stat_col = "stat",
                                  agg = c("max_abs", "max", "mean")) {
  agg <- match.arg(agg)

  x <- res_df %>%
    mutate(gene_id = clean_id(id)) %>%
    filter(!is.na(gene_id), !is.na(.data[[stat_col]]), is.finite(.data[[stat_col]])) %>%
    transmute(gene_id, stat = .data[[stat_col]])

  # Resolve duplicate gene IDs by the maximum absolute Wald statistic.
  x2 <- x %>%
    group_by(gene_id) %>%
    summarise(
      stat = dplyr::case_when(
        agg == "max_abs" ~ stat[which.max(abs(stat))],
        agg == "max"     ~ max(stat),
        agg == "mean"    ~ mean(stat),
        TRUE             ~ stat[which.max(abs(stat))]
      ),
      .groups = "drop"
    )

  ranks <- x2$stat
  names(ranks) <- x2$gene_id
  ranks <- ranks[is.finite(ranks)]
  sort(ranks, decreasing = TRUE)
}

# Use BP sets of size 10-500 and nPermSimple=10000; no RNG seed is set here.
run_fgsea_bp <- function(res_df, go_sets,
                         minSize = 10, maxSize = 500, nPermSimple = 10000,
                         stat_col = "stat", agg = "max_abs") {
  ranks <- make_ranks_from_deseq(res_df, stat_col = stat_col, agg = agg)

  fg <- fgsea(
    pathways = go_sets,
    stats    = ranks,
    minSize  = minSize,
    nPermSimple = nPermSimple,
    maxSize  = maxSize
  ) %>%
    as_tibble() %>%
    arrange(padj)

  fg
}

plot_fgsea <- function(fg, padj_cut = 0.05, top_n = 10, title = "") {
  d <- fg %>%
    filter(padj < padj_cut) %>%
    mutate(label = ifelse(!is.na(TERM), TERM, pathway))

  up <- d %>% filter(NES > 0) %>% arrange(padj) %>% slice_head(n = top_n)
  dn <- d %>% filter(NES < 0) %>% arrange(padj) %>% slice_head(n = top_n)

  dd <- bind_rows(up, dn) %>%
    mutate(label = forcats::fct_reorder(label, NES))

  ggplot(dd, aes(x = label, y = NES)) +
    geom_col() +
    coord_flip() +
    theme_classic() +
    labs(x = NULL, y = "NES", title = title)
}

gsea_barplots <- fglist <- list()
annotation_coverage <- list()
for(s in 1:3){
  if(s == 1){
    sp <- "Cm"
  }
  if(s == 2){
    sp <- "Ls"
  }
  if(s == 3){
    sp <- "Ta"
  }

  gene2go <- read_eggnog_gene2go(eggnog_files[[sp]])
  gene2go_bp <- keep_go_bp(gene2go)
  go_bp_sets <- make_go_sets(gene2go_bp)

  for(j in 1:2){
    if(j == 1){
      testname <- paste0(sp, "_24vs00")
    }

    if(j == 2){
      testname <- paste0(sp, "_48vs24")
    }

    if (is.null(results[[testname]])) {
      stop(
        "Cannot find DESeq2 result: ", testname,
        "\nAvailable results: ",
        paste(names(results), collapse = ", ")
      )
    }

    res_df <- results[[testname]]

    ranked_gene_ids <- res_df %>%
      filter(!is.na(stat), is.finite(stat)) %>%
      transmute(gene_id = clean_id(id)) %>%
      distinct(gene_id) %>%
      pull(gene_id)
    bp_gene_ids <- unique(gene2go_bp$gene_id)
    annotated_ranked_genes <- sum(ranked_gene_ids %in% bp_gene_ids)
    annotation_coverage[[testname]] <- data.frame(
      comparison = testname,
      ranked_genes = length(ranked_gene_ids),
      ranked_genes_with_eggNOG_BP = annotated_ranked_genes,
      coverage_percent = round(100 * annotated_ranked_genes / length(ranked_gene_ids), 2),
      stringsAsFactors = FALSE
    )

    gsea_result <- run_fgsea_bp(
      res_df   = res_df,
      go_sets  = go_bp_sets,
      agg      = "max_abs"
    ) %>% add_go_term()

    gsea_barplots[[testname]] <- plot_fgsea(gsea_result, title = testname)
    fglist[[testname]] <- gsea_result
  }
}

annotation_coverage <- bind_rows(annotation_coverage)
write.table(
  annotation_coverage,
  "eggNOG_GO_BP_annotation_coverage.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

tag_fg <- function(fg, species, contrast) {
  fg %>%
    mutate(
      species  = species,
      contrast = contrast,
      term     = ifelse(!is.na(TERM), TERM, pathway)
    ) %>%
    dplyr::select(species, contrast, pathway, term, NES, padj)
}

fg_all <- bind_rows(
  tag_fg(fglist$Cm_24vs00, "Cm", "0-24"),
  tag_fg(fglist$Ls_24vs00, "Ls", "0-24"),
  tag_fg(fglist$Ta_24vs00, "Ta", "0-24"),

  tag_fg(fglist$Cm_48vs24, "Cm", "24-48"),
  tag_fg(fglist$Ls_48vs24, "Ls", "24-48"),
  tag_fg(fglist$Ta_48vs24, "Ta", "24-48")
)

# Preserve complete results before term selection for plotting.
fg_all_raw <- bind_rows(
  fglist$Cm_24vs00 %>% mutate(species = "Cm", contrast = "0-24"),
  fglist$Ls_24vs00 %>% mutate(species = "Ls", contrast = "0-24"),
  fglist$Ta_24vs00 %>% mutate(species = "Ta", contrast = "0-24"),
  fglist$Cm_48vs24 %>% mutate(species = "Cm", contrast = "24-48"),
  fglist$Ls_48vs24 %>% mutate(species = "Ls", contrast = "24-48"),
  fglist$Ta_48vs24 %>% mutate(species = "Ta", contrast = "24-48")
)

saveRDS(
  fglist,
  file.path(gsea_dir, "fgsea_raw_FINAL_fglist.rds")
)
saveRDS(
  fg_all_raw,
  file.path(gsea_dir, "fgsea_raw_FINAL_fg_all.rds")
)

pick_terms_from_Cm <- function(df, top_pos = 10, top_neg = 10, top_abs = 10, Cm_only = FALSE,
                               padj_cap = 1, require_sig = FALSE, padj_cut = 0.05) {
  if(Cm_only){
    d <- df %>%
      filter(species == "Cm") %>%
      mutate(padj2 = pmin(padj, padj_cap)) %>%
      filter(!is.na(NES))

  } else {
    d <- df %>%
      mutate(padj2 = pmin(padj, padj_cap)) %>%
      filter(!is.na(NES))
  }

  if (require_sig) d <- d %>% filter(!is.na(padj) & padj < padj_cut)

  picked <- d %>%
    group_by(species, contrast) %>%
    summarise(
      abs = list(pathway[order(-abs(NES))][seq_len(min(top_abs, sum(abs(NES))))]),
      pos = list(pathway[order(-NES)][seq_len(min(top_pos, sum(NES > 0)))]),
      neg = list(pathway[order( NES)][seq_len(min(top_neg, sum(NES < 0)))]),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(abs, pos, neg),
      names_to = "direction",
      values_to = "GO"
    ) %>%
    tidyr::unnest(GO)

  unique(picked)
}

collect_heatmap_data <- function(fg_all, Cm_only) {
  hm_list <- list()
  contrast_names <- c("0-24", "24-48")

  for (contrast_name in contrast_names) {
    fg_all2 <- fg_all %>%
      filter(.data$contrast == contrast_name)

    terms_keep <- pick_terms_from_Cm(
      fg_all2,
      Cm_only = Cm_only,
      top_abs = 0,
      top_pos = 10,
      top_neg = 10,
      require_sig = TRUE
    )

    hm <- fg_all2 %>%
      filter(pathway %in% terms_keep$GO) %>%
      dplyr::select(species, contrast, pathway, term, NES, padj) %>%
      distinct()

    term_order <- hm %>%
      arrange(species, NES) %>%
      distinct(term) %>%
      pull(term)

    hm_list[[contrast_name]] <- hm %>%
      mutate(
        term = factor(term, levels = rev(term_order)),
        sig = ifelse(!is.na(padj) & padj < 0.05, "*", "")
      )
  }

  bind_rows(hm_list) %>%
    dplyr::mutate(term_rank = as.integer(term))
}

# Main plot: Cm terms with padj < 0.05, positive and negative NES top 10 per contrast.
hm_collect_rankByCm <- collect_heatmap_data(fg_all, Cm_only = TRUE)

# Supplementary plot: apply the same selection independently to all three species.
hm_collect_rankByAllSP <- collect_heatmap_data(fg_all, Cm_only = FALSE)

write.table(hm_collect_rankByCm, "Data_for_GSEA_plot_rankByCm.txt", quote = FALSE,
            col.names = TRUE, row.names = FALSE, sep = "\t")

write.table(hm_collect_rankByAllSP, "Data_for_GSEA_plot_rankByAllSP.txt", quote = FALSE,
            col.names = TRUE, row.names = FALSE, sep = "\t")

}))

.run_release_stage(quote({
RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

project_root <- RELEASE_ROOT
gsea_dir <- file.path(project_root, "2.global_expression/Fig2B_GSEA")
results_dir <- file.path(gsea_dir, "results")

rank_by_cm_file <- file.path(
  results_dir,
  "Data_for_GSEA_plot_rankByCm.txt"
)
rank_by_allsp_file <- file.path(
  results_dir,
  "Data_for_GSEA_plot_rankByAllSP.txt"
)

plot_input_files <- c(rank_by_cm_file, rank_by_allsp_file)
if (!all(file.exists(plot_input_files))) {
  stop(
    "Missing plotting table: \n",
    paste(plot_input_files[!file.exists(plot_input_files)], collapse = "\n")
  )
}

required_columns <- c(
  "species", "contrast", "pathway", "term",
  "NES", "padj", "sig", "term_rank"
)

read_plot_table <- function(path) {
  x <- read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  missing_columns <- setdiff(required_columns, colnames(x))
  if (length(missing_columns) > 0) {
    stop(
      "Plotting table lacks columns: ",
      paste(missing_columns, collapse = ", "),
      "\nFile: ", path
    )
  }

  x %>%
    mutate(
      species = as.character(species),
      contrast = as.character(contrast),
      term = as.character(term),
      NES = as.numeric(NES),
      padj = as.numeric(padj),
      term_rank = as.numeric(term_rank),
      sig = ifelse(is.na(sig), "", as.character(sig))
    )
}

hm_collect_rank_by_cm <- read_plot_table(rank_by_cm_file)
hm_collect_rank_by_allsp <- read_plot_table(rank_by_allsp_file)

contrast_p1 <- c("0-24", "24-48")

make_term_levels <- function(x) {
  x %>%
    arrange(term_rank) %>%
    distinct(term) %>%
    pull(term)
}

# Keep GO labels separate across contrasts when ordering by Cm NES.
hm_collect_rank_by_cm <- hm_collect_rank_by_cm %>%
  mutate(
    term_label = str_wrap(term, width = 50),
    term_unique = paste(contrast, pathway, sep = "___")
  )

term_levels_cm <- hm_collect_rank_by_cm %>%
  filter(species == "Cm", contrast %in% contrast_p1) %>%
  arrange(
    factor(contrast, levels = contrast_p1),
    desc(NES),
    pathway
  ) %>%
  pull(term_unique)

hm_collect_rank_by_cm <- hm_collect_rank_by_cm %>%
  mutate(
    term_unique = factor(
      term_unique,
      levels = rev(term_levels_cm)
    )
  )

g_main <- ggplot(
  filter(hm_collect_rank_by_cm, contrast %in% contrast_p1),
  aes(x = species, y = term_unique, fill = NES)
) +
  geom_tile(color = "white", size = 0.2) +
  geom_text(
    aes(label = sig),
    size = 5,
    color = "black",
    fontface = "bold",
    vjust = 0.7,
    hjust = 0.5,
    na.rm = TRUE
  ) +
  scale_fill_gradient2(
    name = "NES",
    midpoint = 0,
    low = "#2166ac",
    mid = "#fefefe",
    high = "#b2182b",
    limits = c(-2, 2),
    oob = squish,
    breaks = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")
  ) +
  facet_grid(
    contrast ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) +
  scale_y_discrete(
    labels = function(x) {
      ids <- sub("^.*?___", "", x)
      str_wrap(
        hm_collect_rank_by_cm$term[
          match(ids, hm_collect_rank_by_cm$pathway)
        ],
        width = 50
      )
    }
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(
      size = 10, hjust = 1, color = "black", lineheight = 0.55
    ),
    axis.text.x = element_text(
      size = 10, color = "black", face = "italic"
    ),
    strip.text.x = element_text(
      size = 10, hjust = 0, color = "black"
    ),
    strip.background = element_blank(),
    legend.position = "none",
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.1, "cm"),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
ggsave(
  filename = file.path(results_dir, "GSEA_GO_BP_rankByCm_htmap.pdf"),
  plot = g_main,
  width = 5,
  height = 7.8,
  device = "pdf",
  dpi = 300
)

hm_collect_rank_by_allsp <- hm_collect_rank_by_allsp %>%
  mutate(term = str_wrap(term, width = 65))

term_levels_allsp <- make_term_levels(hm_collect_rank_by_allsp)
hm_collect_rank_by_allsp$term <- factor(
  hm_collect_rank_by_allsp$term,
  levels = term_levels_allsp
)

g_supp <- ggplot(
  filter(hm_collect_rank_by_allsp, contrast %in% contrast_p1),
  aes(x = species, y = term, fill = NES)
) +
  geom_tile(color = "white", size = 0.2) +
  geom_text(
    aes(label = sig),
    size = 4,
    color = "black",
    fontface = "bold",
    vjust = 0.7,
    hjust = 0.5,
    na.rm = TRUE
  ) +
  scale_fill_gradient2(
    name = "NES",
    midpoint = 0,
    low = "#2166ac",
    mid = "#fefefe",
    high = "#b2182b",
    limits = c(-2, 2),
    oob = squish,
    breaks = c(-2, -1, 0, 1, 2)
  ) +
  facet_grid(
    . ~ contrast,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_discrete(limits = rev(term_levels_allsp)) +
  scale_x_discrete(
    limits = unique(hm_collect_rank_by_allsp$species)
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(
      size = 10, hjust = 1, color = "black", lineheight = 0.6
    ),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.x = element_text(size = 10, color = "black",face = "italic"),
    strip.text.x = element_text(size = 12, color = "black"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

ggsave(
  filename = file.path(results_dir, "GSEA_GO_BP_rankByAllSP_htmap.pdf"),
  plot = g_supp,
  width = 8,
  height = 14,
  device = "pdf",
  dpi = 300
)

message("Main figure written to: ", file.path(results_dir, "GSEA_GO_BP_rankByCm_htmap.pdf"))
message("Supplementary figure written to: ", file.path(results_dir, "GSEA_GO_BP_rankByAllSP_htmap.pdf"))
}))
