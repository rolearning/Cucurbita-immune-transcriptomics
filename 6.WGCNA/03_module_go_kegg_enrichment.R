#!/usr/bin/env Rscript

# Test GO and KEGG over-representation in trait-associated modules of final model 0.5_22_2_50.
# Read the Stage 1 workspace, final-model membership/summary and Cm eggNOG annotations.
# Write enrichment CSVs, coverage, diagnostic plots and audit records to stage2_GO_kegg/.
# Use annotated network genes as the universe, BH correction and gene sets of 10-500 genes.

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(GO.db)
})

stage1_workspace <- file.path(RELEASE_ROOT, "6.WGCNA/stage1/stage1_workspace.rds")
final_stage2_dir <- file.path(RELEASE_ROOT, "6.WGCNA/stage2_only_0.5_22_2_50/0.5_22_2_50")
eggnog_annotation_file <- file.path(RELEASE_ROOT, "eggnog_Cm_Ls_Ta/Cmos_eggnog_anno.emapper.annotations")

output_dir <- file.path(RELEASE_ROOT, "6.WGCNA/stage2_GO_kegg/0.5_22_2_50/")
go_dir <- file.path(output_dir, "GO")
kegg_dir <- file.path(output_dir, "KEGG")

final_model <- "0.5_22_2_50"
padj_cutoff <- 0.05
min_gs_size <- 10L
max_gs_size <- 500L

required_inputs <- c(
  stage1_workspace,
  eggnog_annotation_file,
  file.path(final_stage2_dir, "formal_results.rds"),
  file.path(final_stage2_dir, "module_summary.csv")
)
if (!all(file.exists(required_inputs))) {
  stop(
    "Missing required input files: \n",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "\n")
  )
}

dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)

clean_gene_id <- function(x) {
  # Remove gene/transcript version suffixes before mapping annotations.
  sub("\\.[0-9]+$", "", as.character(x))
}

read_cm_eggnog_annotation <- function(annotation_file) {
  # Search the complete file because the #query header may follow many metadata lines.
  annotation_lines <- readLines(annotation_file, warn = FALSE)
  header_line <- grep("^#query\\t", annotation_lines)[1]
  if (length(header_line) != 1L || is.na(header_line)) {
    stop("No #query header found in the eggNOG annotation: ", annotation_file)
  }

  column_names <- strsplit(
    sub("^#", "", annotation_lines[header_line]),
    "\\t"
  )[[1]]

  annotation_df <- read.delim(
    annotation_file,
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

  required_columns <- c("query", "GOs", "KEGG_Pathway", "KEGG_ko")
  missing_columns <- setdiff(required_columns, colnames(annotation_df))
  if (length(missing_columns)) {
    stop(
      "The eggNOG annotation is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  annotation_df %>%
    transmute(
      Gene = clean_gene_id(query),
      GOs = as.character(GOs),
      KEGG_Pathway = as.character(KEGG_Pathway),
      KEGG_ko = as.character(KEGG_ko)
    ) %>%
    filter(!is.na(Gene), Gene != "") %>%
    distinct()
}

ann <- read_cm_eggnog_annotation(eggnog_annotation_file)

stage1 <- readRDS(stage1_workspace)
formal_results <- readRDS(file.path(final_stage2_dir, "formal_results.rds"))
if (!identical(as.character(formal_results$run_id), final_model)) {
  stop(
    "The run_id in formal_results.rds does not match final_model: ",
    formal_results$run_id
  )
}

module_summary <- read.csv(
  file.path(final_stage2_dir, "module_summary.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("module", "top_trait") %in% colnames(module_summary))) {
  stop("module_summary.csv is missing module or top_trait columns.")
}

run <- stage1$runs[[final_model]]
if (is.null(run)) stop("The final model was not found in the Stage 1 workspace: ", final_model)

pre <- stage1$preprocess[[run$keep_key]]
network_genes <- colnames(pre$datExpr)
colors <- run$net$colors

if (is.null(names(colors))) {
  stop("The final model net$colors vector lacks gene names.")
}

if (!setequal(network_genes, names(colors))) {
  stop("The final model network genes do not match the gene names in net$colors.")
}

network_genes <- as.character(network_genes)
names(colors) <- as.character(names(colors))

significant_modules <- module_summary$module[
  !is.na(module_summary$top_trait) &
    module_summary$top_trait != "" &
    module_summary$top_trait != "NA"
]
significant_modules <- intersect(significant_modules, unique(paste0("ME", colors)))

if (!length(significant_modules)) {
  stop("The final model has no trait-associated modules available for enrichment.")
}

# Keep GO ontologies, KEGG pathways and KO terms as separate annotation sources.

go_terms <- ann %>%
  dplyr::select(Gene, GOs) %>%
  filter(!is.na(GOs), GOs != "", GOs != "-") %>%
  separate_rows(GOs, sep = ",") %>%
  transmute(GO = trimws(GOs), Gene) %>%
  filter(grepl("^GO:[0-9]{7}$", GO)) %>%
  distinct()

go_info <- AnnotationDbi::select(
  GO.db::GO.db,
  keys = unique(go_terms$GO),
  keytype = "GOID",
  columns = c("TERM", "ONTOLOGY")
) %>%
  as_tibble() %>%
  filter(!is.na(GOID), !is.na(TERM)) %>%
  distinct(GOID, .keep_all = TRUE)

go_terms <- go_terms %>%
  inner_join(go_info, by = c("GO" = "GOID")) %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC"))

kegg_pathway_terms <- ann %>%
  dplyr::select(Gene, KEGG_Pathway) %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "", KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  transmute(term = trimws(KEGG_Pathway), gene = Gene) %>%
  filter(term != "", term != "-") %>%
  distinct()

kegg_pathway_names <- kegg_pathway_terms %>%
  distinct(term) %>%
  mutate(name = term)

kegg_ko_terms <- ann %>%
  dplyr::select(Gene, KEGG_ko) %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "", KEGG_ko != "-") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  transmute(term = trimws(KEGG_ko), gene = Gene) %>%
  filter(grepl("^ko:[Kk][0-9]{5}$|^[Kk][0-9]{5}$", term)) %>%
  mutate(term = ifelse(grepl("^ko:", term), term, paste0("ko:", term))) %>%
  distinct()

kegg_ko_names <- kegg_ko_terms %>%
  distinct(term) %>%
  mutate(name = term)

# Intersect foreground and universe with annotations for each tested source.

run_ora <- function(foreground, universe, term2gene, term2name,
                    module, source) {
  foreground <- intersect(unique(as.character(foreground)),
                          unique(as.character(universe)))
  universe <- unique(as.character(universe))

  term2gene <- term2gene %>%
    transmute(term = as.character(term), gene = as.character(gene)) %>%
    filter(!is.na(term), term != "", !is.na(gene), gene != "") %>%
    distinct()

  annotated_universe <- intersect(universe, unique(term2gene$gene))
  annotated_foreground <- intersect(foreground, annotated_universe)

  if (!length(annotated_foreground) || !length(annotated_universe) ||
      !nrow(term2gene)) {
    return(tibble())
  }

  result <- clusterProfiler::enricher(
    gene = annotated_foreground,
    universe = annotated_universe,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  )

  result_df <- as.data.frame(result)
  if (!nrow(result_df)) return(tibble())

  as_tibble(result_df) %>%
    mutate(
      final_model = final_model,
      module = module,
      source = source,
      n_foreground_annotated = length(annotated_foreground),
      n_universe_annotated = length(annotated_universe),
      .before = 1
    )
}

network_with_annotation <- list(
  GO = intersect(network_genes, unique(go_terms$Gene)),
  KEGG_Pathway = intersect(network_genes, unique(kegg_pathway_terms$gene)),
  KEGG_KO = intersect(network_genes, unique(kegg_ko_terms$gene))
)

all_go <- list()
all_kegg <- list()
coverage <- list()
go_index <- 1L
kegg_index <- 1L
coverage_index <- 1L

for (me in significant_modules) {
  module_name <- sub("^ME", "", me)
  module_genes <- network_genes[as.character(colors) == module_name]

  coverage[[coverage_index]] <- tibble(
    final_model = final_model,
    module = me,
    n_module_gene = length(module_genes),
    n_module_gene_with_GO = length(intersect(module_genes, network_with_annotation$GO)),
    n_module_gene_with_KEGG_Pathway = length(intersect(module_genes, network_with_annotation$KEGG_Pathway)),
    n_module_gene_with_KEGG_KO = length(intersect(module_genes, network_with_annotation$KEGG_KO)),
    n_universe_gene = length(network_genes),
    n_universe_gene_with_GO = length(network_with_annotation$GO),
    n_universe_gene_with_KEGG_Pathway = length(network_with_annotation$KEGG_Pathway),
    n_universe_gene_with_KEGG_KO = length(network_with_annotation$KEGG_KO)
  )
  coverage_index <- coverage_index + 1L

  for (ont in c("BP", "MF", "CC")) {
    ont_terms <- go_terms %>%
      filter(ONTOLOGY == ont) %>%
      transmute(term = GO, gene = Gene, name = TERM)
    go_result <- run_ora(
      module_genes,
      network_genes,
      ont_terms %>% dplyr::select(term, gene),
      ont_terms %>% dplyr::select(term, name) %>% distinct(),
      me,
      paste0("GO_", ont)
    )
    if (nrow(go_result)) {
      all_go[[go_index]] <- go_result
      go_index <- go_index + 1L
    }
  }

  pathway_result <- run_ora(
    module_genes, network_genes,
    kegg_pathway_terms, kegg_pathway_names,
    me, "KEGG_Pathway"
  )
  if (nrow(pathway_result)) {
    all_kegg[[kegg_index]] <- pathway_result
    kegg_index <- kegg_index + 1L
  }

  ko_result <- run_ora(
    module_genes, network_genes,
    kegg_ko_terms, kegg_ko_names,
    me, "KEGG_KO"
  )
  if (nrow(ko_result)) {
    all_kegg[[kegg_index]] <- ko_result
    kegg_index <- kegg_index + 1L
  }
}

go_all <- bind_rows(all_go)
kegg_all <- bind_rows(all_kegg)
coverage_df <- bind_rows(coverage)

go_significant <- go_all %>%
  filter(!is.na(p.adjust), p.adjust <= padj_cutoff)
kegg_significant <- kegg_all %>%
  filter(!is.na(p.adjust), p.adjust <= padj_cutoff)

write.csv(go_all, file.path(go_dir, "GO_enrichment_all_modules.csv"), row.names = FALSE)
write.csv(go_significant, file.path(go_dir, "GO_enrichment_significant.csv"), row.names = FALSE)
write.csv(kegg_all, file.path(kegg_dir, "KEGG_enrichment_all_modules.csv"), row.names = FALSE)
write.csv(kegg_significant, file.path(kegg_dir, "KEGG_enrichment_significant.csv"), row.names = FALSE)
write.csv(coverage_df, file.path(output_dir, "module_annotation_coverage.csv"), row.names = FALSE)

# These overview plots are diagnostic; publication panels are generated downstream.

gene_ratio_numeric <- function(x) {
  x <- as.character(x)
  numerator <- suppressWarnings(as.numeric(sub("/.*$", "", x)))
  denominator <- suppressWarnings(as.numeric(sub("^.*/", "", x)))
  numerator / denominator
}

plot_enrichment <- function(tbl, out_file, title_text, source_filter = NULL) {
  plot_tbl <- tbl
  if (!is.null(source_filter)) plot_tbl <- plot_tbl %>% filter(source %in% source_filter)
  plot_tbl <- plot_tbl %>% filter(!is.na(p.adjust))

  if (!nrow(plot_tbl)) {
    grDevices::pdf(out_file, width = 8, height = 6)
    plot.new()
    title(main = paste0(title_text, "\nNo enrichment terms returned"))
    grDevices::dev.off()
    return(invisible(NULL))
  }

  plot_tbl <- plot_tbl %>%
    mutate(
      GeneRatio_numeric = gene_ratio_numeric(GeneRatio),
      term_label = ifelse(
        source == "GO_BP" | source == "GO_MF" | source == "GO_CC",
        paste0(Description, " (", ID, ")"),
        Description
      )
    ) %>%
    group_by(module, source) %>%
    slice_min(order_by = p.adjust, n = 10L, with_ties = FALSE) %>%
    ungroup()

  p <- ggplot(
    plot_tbl,
    aes(x = GeneRatio_numeric,
        y = reorder(term_label, GeneRatio_numeric),
        size = Count,
        color = -log10(pmax(p.adjust, .Machine$double.xmin)))
  ) +
    geom_point() +
    facet_grid(source ~ module, scales = "free_y", space = "free_y") +
    scale_color_gradient(
      low = "#abdcff",
      high = "darkblue",
      name = "-log10 adjusted P"
    ) +
    labs(
      x = "Gene ratio",
      y = NULL,
      color = "-log10 adjusted P",
      size = "Gene count",
      title = title_text
    ) +
    theme_classic(base_size = 9, base_family = "Arial") +
    theme(
      axis.title.x = element_text(family = "Arial", face = "plain", size = 9),
      axis.title.y = element_text(family = "Arial", face = "plain", size = 9),
      axis.text.y = element_text(
        family = "Arial", face = "plain", size = 7.2, colour = "#222222"
      ),
      axis.text.x = element_text(
        family = "Arial", face = "plain", size = 8,
        angle = 45, hjust = 1, vjust = 1, colour = "#333333"
      ),
      axis.ticks = element_line(colour = "#BDB9B2", linewidth = 0.28),
      axis.line = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(
        colour = "#F0EEE9", fill = NA, linewidth = 0.32
      ),
      panel.grid.major = element_line(
        colour = "#F4F3F0", linewidth = 0.22
      ),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#F0EEE9", colour = "#E4E1DB"),
      strip.text = element_text(
        family = "Arial", face = "plain", size = 7.5, colour = "#333333"
      ),
      legend.title = element_text(
        family = "Arial", face = "plain", size = 8
      ),
      legend.text = element_text(
        family = "Arial", face = "plain", size = 7
      ),
      legend.key = element_blank(),
      legend.position = "right",
      plot.title = element_text(
        family = "Arial", face = "plain", size = 10, hjust = 0.5
      )
    )

  ggsave(
    out_file,
    p,
    width = max(9, 2.5 * length(unique(plot_tbl$module))),
    height = max(6, 0.20 * nrow(plot_tbl) + 3),
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
}

plot_enrichment(
  go_all,
  file.path(go_dir, "GO_enrichment_plot.pdf"),
  paste0("eggNOG GO enrichment: ", final_model)
)
plot_enrichment(
  kegg_all,
  file.path(kegg_dir, "KEGG_enrichment_plot.pdf"),
  paste0("eggNOG KEGG enrichment: ", final_model),
  c("KEGG_Pathway", "KEGG_KO")
)

audit_lines <- c(
  paste0("final_model: ", final_model),
  paste0("stage1_workspace: ", stage1_workspace),
  paste0("final_stage2_dir: ", final_stage2_dir),
  paste0("eggnog_annotation_file: ", eggnog_annotation_file),
  paste0("n_WGCNA_network_genes: ", length(network_genes)),
  paste0("n_GO_annotated_network_genes: ", length(network_with_annotation$GO)),
  paste0("n_KEGG_Pathway_annotated_network_genes: ", length(network_with_annotation$KEGG_Pathway)),
  paste0("n_KEGG_KO_annotated_network_genes: ", length(network_with_annotation$KEGG_KO)),
  paste0("GO_background_definition: final WGCNA network genes intersected with eggNOG GO-annotated genes"),
  paste0("KEGG_Pathway_background_definition: final WGCNA network genes intersected with eggNOG KEGG_Pathway-annotated genes"),
  paste0("KEGG_KO_background_definition: final WGCNA network genes intersected with eggNOG KEGG_ko-annotated genes"),
  paste0("n_GO_terms_all: ", nrow(go_all)),
  paste0("n_GO_terms_significant: ", nrow(go_significant)),
  paste0("n_KEGG_terms_all: ", nrow(kegg_all)),
  paste0("n_KEGG_terms_significant: ", nrow(kegg_significant)),
  paste0("ORA_method: clusterProfiler::enricher with custom TERM2GENE/TERM2NAME"),
  paste0("multiple_testing_correction: BH"),
  paste0("GO/KEGG_significance_threshold: adjusted P <= ", padj_cutoff),
  paste0("minGSSize: ", min_gs_size),
  paste0("maxGSSize: ", max_gs_size),
  "gene_ID_rule: eggNOG query IDs were cleaned by removing terminal numeric version suffix",
  "GO_source: eggNOG-mapper GOs column, ontology and term names resolved with GO.db",
  "KEGG_pathway_source: eggNOG-mapper KEGG_Pathway column",
  "KEGG_KO_source: eggNOG-mapper KEGG_ko column",
  "foreground_definition: genes in each final-model trait-associated module intersected with the relevant annotated universe",
  "module_assignment: read from existing final-model Stage 1 network; not recomputed",
  "analysis_scope: final model 0.5_22_2_50 only",
  "",
  "sessionInfo:",
  capture.output(sessionInfo())
)
writeLines(audit_lines, file.path(output_dir, "enrichment_audit.txt"))

message("Final-model eggNOG GO/KEGG enrichment script completed.")
message("Output directory: ", output_dir)
