#!/usr/bin/env Rscript

# Enrich Fig. 2a clusters using species-specific eggNOG GO annotations.
# Read cluster membership, centroids, and the syngroup ID map from clustering.
# Test BP, CC, and MF with BH adjustment for the TMM and median branches.
# Write enrichment tables, category summaries, BP plots, and separate legends.

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

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(GO.db)
})

project_dir <- RELEASE_ROOT
go_dir <- file.path(project_dir, "2.global_expression/Fig2A_GO/GO")
kmeans_dir <- file.path(project_dir, "2.global_expression/Fig2A_GO/kmeans")
eggnog_dir <- file.path(project_dir, "eggnog_Cm_Ls_Ta")
id_file <- file.path(kmeans_dir, "All_syngroup_ID_list.csv")

METHODS <- c("TMM", "median")
HOURS <- c("00", "24", "48")
SPECIES <- c("Cm", "Ls", "Ta")
ONTOLOGIES <- c("BP", "CC", "MF")
PADJ_CUTOFF <- 0.05
MIN_GS_SIZE <- 10L
MAX_GS_SIZE <- 500L

dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(go_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(go_dir, "plot"), recursive = TRUE, showWarnings = FALSE)

eggnog_files <- c(
  Cm = file.path(eggnog_dir, "Cmos_eggnog_anno.emapper.annotations"),
  Ls = file.path(eggnog_dir, "Lsic_eggnog_anno.emapper.annotations"),
  Ta = file.path(eggnog_dir, "Tang_eggnog_anno.emapper.annotations")
)

required_files <- c(id_file, eggnog_files)
if (!all(file.exists(required_files))) {
  stop("Missing GO input files: ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}

clean_gene_id <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

read_eggnog_gene2go <- function(annotation_file) {
  annotation_lines <- readLines(annotation_file, warn = FALSE)
  header_line <- grep("^#query\\t", annotation_lines)[1]
  if (length(header_line) != 1L || is.na(header_line)) {
    stop("The eggNOG annotation lacks a #query header: ", annotation_file)
  }
  column_names <- strsplit(sub("^#", "", annotation_lines[header_line]), "\\t")[[1]]
  annotation_df <- read.delim(
    annotation_file, sep = "\t", header = FALSE, skip = header_line,
    comment.char = "#", quote = "", fill = TRUE,
    col.names = column_names, check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!all(c("query", "GOs") %in% colnames(annotation_df))) {
    stop("The eggNOG annotation lacks query or GOs columns: ", annotation_file)
  }

  gene2go <- annotation_df %>%
    transmute(gene_id = clean_gene_id(query), GO = GOs) %>%
    filter(!is.na(gene_id), gene_id != "", !is.na(GO), GO != "", GO != "-") %>%
    separate_rows(GO, sep = ",") %>%
    filter(grepl("^GO:[0-9]{7}$", GO)) %>%
    distinct(gene_id, GO)

  go_info <- AnnotationDbi::select(
    GO.db::GO.db, keys = unique(gene2go$GO), keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
  ) %>%
    as_tibble() %>%
    filter(!is.na(GOID), !is.na(TERM), ONTOLOGY %in% ONTOLOGIES) %>%
    distinct(GOID, .keep_all = TRUE)

  out <- gene2go %>%
    inner_join(go_info, by = c("GO" = "GOID")) %>%
    distinct(gene_id, GO, TERM, ONTOLOGY)
  if (nrow(out) == 0L) stop("No usable eggNOG GO annotations found: ", annotation_file)
  out
}

syngroup_to_species_genes <- function(syngroup_ids, id_map, species) {
  species_columns <- list(Cm = c("Cm1", "Cm2"), Ls = "Lsi", Ta = "Tan")
  if (!species %in% names(species_columns)) stop("Unknown species: ", species)
  missing_ids <- setdiff(syngroup_ids, id_map$unique_id)
  if (length(missing_ids) > 0L) {
    stop("Cluster syngroups are absent from the ID table: ", paste(head(missing_ids, 10), collapse = ", "))
  }
  selected <- id_map[id_map$unique_id %in% syngroup_ids, species_columns[[species]], drop = FALSE]
  genes <- clean_gene_id(unlist(selected, use.names = FALSE))
  sort(unique(genes[!is.na(genes) & genes != ""]))
}

select_species_high_clusters <- function(pattern_df, species) {
  centroid_columns <- c("Cm_centroid", "Ls_centroid", "Ta_centroid")
  if (!all(c("cluster", centroid_columns) %in% colnames(pattern_df))) {
    stop("cluster_pattern lacks cluster or centroid columns.")
  }
  centroid_matrix <- as.matrix(pattern_df[, centroid_columns, drop = FALSE])
  storage.mode(centroid_matrix) <- "numeric"
  # Assign each cluster to its highest-centroid species; ties select the first species.
  highest_species <- c("Cm", "Ls", "Ta")[max.col(centroid_matrix, ties.method = "first")]
  as.integer(pattern_df$cluster[highest_species == species])
}

run_one_go_ora <- function(foreground_genes, universe_genes, gene2go, ontology) {
  annotation_ont <- gene2go %>%
    filter(ONTOLOGY == ontology) %>%
    distinct(gene_id, GO, TERM)
  # Restrict the clustered-gene universe and foreground to species-specific GO annotation.
  annotated_universe <- intersect(unique(universe_genes), unique(annotation_ont$gene_id))
  annotated_foreground <- intersect(unique(foreground_genes), annotated_universe)
  if (length(annotated_foreground) == 0L || length(annotated_universe) == 0L) {
    return(data.frame())
  }

  # Test gene sets of size 10-500; apply BH p.adjust <= 0.05 downstream.
  result <- clusterProfiler::enricher(
    gene = annotated_foreground,
    universe = annotated_universe,
    TERM2GENE = annotation_ont %>% dplyr::select(GO, gene_id),
    TERM2NAME = annotation_ont %>% dplyr::select(GO, TERM) %>% dplyr::distinct(),
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE
  )
  result <- as.data.frame(result)
  if (nrow(result) == 0L) return(data.frame())
  result$ONTOLOGY <- ontology
  result
}

id_map <- read.csv(id_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("unique_id", "Cm1", "Cm2", "Lsi", "Tan") %in% colnames(id_map))) {
  stop("All_syngroup_ID_list.csv lacks required columns.")
}
gene2go_by_species <- lapply(eggnog_files, read_eggnog_gene2go)

all_enrichment <- list()
coverage_rows <- list()
selection_rows <- list()
centroid_audit_rows <- list()
result_index <- 1L
coverage_index <- 1L
selection_index <- 1L
centroid_index <- 1L

for (method in METHODS) {
  method_kmeans_dir <- file.path(kmeans_dir, method)
  method_result_dir <- file.path(go_dir, "results", method)
  dir.create(method_result_dir, recursive = TRUE, showWarnings = FALSE)

  for (hour in HOURS) {
    cluster_file <- file.path(method_kmeans_dir,
                              paste0("all_syngroup_expr_", method, "_", hour, "_cluster.txt"))
    pattern_file <- file.path(method_kmeans_dir,
                              paste0("cluster_pattern_", method, "_", hour, ".tsv"))
    if (!all(file.exists(c(cluster_file, pattern_file)))) {
      stop("Missing cluster or cluster_pattern input: ", method, " ", hour)
    }

    cluster_df <- read.table(cluster_file, header = TRUE, row.names = 1, sep = "\t",
                             quote = "", comment.char = "", check.names = FALSE)
    pattern_df <- read.table(pattern_file, header = TRUE, sep = "\t",
                             quote = "", comment.char = "", check.names = FALSE)
    if (!"cluster" %in% colnames(cluster_df)) stop("Cluster file lacks the cluster column: ", cluster_file)
    universe_syngroups <- rownames(cluster_df)

    for (species in SPECIES) {
      high_clusters <- select_species_high_clusters(pattern_df, species)
      foreground_syngroups <- universe_syngroups[cluster_df$cluster %in% high_clusters]
      universe_genes <- syngroup_to_species_genes(universe_syngroups, id_map, species)
      foreground_genes <- syngroup_to_species_genes(foreground_syngroups, id_map, species)
      gene2go <- gene2go_by_species[[species]]
      annotated_gene_ids <- unique(gene2go$gene_id)

      annotated_universe_genes <- intersect(universe_genes, annotated_gene_ids)
      annotated_foreground_genes <- intersect(foreground_genes, annotated_universe_genes)

      centroid_audit_rows[[centroid_index]] <- pattern_df %>%
        transmute(
          method = method,
          hour = hour,
          species = species,
          cluster = as.integer(cluster),
          expression_order,
          Cm_centroid,
          Ls_centroid,
          Ta_centroid,
          selected_for_species = as.integer(cluster) %in% high_clusters
        )
      centroid_index <- centroid_index + 1L

      gene_set_dir <- file.path(method_result_dir, paste0(hour, "hpi"), species)
      dir.create(gene_set_dir, recursive = TRUE, showWarnings = FALSE)
      write.table(data.frame(gene_id = annotated_foreground_genes),
                  file.path(gene_set_dir, "foreground_genes.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(data.frame(gene_id = annotated_universe_genes),
                  file.path(gene_set_dir, "universe_genes.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      selection_rows[[selection_index]] <- data.frame(
        method = method, role = "pending_final_decision", hour = hour, species = species,
        high_clusters = paste(sort(high_clusters), collapse = ","),
        n_foreground_syngroup = length(foreground_syngroups),
        n_universe_syngroup = length(universe_syngroups), stringsAsFactors = FALSE
      )
      selection_index <- selection_index + 1L

      coverage_rows[[coverage_index]] <- data.frame(
        method = method, hour = hour, species = species,
        n_foreground_gene = length(foreground_genes),
        n_foreground_with_eggNOG_GO = length(annotated_foreground_genes),
        foreground_GO_coverage = ifelse(length(foreground_genes) > 0L,
                                        length(annotated_foreground_genes) / length(foreground_genes), NA_real_),
        n_universe_gene = length(universe_genes),
        n_universe_with_eggNOG_GO = length(annotated_universe_genes),
        universe_GO_coverage = ifelse(length(universe_genes) > 0L,
                                      length(annotated_universe_genes) / length(universe_genes), NA_real_),
        stringsAsFactors = FALSE
      )
      coverage_index <- coverage_index + 1L

      for (ontology in ONTOLOGIES) {
        enrichment_df <- run_one_go_ora(
          annotated_foreground_genes, annotated_universe_genes, gene2go, ontology
        )
        if (nrow(enrichment_df) == 0L) next
        enrichment_df <- enrichment_df %>% mutate(
          method = method, role = "pending_final_decision", hour = hour, species = species,
          high_clusters = paste(sort(high_clusters), collapse = ","), .before = 1
        )
        all_enrichment[[result_index]] <- enrichment_df
        result_index <- result_index + 1L
      }
    }
  }
}

combined_enrichment <- bind_rows(all_enrichment)
combined_significant <- combined_enrichment %>%
  filter(!is.na(p.adjust), p.adjust <= PADJ_CUTOFF)
write.table(combined_enrichment, file.path(go_dir, "GO_enrichment_all_terms.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(combined_significant, file.path(go_dir, "GO_enrichment_padj0.05.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(selection_rows), file.path(go_dir, "high_cluster_selection.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(coverage_rows), file.path(go_dir, "eggNOG_GO_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(centroid_audit_rows),
            file.path(go_dir, "cluster_centroid_selection_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Functional categories affect display only, leaving enrichment statistics unchanged.
classify_go_term <- function(description) {
  term <- str_to_lower(description)
  case_when(
    str_detect(term, "pathogen|microbe|elicitor|pattern recognition|prr|pamp|mamp|damp|pti|eti|effector|avr|chitin|flagellin|receptor-like kinase|receptor-like protein|lysm|lectin|lrr|leucine-rich repeat|nb-arc|nucleotide-binding|nbs|tir") ~ "Pathogen recognition",
    str_detect(term, "immune|defense|resistance|disease resistance|pathogen response|hypersensitive response|systemic acquired resistance|antimicrobial|pathogenesis-related|phytoalexin|callose|lignin|phenylpropanoid|flavonoid|glucosinolate|camalexin|cell wall|biotic stress|reactive oxygen species|oxidative stress|programmed cell death|cell death|ethylene|jasmon|salicy|abscisic") ~ "Immune defense response",
    str_detect(term, "signal|signaling pathway|signal transduction|transduction|activation|phosphorylation|kinase|phosphatase|phosphorelay|cascade|mapk|calcium|calmodulin|cdpk|cbl|cipk|g protein|gtpase|small gtpase|phospholipase|phosphatidylinositol|ubiquitin|sumoylation") ~ "Signal transduction or activation",
    str_detect(term, "morphogenesis|development|growth|organ formation|organogenesis|embryogenesis|germination|senescence|differentiation|pattern formation|histogenesis|tropism|meristem|root|shoot|leaf|flower|pollen|seed|embryo|vascular|xylem|phloem|cell cycle|mitotic|cell division|cell proliferation|cytokinesis|auxin|cytokinin|brassinosteroid|strigolactone") ~ "Growth and development",
    str_detect(term, "biosynthetic|catabolic process|amino acid|nucleotide|purine|pyrimidine|vitamin|cofactor|lipid|carbohydrate|photosynthetic membrane|thylakoid|plastid|chloroplast|ribosome|nucleolus|endosome|peroxisome|membrane|organelle|protein complex|chromatin|mitochondrion|binding|protein|acetylation|methylation") ~ "Primary metabolism",
    TRUE ~ "Other"
  )
}

category_levels <- c("Growth and development", "Primary metabolism",
                     "Pathogen recognition", "Signal transduction or activation",
                     "Immune defense response")
category_colors <- c(
  "Growth and development" = "#659157",
  "Primary metabolism" = "#9FBD82",
  "Pathogen recognition" = "#69A2B0",
  "Signal transduction or activation" = "#FFC1CF",
  "Immune defense response" = "#DD5162"
)

gene_ratio_numeric <- function(x) {
  x <- as.character(x)
  direct <- suppressWarnings(as.numeric(x))
  numerator <- suppressWarnings(as.numeric(sub("/.*$", "", x)))
  denominator <- suppressWarnings(as.numeric(sub("^.*/", "", x)))
  ifelse(is.finite(direct), direct, numerator / denominator)
}

go_bp <- combined_enrichment %>%
  filter(ONTOLOGY == "BP", !is.na(p.adjust), p.adjust <= PADJ_CUTOFF) %>%
  mutate(
    Description_group = classify_go_term(Description),
    GeneRatio_numeric = gene_ratio_numeric(GeneRatio),
    term_label = paste0(Description, " (", ID, ")")
  ) %>%
  group_by(method, hour, species, ID) %>%
  slice_min(order_by = p.adjust, n = 1L, with_ties = FALSE) %>%
  ungroup()

write.table(go_bp, file.path(go_dir, "GO_BP_classified_all_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Count unique foreground genes across significant terms within each functional category.
category_summary <- go_bp %>%
  filter(Description_group %in% category_levels) %>%
  dplyr::select(method, hour, species, Description_group, ID, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  filter(!is.na(geneID), geneID != "") %>%
  group_by(method, hour, species, Description_group) %>%
  summarise(
    n_terms = n_distinct(ID),
    n_unique_gene = n_distinct(geneID),
    .groups = "drop"
  )

summary_grid <- expand_grid(
  method = METHODS, hour = HOURS, species = SPECIES,
  Description_group = category_levels
) %>%
  left_join(category_summary,
            by = c("method", "hour", "species", "Description_group")) %>%
  mutate(n_terms = replace_na(n_terms, 0L),
         n_unique_gene = replace_na(n_unique_gene, 0L)) %>%
  left_join(bind_rows(coverage_rows) %>%
              dplyr::select(method, hour, species, n_foreground_with_eggNOG_GO),
            by = c("method", "hour", "species")) %>%
  left_join(bind_rows(selection_rows) %>%
              dplyr::select(method, hour, species, high_clusters),
            by = c("method", "hour", "species")) %>%
  mutate(GeneRatio = ifelse(n_foreground_with_eggNOG_GO > 0,
                            n_unique_gene / n_foreground_with_eggNOG_GO, 0),
         Genome = paste0(species, "C", gsub(",", "C", high_clusters)))

write.table(summary_grid, file.path(go_dir, "GO_category_summary_source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(go_dir, "GO_sessionInfo.txt"))
message("GO enrichment and category-summary source data are ready for the plotting stage.")
}))

.run_release_stage(quote({
RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggh4x)
})

project_dir <- RELEASE_ROOT
go_dir <- file.path(project_dir, "2.global_expression/Fig2A_GO/GO")
plot_dir <- file.path(go_dir, "plot")
source_file <- file.path(go_dir, "GO_category_summary_source_data.tsv")
bp_file <- file.path(go_dir, "GO_BP_classified_all_significant.tsv")

required_files <- c(source_file, bp_file)
if (!all(file.exists(required_files))) {
  stop("Cannot find GO plotting inputs: ",
       paste(required_files[!file.exists(required_files)], collapse = ", "))
}
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

summary_grid <- read.delim(
  source_file, header = TRUE, sep = "\t", quote = "",
  comment.char = "", check.names = FALSE, stringsAsFactors = FALSE
)

required_columns <- c(
  "method", "hour", "species", "Description_group", "n_terms",
  "n_unique_gene", "n_foreground_with_eggNOG_GO", "high_clusters",
  "GeneRatio", "Genome"
)
if (!all(required_columns %in% colnames(summary_grid))) {
  stop("GO category-summary source data lack required columns: ",
       paste(setdiff(required_columns, colnames(summary_grid)), collapse = ", "))
}

METHODS <- c("TMM", "median")
HOURS <- c("00", "24", "48")
SPECIES <- c("Cm", "Ls", "Ta")
CATEGORY_LEVELS <- c(
  "Growth and development",
  "Primary metabolism",
  "Pathogen recognition",
  "Signal transduction or activation",
  "Immune defense response"
)
CATEGORY_COLORS <- c(
  "Growth and development" = "#659157",
  "Primary metabolism" = "#9FBD82",
  "Pathogen recognition" = "#69A2B0",
  "Signal transduction or activation" = "#FFC1CF",
  "Immune defense response" = "#DD5162"
)

summary_grid$hour <- sprintf("%02d", as.integer(summary_grid$hour))
summary_grid$GeneRatio <- as.numeric(summary_grid$GeneRatio)
summary_grid <- summary_grid %>%
  filter(method %in% METHODS, hour %in% HOURS, species %in% SPECIES) %>%
  mutate(
    Description_group = factor(Description_group, levels = CATEGORY_LEVELS),
    x_key = paste(hour, species, sep = "_"),
    hour_label = factor(hour, levels = HOURS, labels = c("0", "24", "48"))
  )

expected_n <- length(METHODS) * length(HOURS) * length(SPECIES) * length(CATEGORY_LEVELS)
if (nrow(summary_grid) != expected_n) {
  stop("Unexpected GO category-summary row count; expected one row per method, time point, species, and category.")
}
if (any(!is.finite(summary_grid$GeneRatio)) || any(summary_grid$GeneRatio < 0)) {
  stop("GeneRatio contains non-finite or negative values.")
}

make_category_plot <- function(method_name) {
  plot_data <- summary_grid %>%
    filter(method == method_name) %>%
    mutate(
      Description_group = factor(Description_group, levels = rev(CATEGORY_LEVELS)),
      x_key = factor(
        x_key,
        levels = as.vector(outer(HOURS, SPECIES, paste, sep = "_"))
      )
    )

  x_label_data <- plot_data %>%
    distinct(x_key, Genome)
  x_labels <- stats::setNames(x_label_data$Genome, as.character(x_label_data$x_key))

  max_ratio <- max(plot_data$GeneRatio, na.rm = TRUE)
  size_upper <- max(0.20, max_ratio)

  ggplot(plot_data, aes(x = x_key, y = Description_group)) +
    geom_blank() +
    # Empty combinations retain their axes but receive no artificial zero-ratio points.
    geom_point(
      data = plot_data %>% filter(GeneRatio > 0),
      aes(size = GeneRatio, fill = Description_group),
      shape = 21, color = "transparent", stroke = 0
    ) +
    facet_grid(
      . ~ hour_label,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_x_discrete(labels = x_labels, drop = TRUE) +
    scale_fill_manual(
      values = CATEGORY_COLORS,
      breaks = CATEGORY_LEVELS,
      drop = FALSE,
      name = "Description group"
    ) +
    scale_size_continuous(
      range = c(2.5, 10),
      limits = c(0, size_upper),
      breaks = c(0.05, 0.10, 0.15, 0.20),
      name = "GeneRatio"
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 14, base_family = "Arial") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "#E5E5E5"),
      panel.border = element_rect(linewidth = 0.45, color = "#333333"),
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1, size = 14, color = "#333333"
      ),
      axis.text.y = element_text(size = 14, color = "#333333"),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 16, color = "#222222"),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(size = 5)),
      size = guide_legend(
        order = 2,
        override.aes = list(shape = 16, colour = "black", fill = "black")
      )
    )
}

save_category_plot <- function(plot_object, output_stem) {
  grDevices::cairo_pdf(
    paste0(output_stem, ".pdf"),
    width = 11.5,
    height = 4.8,
    family = "Arial"
  )
  print(plot_object)
  grDevices::dev.off()
}

save_category_plot(
  make_category_plot("TMM"),
  file.path(plot_dir, "GO_category_summary_TMM")
)
save_category_plot(
  make_category_plot("median"),
  file.path(plot_dir, "GO_category_summary_median")
)

bp_result <- read.delim(
  bp_file, header = TRUE, sep = "\t", quote = "",
  comment.char = "", check.names = FALSE, stringsAsFactors = FALSE
)
required_bp_columns <- c(
  "method", "hour", "species", "ID", "Description", "GeneRatio",
  "p.adjust", "high_clusters", "Description_group"
)
if (!all(required_bp_columns %in% colnames(bp_result))) {
  stop("BP enrichment plotting input lacks required columns: ",
       paste(setdiff(required_bp_columns, colnames(bp_result)), collapse = ", "))
}

bp_result$hour <- sprintf("%02d", as.integer(bp_result$hour))
bp_result$GeneRatio_numeric <- {
  ratio_text <- as.character(bp_result$GeneRatio)
  direct <- suppressWarnings(as.numeric(ratio_text))
  numerator <- suppressWarnings(as.numeric(sub("/.*$", "", ratio_text)))
  denominator <- suppressWarnings(as.numeric(sub("^.*/", "", ratio_text)))
  ifelse(is.finite(direct), direct, numerator / denominator)
}
bp_result$term_label <- paste0(bp_result$Description, " (", bp_result$ID, ")")
bp_result <- bp_result %>%
  filter(method %in% METHODS, hour %in% HOURS, species %in% SPECIES,
         is.finite(GeneRatio_numeric), GeneRatio_numeric > 0,
         Description_group %in% CATEGORY_LEVELS) %>%
  group_by(method, hour, species, ID) %>%
  slice_min(order_by = p.adjust, n = 1L, with_ties = FALSE) %>%
  ungroup()

# Share numerical scales and fixed GO-term row heights across methods and time points.
bp_row_height_in <- 0.34
bp_panel_width_in <- 1.8
bp_ratio_limits <- c(0, max(c(0, bp_result$GeneRatio_numeric), na.rm = TRUE))
if (bp_ratio_limits[2] <= 0) bp_ratio_limits[2] <- 1

bp_ratio_breaks <- c(0.03, 0.06, 0.09)
bp_ratio_breaks <- bp_ratio_breaks[
  bp_ratio_breaks > 0 & bp_ratio_breaks <= bp_ratio_limits[2]
]

make_bp_dotplot <- function(method_name, hour_name) {
  plot_data <- bp_result %>%
    filter(method == method_name, hour == hour_name) %>%
    mutate(
      Description_group = factor(Description_group, levels = CATEGORY_LEVELS),
      Genome = paste0(species, "C", gsub(",", "C", high_clusters))
    )
  if (nrow(plot_data) == 0L) return(NULL)

  genome_levels <- plot_data %>%
    distinct(species, Genome) %>%
    mutate(species = factor(species, levels = SPECIES)) %>%
    arrange(species) %>%
    pull(Genome)
  term_order <- plot_data %>%
    group_by(Description_group, ID, term_label) %>%
    summarise(n_species = n_distinct(species), .groups = "drop") %>%
    arrange(Description_group, desc(n_species), desc(term_label)) %>%
    pull(term_label) %>% unique()
  plot_data$Genome <- factor(plot_data$Genome, levels = genome_levels)
  plot_data$term_label <- factor(plot_data$term_label, levels = term_order)
  present_categories <- levels(droplevels(plot_data$Description_group))
  strip_theme <- ggh4x::strip_themed(
    background_y = ggh4x::elem_list_rect(fill = unname(CATEGORY_COLORS[present_categories]))
  )

  ggplot(plot_data, aes(x = Genome, y = term_label)) +
    geom_point(aes(size = GeneRatio_numeric), fill = "#444444",
               shape = 21, stroke = 0) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(position = "right", expand = expansion(add = 0.5)) +
    scale_size_continuous(
      range = c(4, 12), limits = bp_ratio_limits, breaks = bp_ratio_breaks
    ) +
    labs(x = NULL, y = NULL,
         title = paste0(ifelse(hour_name == "00", "0", hour_name), " hpi"),
         size = "Gene ratio") +
    theme_bw(base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5,
                                 size = 15, face = "italic", color = "black"),
      axis.text.y.left = element_blank(),
      axis.text.y.right = element_text(size = 15, color = "black"),
      panel.grid = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(color = "transparent"),
      plot.title = element_text(hjust = 0.5, size = 16),
      legend.position = "none"
    ) +
    ggh4x::facet_grid2(
      rows = vars(Description_group), scales = "free_y", space = "free",
      switch = "y", strip = strip_theme
    )
}

bp_grobs <- list()
for (method_name in METHODS) {
  for (hour_name in HOURS) {
    bp_plot <- make_bp_dotplot(method_name, hour_name)
    if (is.null(bp_plot)) next
    bp_build <- ggplot_build(bp_plot)
    bp_grob <- ggplotGrob(bp_plot)
    panel_cells <- bp_grob$layout[grepl("^panel", bp_grob$layout$name), ]
    panel_rows <- sort(unique(panel_cells$t))
    panel_cols <- sort(unique(panel_cells$l))
    facet_layout <- bp_build$layout$layout
    facet_rows <- sort(unique(facet_layout$ROW))
    if (length(panel_rows) != length(facet_rows)) {
      stop("BP facet rows and plot panels differ; cannot set a fixed row height.")
    }
    for (i in seq_along(panel_rows)) {
      scale_id <- unique(facet_layout$SCALE_Y[facet_layout$ROW == facet_rows[i]])
      if (length(scale_id) != 1L) stop("A BP facet contains multiple Y scales.")
      n_terms <- length(bp_build$layout$panel_scales_y[[scale_id]]$get_limits())
      bp_grob$heights[panel_rows[i]] <- grid::unit(n_terms * bp_row_height_in, "in")
    }
    bp_grob$widths[panel_cols] <- grid::unit(bp_panel_width_in, "in")
    bp_grobs[[paste0("GO_BP_dotplot_", method_name, "_", hour_name)]] <- bp_grob
  }
}

if (length(bp_grobs)) {
  width_lengths <- vapply(bp_grobs, function(g) length(g$widths), integer(1))
  if (length(unique(width_lengths)) != 1L) {
    stop("BP plots have different horizontal layouts; widths cannot be aligned safely.")
  }
  # Equal label widths and bottom padding align panels without stretching their rows.
  common_widths <- do.call(grid::unit.pmax, lapply(bp_grobs, function(g) g$widths))
  content_heights <- vapply(bp_grobs, function(g) {
    grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  }, numeric(1))
  page_height <- max(4, content_heights)
  page_width <- max(12, grid::convertWidth(sum(common_widths), "in", valueOnly = TRUE))
  for (plot_name in names(bp_grobs)) {
    bp_grob <- bp_grobs[[plot_name]]
    bp_grob$widths <- common_widths
    padding <- max(0, page_height - content_heights[[plot_name]])
    bp_grob <- gtable::gtable_add_rows(bp_grob, grid::unit(padding, "in"), pos = -1)
    ggsave(
      file.path(plot_dir, paste0(plot_name, ".pdf")),
      bp_grob, width = page_width, height = page_height, units = "in",
      limitsize = FALSE, device = grDevices::cairo_pdf, bg = "white"
    )
  }
}

bp_legend_data <- data.frame(
  Description_group = factor(CATEGORY_LEVELS, levels = CATEGORY_LEVELS),
  GeneRatio_numeric = seq(
    max(bp_ratio_limits[2] / length(CATEGORY_LEVELS), .Machine$double.eps),
    bp_ratio_limits[2],
    length.out = length(CATEGORY_LEVELS)
  )
)

bp_legend_plot <- ggplot(
  bp_legend_data,
  aes(x = Description_group, y = 1,
      size = GeneRatio_numeric, fill = Description_group)
) +
  geom_point(shape = 21, color = "transparent", stroke = 0) +
  scale_fill_manual(
    values = CATEGORY_COLORS,
    breaks = CATEGORY_LEVELS,
    drop = FALSE,
    name = "Description group"
  ) +
  scale_size_continuous(
    range = c(4, 12),
    limits = bp_ratio_limits,
    breaks = bp_ratio_breaks,
    name = "Gene ratio"
  ) +
  guides(
    size = guide_legend(
      order = 1,
      direction = "vertical",
      title.position = "left",
      override.aes = list(shape = 16, colour = "#444444", fill = "#444444")
    ),
    fill = guide_legend(
      order = 2,
      direction = "vertical",
      title.position = "left",
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 5, colour = NA)
    )
  ) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.key.width = grid::unit(0.55, "in"),
    legend.spacing.y = grid::unit(0.10, "in"),
    legend.margin = margin(0, 0, 0, 0)
  )

bp_legend_grob <- ggplotGrob(bp_legend_plot)
bp_guide_index <- which(bp_legend_grob$layout$name == "guide-box-bottom")
if (length(bp_guide_index) != 1L) {
  stop("Cannot extract a separate legend from the BP legend object.")
}
bp_guide_grob <- bp_legend_grob$grobs[[bp_guide_index]]

ggsave(
  file.path(plot_dir, "GO_BP_dotplot_legend.pdf"),
  bp_guide_grob, width = 15.5, height = 2.2, units = "in",
  device = grDevices::cairo_pdf, bg = "white"
)
ggsave(
  file.path(plot_dir, "GO_BP_dotplot_legend.png"),
  bp_guide_grob, width = 15.5, height = 2.2, units = "in",
  dpi = 600, bg = "white"
)

message("GO category summaries, BP dot plots, and the separate legend have been rendered.")
}))
