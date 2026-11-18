#!/usr/bin/env Rscript

# Select representative GO BP terms from shared enrichment tables and draw summary/S11-style plots.
# Read the final-model GO and KEGG CSVs; write selected terms, independent figures and legends.
# Keep up to three terms per module at adjusted P < 0.05, excluding ancestor/descendant redundancy.
# Draw selected terms without repeating enrichment or selection; keep empty modules in the summary.

.run_release_stage <- function(code) {
  prior_wd <- getwd()
  on.exit(setwd(prior_wd), add = TRUE)
  eval(code, envir = new.env(parent = globalenv()))
}

# Isolate selection and plotting so each stage can reset its own workspace.
.run_release_stage(quote({
RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(AnnotationDbi)
  library(GO.db)
})

project_dir <- RELEASE_ROOT
enrichment_dir <- file.path(
  project_dir,
  "6.WGCNA/stage2_GO_kegg/0.5_22_2_50"
)

go_table <- file.path(enrichment_dir, "GO/GO_enrichment_all_modules.csv")
kegg_table <- file.path(enrichment_dir, "KEGG/KEGG_enrichment_all_modules.csv")

publication_dir <- file.path(enrichment_dir, "Main_plot")
dir.create(publication_dir, recursive = TRUE, showWarnings = FALSE)

go_output_pdf <- file.path(publication_dir, "GO_BP_summary_plot.pdf")
go_output_csv <- file.path(publication_dir, "GO_BP_selected_terms.csv")
kegg_mapping_csv <- file.path(
  publication_dir,
  "KEGG_pathway_ID_name_mapping.csv"
)
kegg_table_named_csv <- file.path(
  publication_dir,
  "KEGG_pathway_enrichment_with_names.csv"
)

final_model <- "0.5_22_2_50"
model_id <- final_model
padj_cutoff <- 0.05
max_terms_per_module <- 3L

# KEGG name queries are disabled; the GO panel does not require the KEGG service.
run_kegg_mapping <- FALSE

module_order <- c(
  "MEgrey60",
  "MEpurple",
  "MEturquoise",
  "MEred",
  "MEgreen",
  "MEpink",
  "MEmagenta",
  "MEbrown"
)

if (!file.exists(go_table)) stop("GO enrichment table not found: ", go_table)
if (!file.exists(kegg_table)) stop("KEGG enrichment table not found: ", kegg_table)

go_all <- read.csv(go_table, check.names = FALSE, stringsAsFactors = FALSE)
kegg_all <- read.csv(kegg_table, check.names = FALSE, stringsAsFactors = FALSE)

required_go_columns <- c(
  "module", "source", "ID", "Description", "GeneRatio",
  "p.adjust", "Count"
)
required_kegg_columns <- required_go_columns

if (!all(required_go_columns %in% colnames(go_all))) {
  stop("The GO table is missing required columns: ", paste(
    setdiff(required_go_columns, colnames(go_all)), collapse = ", "
  ))
}
if (!all(required_kegg_columns %in% colnames(kegg_all))) {
  stop("The KEGG table is missing required columns: ", paste(
    setdiff(required_kegg_columns, colnames(kegg_all)), collapse = ", "
  ))
}

if ("final_model" %in% colnames(go_all)) {
  go_all <- go_all %>% filter(.data$final_model == model_id)
}
if ("final_model" %in% colnames(kegg_all)) {
  kegg_all <- kegg_all %>% filter(.data$final_model == model_id)
}

go_bp <- go_all %>%
  filter(source == "GO_BP", !is.na(p.adjust), p.adjust < padj_cutoff) %>%
  mutate(
    module = as.character(module),
    p.adjust = as.numeric(p.adjust),
    Count = as.numeric(Count),
    GeneRatio_numeric = {
      numerator <- suppressWarnings(as.numeric(sub("/.*$", "", GeneRatio)))
      denominator <- suppressWarnings(as.numeric(sub("^.*/", "", GeneRatio)))
      numerator / denominator
    }
  )

# Rank by adjusted P, raw P, decreasing Count and GeneRatio, then ID before removing redundancy.

get_bp_ancestors <- function(go_id) {
  ancestor_values <- as.list(GOBPANCESTOR)[[go_id]]
  if (is.null(ancestor_values)) return(character())
  unique(as.character(ancestor_values))
}

is_parent_child_related <- function(term_a, term_b) {
  if (identical(term_a, term_b)) return(TRUE)
  term_a %in% get_bp_ancestors(term_b) ||
    term_b %in% get_bp_ancestors(term_a)
}

select_representative_terms <- function(module_df) {
  if (!nrow(module_df)) {
    return(tibble(
      module = unique(module_df$module),
      selection_status = "NO_SIGNIFICANT_GO_BP",
      selection_reason = "No GO_BP term with adjusted P < 0.05",
      rank_by_padj = NA_integer_,
      selected_term = FALSE,
      redundant_with_term = NA_character_,
      ID = NA_character_,
      Description = NA_character_,
      p.adjust = NA_real_,
      GeneRatio = NA_character_,
      GeneRatio_numeric = NA_real_,
      Count = NA_real_
    ))
  }

  ordered <- module_df %>%
    arrange(p.adjust, pvalue, desc(Count), desc(GeneRatio_numeric), ID) %>%
    mutate(
      rank_by_padj = row_number(),
      selected_term = FALSE,
      selection_status = "REVIEW_CANDIDATE",
      selection_reason = "Not yet evaluated",
      redundant_with_term = NA_character_
    )

  selected_ids <- character()
  selected_rows <- integer()

  for (i in seq_len(nrow(ordered))) {
    current_id <- ordered$ID[i]

    if (length(selected_ids) >= max_terms_per_module) {
      ordered$selection_status[i] <- "NOT_SELECTED"
      ordered$selection_reason[i] <- "Module already has 3 selected terms"
      next
    }

    related_selected <- selected_ids[vapply(
      selected_ids,
      function(x) is_parent_child_related(current_id, x),
      logical(1)
    )]

    if (length(related_selected)) {
      ordered$selection_status[i] <- "NOT_SELECTED"
      ordered$selection_reason[i] <-
        "Redundant parent/child relationship with a more highly ranked term"
      ordered$redundant_with_term[i] <- related_selected[1]
    } else {
      ordered$selected_term[i] <- TRUE
      ordered$selection_status[i] <- "SELECTED"
      ordered$selection_reason[i] <-
        "Highest adjusted-P-ranked non-redundant term within module"
      selected_ids <- c(selected_ids, current_id)
      selected_rows <- c(selected_rows, i)
    }
  }

  ordered
}

selection_table <- bind_rows(lapply(
  module_order,
  function(module_name) {
    select_representative_terms(go_bp %>% filter(module == module_name))
  }
)) %>%
  mutate(module = factor(module, levels = module_order)) %>%
  arrange(module, rank_by_padj)

# Export all significant candidates, including terms rejected by redundancy or the three-term limit.
write.csv(
  selection_table %>% mutate(module = as.character(module)),
  go_output_csv,
  row.names = FALSE,
  na = ""
)

selected_terms <- selection_table %>%
  filter(selected_term == TRUE) %>%
  mutate(
    module = factor(as.character(module), levels = module_order),
    term_label = Description,
    neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
  )

missing_modules <- setdiff(module_order, as.character(selected_terms$module))
placeholder_terms <- tibble(
  module = factor(missing_modules, levels = module_order),
  term_label = "No significant GO_BP enrichment",
  GeneRatio_numeric = NA_real_,
  Count = NA_real_,
  neg_log10_padj = NA_real_,
  plot_placeholder = TRUE
)

plot_terms <- selected_terms %>%
  transmute(
    module,
    term_label,
    GeneRatio_numeric,
    Count,
    neg_log10_padj,
    plot_placeholder = FALSE
  ) %>%
  bind_rows(placeholder_terms) %>%
  mutate(module = factor(as.character(module), levels = module_order))

plot_terms$term_label <- factor(
  plot_terms$term_label,
  levels = rev(unique(plot_terms$term_label))
)

p <- ggplot(plot_terms) +
  geom_point(
    data = filter(plot_terms, !plot_placeholder),
    aes(
      x = GeneRatio_numeric,
      y = term_label,
      size = Count,
      color = neg_log10_padj
    )
  ) +
  geom_text(
    data = filter(plot_terms, plot_placeholder),
    aes(x = 0, y = term_label, label = term_label),
    size = 3.0,
    color = "grey45",
    hjust = 0
  ) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y", drop = FALSE) +
  scale_color_gradient(
    low = "#abdcff",
    high = "darkblue",
    name = "-log10 adjusted P",
    guide = guide_colorbar(
      frame.colour = "#F0EEE9",
      ticks.colour = "#C7C4BE",
      title.position = "top"
    )
  ) +
  scale_size_continuous(
    name = "Gene count",
    range = c(2.6, 8.0),
    guide = guide_legend(
      override.aes = list(color = "#657A95", alpha = 0.9)
    )
  ) +
  labs(
    x = "Gene ratio",
    y = NULL,
    title = paste0("GO Biological Process enrichment: ", final_model)
  ) +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.title.x = element_text(
      family = "Arial", face = "plain", size = 9
    ),
    axis.title.y = element_text(
      family = "Arial", face = "plain", size = 9
    ),
    axis.text.y = element_text(
      family = "Arial", face = "plain", size = 8.2, colour = "#222222"
    ),
    axis.text.x = element_text(
      family = "Arial", face = "plain", size = 8,
      colour = "#333333"
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
    strip.text.y = element_text(
      family = "Arial", face = "plain", size = 8, colour = "#333333"
    ),
    plot.title = element_text(
      family = "Arial", face = "plain", size = 10, hjust = 0.5
    ),
    legend.title = element_text(
      family = "Arial", face = "plain", size = 8
    ),
    legend.text = element_text(
      family = "Arial", face = "plain", size = 7
    ),
    legend.key = element_blank(),
    legend.position = "right"
  )

ggsave(
  go_output_pdf,
  p,
  width = 9,
  height = 10,
  units = "in",
  device = grDevices::cairo_pdf,
  limitsize = FALSE
)

# Optional name lookup preserves pathway IDs and records failed queries.
kegg_pathway_ids <- kegg_all %>%
  filter(source == "KEGG_Pathway") %>%
  distinct(ID) %>%
  filter(!is.na(ID), ID != "")

if (isTRUE(run_kegg_mapping) && requireNamespace("KEGGREST", quietly = TRUE)) {
  message("Starting KEGG pathway name mapping: ", nrow(kegg_pathway_ids), " IDs")
  mapping_rows <- lapply(seq_len(nrow(kegg_pathway_ids)), function(i) {
    pathway_id <- kegg_pathway_ids$ID[i]
    message("KEGG mapping ", i, "/", nrow(kegg_pathway_ids), ": ", pathway_id)
    query_id <- pathway_id
    result <- tryCatch(
      {
        setTimeLimit(elapsed = 15, transient = TRUE)
        KEGGREST::keggGet(query_id)
      },
      error = function(e) NULL,
      finally = setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    )

    if (length(result) && !is.null(result[[1]]$NAME)) {
      tibble(
        pathway_id = pathway_id,
        pathway_name = as.character(result[[1]]$NAME[1]),
        mapping_status = "KEGGREST_success"
      )
    } else {
      tibble(
        pathway_id = pathway_id,
        pathway_name = NA_character_,
        mapping_status = "KEGGREST_failed_or_unavailable"
      )
    }
  })
  kegg_mapping <- bind_rows(mapping_rows)
} else if (isTRUE(run_kegg_mapping)) {
  kegg_mapping <- kegg_pathway_ids %>%
    transmute(
      pathway_id = ID,
      pathway_name = NA_character_,
      mapping_status = "KEGGREST_package_unavailable"
    )
} else {
  message("KEGG pathway name mapping skipped; GO_BP outputs are complete.")
  kegg_mapping <- kegg_pathway_ids %>%
    transmute(
      pathway_id = ID,
      pathway_name = NA_character_,
      mapping_status = "KEGGREST_skipped_by_default"
    )
}

write.csv(kegg_mapping, kegg_mapping_csv, row.names = FALSE, na = "")

kegg_named <- kegg_all %>%
  filter(source == "KEGG_Pathway") %>%
  left_join(kegg_mapping, by = c("ID" = "pathway_id")) %>%
  mutate(
    display_name = ifelse(
      is.na(pathway_name) | pathway_name == "",
      Description,
      pathway_name
    )
  )
write.csv(kegg_named, kegg_table_named_csv, row.names = FALSE, na = "")

audit_lines <- c(
  paste0("final_model: ", final_model),
  paste0("GO_input_table: ", go_table),
  paste0("KEGG_input_table: ", kegg_table),
  paste0("module_order: ", paste(module_order, collapse = " > ")),
  "GO_main_figure_scope: GO_BP only",
  "GO_CC_and_GO_MF: excluded from Main Fig. 7B candidate",
  paste0("GO_padj_cutoff: ", padj_cutoff),
  paste0("max_terms_per_module: ", max_terms_per_module),
  "term_ranking: adjusted P, then raw P, Count and GeneRatio",
  "redundancy_rule: GO.db BP ancestor parent/child relationship",
  "empty_modules: retained as placeholder panels",
  paste0("KEGG_pathway_name_mapping_enabled: ", run_kegg_mapping),
  "enrichment_recalculated: FALSE",
  "existing_overview_plots_overwritten: FALSE",
  "",
  "sessionInfo:",
  capture.output(sessionInfo())
)
writeLines(audit_lines, file.path(publication_dir, "GO_BP_summary_plot_audit.txt"))

message("GO_BP publication-oriented summary plot and candidate table generated.")
message("Output directory: ", publication_dir)
}))

.run_release_stage(quote({
RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

project_dir <- RELEASE_ROOT
publication_dir <- file.path(
  project_dir,
  "6.WGCNA/stage2_GO_kegg/0.5_22_2_50/Main_plot"
)

selected_terms_file <- file.path(
  publication_dir,
  "GO_BP_selected_terms.csv"
)

output_pdf <- file.path(
  publication_dir,
  "GO_BP_figS11.pdf"
)
output_png <- file.path(
  publication_dir,
  "GO_BP_figS11.png"
)

module_order <- c(
  "MEgrey60",
  "MEpurple",
  "MEturquoise",
  "MEred",
  "MEgreen",
  "MEpink",
  "MEmagenta",
  "MEbrown"
)

if (!file.exists(selected_terms_file)) {
  stop("Selected-terms file not found: ", selected_terms_file)
}

selected_all <- read.csv(
  selected_terms_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "module", "selected_term", "ID", "Description",
  "p.adjust", "GeneRatio", "Count"
)

missing_columns <- setdiff(required_columns, colnames(selected_all))
if (length(missing_columns)) {
  stop("The selected-terms file is missing required columns: ", paste(missing_columns, collapse = ", "))
}

# Accept GeneRatio as either a fraction or a numeric value.
parse_gene_ratio <- function(x) {
  x <- as.character(x)

  out <- suppressWarnings(as.numeric(x))

  frac_idx <- grepl("^[[:space:]]*[0-9.]+[[:space:]]*/[[:space:]]*[0-9.]+[[:space:]]*$", x)
  if (any(frac_idx, na.rm = TRUE)) {
    parts <- strsplit(x[frac_idx], "/", fixed = TRUE)
    out[frac_idx] <- vapply(
      parts,
      function(z) {
        num <- suppressWarnings(as.numeric(trimws(z[1])))
        den <- suppressWarnings(as.numeric(trimws(z[2])))
        if (!is.finite(num) || !is.finite(den) || den == 0) return(NA_real_)
        num / den
      },
      numeric(1)
    )
  }

  out
}

selected <- selected_all %>%
  filter(selected_term == TRUE) %>%
  mutate(
    module = as.character(module),
    p.adjust = as.numeric(p.adjust),
    Count = as.numeric(Count),
    GeneRatio_numeric = parse_gene_ratio(GeneRatio)
  ) %>%
  filter(
    module %in% module_order,
    is.finite(p.adjust),
    p.adjust >= 0,
    is.finite(GeneRatio_numeric),
    GeneRatio_numeric >= 0
  )

if (!nrow(selected)) {
  stop("No selected GO BP terms are available to plot in selected_terms.csv.")
}

# Order terms by biological theme without changing the selected set.
term_order_table <- selected %>%
  distinct(ID, Description) %>%
  mutate(
    theme_rank = case_when(
      str_detect(
        str_to_lower(Description),
        "defense|immune|response to biotic|pathogen|chitin"
      ) ~ 1L,
      str_detect(
        str_to_lower(Description),
        "signal|hormone|auxin|jasmon|ethylene|abscisic"
      ) ~ 2L,
      str_detect(
        str_to_lower(Description),
        "protein modification|ubiquitin|proteol|catabol|amino acid"
      ) ~ 3L,
      str_detect(
        str_to_lower(Description),
        "rna|gene expression|transcription|translation|ribosom"
      ) ~ 4L,
      str_detect(
        str_to_lower(Description),
        "cell wall|lipid|cutin|suberin|phenylpropanoid|secondary"
      ) ~ 5L,
      str_detect(
        str_to_lower(Description),
        "development|transport|meristem|morphogenesis"
      ) ~ 6L,
      TRUE ~ 7L
    ),
    original_order = row_number()
  ) %>%
  arrange(theme_rank, original_order)

term_levels <- term_order_table$Description
n_terms <- length(term_levels)

plot_data <- selected %>%
  mutate(
    module = factor(module, levels = module_order),
    term = factor(Description, levels = rev(term_levels))
  )

# Keep four readable legend breaks even when GeneRatio spans a narrow range.
ratio_breaks <- pretty(range(plot_data$GeneRatio_numeric, na.rm = TRUE), n = 4)
ratio_breaks <- ratio_breaks[
  ratio_breaks >= min(plot_data$GeneRatio_numeric, na.rm = TRUE) &
  ratio_breaks <= max(plot_data$GeneRatio_numeric, na.rm = TRUE)
]
if (length(ratio_breaks) < 2) ratio_breaks <- waiver()

p <- ggplot(
  plot_data,
  aes(x = module, y = term)
) +
  geom_point(
    aes(
      size = GeneRatio_numeric,
      fill = p.adjust
    ),
    shape = 21,
    colour = "#222222",
    stroke = 0.55,
    alpha = 0.96
  ) +

  # Map smaller adjusted P values to red and reverse the horizontal legend.
  scale_fill_gradient(
    low = "#E96B6B",
    high = "#3F86C3",
    name = "p.adjust",
    guide = guide_colorbar(
      order = 1,
      direction = "vertical",
      reverse = TRUE,
      title.position = "top",
      title.hjust = 1,
      barwidth = grid::unit(9, "mm"),
      barheight = grid::unit(42, "mm"),
      frame.colour = NA,
      ticks.colour = "white"
    )
  ) +
  scale_size_continuous(
    name = "GeneRatio",
    range = c(4, 12),
    breaks = ratio_breaks,
    labels = label_number(
      accuracy = 0.01,
      trim = TRUE
    ),
    guide = guide_legend(
      order = 2,
      direction = "vertical",
      title.position = "top",
      title.hjust = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        colour = "#222222",
        alpha = 1,
        stroke = 0.5
      )
    )
  ) +

  scale_x_discrete(
    limits = module_order,
    drop = FALSE,
    expand = expansion(add = c(0.45, 0.45))
  ) +
  scale_y_discrete(
    limits = rev(term_levels),
    labels = function(x) str_wrap(x, width = 38),
    expand = expansion(add = c(0.45, 0.45))
  ) +

  labs(x = NULL, y = NULL) +

  theme_bw(base_family = "Arial", base_size = 18) +
  theme(
    axis.text.x = element_text(
      family = "Arial",
      face = "plain",
      size = 18,
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black",
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(
      family = "Arial",
      face = "plain",
      size = 18,
      hjust = 1,
      colour = "black",
      lineheight = 0.75,
      margin = margin(r = 6)
    ),
    axis.title.x = element_text(
      family = "Arial", face = "plain", size = 18
    ),
    axis.title.y = element_text(
      family = "Arial", face = "plain", size = 18
    ),

    panel.background = element_rect(fill = "#FAFAFA", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(
      colour = "#E5E5E5",
      linewidth = 0.42
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "#3F3F3F",
      fill = NA,
      linewidth = 0.65
    ),

    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(
      colour = "black",
      linewidth = 0.45
    ),
    axis.ticks.length.y = grid::unit(2.2, "mm"),

    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.justification = "left",
    legend.direction = "vertical",
    legend.title = element_text(family = "Arial",face = "plain",size = 17,colour = "black"),
    legend.text = element_text(family = "Arial",face = "plain",size = 17,colour = "black"),
    legend.key = element_blank(),
    legend.margin = margin(t = 1, r = 0, b = 1, l = 0),
    legend.spacing.x = grid::unit(4, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),
    legend.key.width = grid::unit(8, "mm"),
    plot.margin = margin(8, 3, 6, 4)
  )

plot_height <- max(7.5, min(12.5, 2 + n_terms * 0.80))
plot_width <- 10.5

ggsave(
  filename = output_pdf,
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  device = grDevices::cairo_pdf,
  limitsize = FALSE
)

ggsave(
  filename = output_png,
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)

message("Reference-style Main Fig. 7B written to: ", publication_dir)
message("Existing selected terms table was read only: ", selected_terms_file)

p_for_legend <- p + theme(
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.title = element_text(
    family = "Arial",
    face = "plain",
    size = 14,
    colour = "black"
  ),
  legend.text = element_text(
    family = "Arial",
    face = "plain",
    size = 12,
    colour = "black"
  ),
)

g <- ggplotGrob(p_for_legend)

legend_index <- which(g$layout$name == "guide-box-bottom")

if (length(legend_index) != 1) {
  stop("Cannot uniquely find bottom legend.")
}

legend_grob <- g$grobs[[legend_index]]

legend_pdf <- file.path(
  publication_dir,
  "GO_BP_figS11_legend.pdf"
)

grDevices::cairo_pdf(
  legend_pdf,
  width = 10,
  height = 1.3
)

grid::grid.newpage()
grid::grid.draw(legend_grob)

dev.off()

legend_png <- file.path(
  publication_dir,
  "GO_BP_MainFig7B_legend.png"
)

png(
  legend_png,
  width = 10,
  height = 1.3,
  units = "in",
  res = 600
)

grid::grid.newpage()
grid::grid.draw(legend_grob)

dev.off()

audit_lines <- c(
  paste0("analysis_date: ", Sys.time()),
  paste0("R_version: ", R.version.string),
  paste0("working_directory: ", getwd()),
  paste0("input_selected_terms: ", selected_terms_file),
  paste0("output_pdf: ", output_pdf),
  paste0("output_png: ", output_png),
  "",
  "Display mapping:",
  "  X = module",
  "  Y = selected GO_BP term",
  "  size = GeneRatio",
  "  fill = p.adjust (small = red, large = blue)",
  "  point border = black",
  "",
  "sessionInfo:",
  capture.output(sessionInfo())
)

writeLines(
  audit_lines,
  file.path(publication_dir, "analysis_audit_reference_style.txt")
)
}))

