#!/usr/bin/env Rscript

# Select and plot KEGG pathway candidates from the shared final-model enrichment results.
# Read the significant KEGG table and cached names; query KEGGREST only for missing names.
# Write mapping/selection tables, kegg.pdf/png and separate legends under Main_plot/kegg/.
# Keep adjusted P < 0.05 map pathways, apply non-plant exclusions and select up to three per module.
# Rank by increasing adjusted P then decreasing Count; plot without repeating enrichment.

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
  library(stringr)
  library(KEGGREST)
})

project_dir <- RELEASE_ROOT
publication_dir <- file.path(
  project_dir,
  "6.WGCNA/stage2_GO_kegg/0.5_22_2_50/Main_plot"
)
dir.create(file.path(RELEASE_ROOT, "6.WGCNA/stage2_GO_kegg/0.5_22_2_50/Main_plot/kegg"), recursive = TRUE, showWarnings = FALSE)
kegg_table <- file.path(
  project_dir,
  "6.WGCNA/stage2_GO_kegg/0.5_22_2_50/KEGG/KEGG_enrichment_significant.csv"
)

old_mapping <- file.path(
  publication_dir,
  "KEGG_pathway_ID_name_mapping.csv"
)

new_mapping <- file.path(
  publication_dir,
  "kegg/KEGG_pathway_ID_name_mapping_KEGGREST.csv"
)

candidate_table <- file.path(
  publication_dir,
  "kegg/MainFig7B_kegg_candidates.csv"
)

selected_table <- file.path(
  publication_dir,
  "kegg/MainFig7B_kegg_terms.csv"
)

selection_audit <- file.path(
  publication_dir,
  "kegg/MainFig7B_kegg_audit.txt"
)

excluded_table <- file.path(
  publication_dir,
  "kegg/KEGG_excluded_nonplant_pathways.csv"
)

module_order <- c(
  "MEgrey60", "MEpurple", "MEturquoise", "MEred",
  "MEgreen", "MEpink", "MEmagenta", "MEbrown"
)

if (!file.exists(kegg_table)) {
  stop("KEGG enrichment table not found: ", kegg_table)
}

kegg_all <- read.csv(kegg_table, check.names = FALSE, stringsAsFactors = FALSE)

required_columns <- c(
  "module", "source", "ID", "Description",
  "p.adjust", "GeneRatio", "Count"
)
missing_columns <- setdiff(required_columns, colnames(kegg_all))
if (length(missing_columns)) {
  stop("The KEGG table is missing required columns: ", paste(missing_columns, collapse = ", "))
}

parse_gene_ratio <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  frac_idx <- grepl(
    "^[[:space:]]*[0-9.]+[[:space:]]*/[[:space:]]*[0-9.]+[[:space:]]*$",
    x
  )
  if (any(frac_idx, na.rm = TRUE)) {
    parts <- strsplit(x[frac_idx], "/", fixed = TRUE)
    out[frac_idx] <- vapply(parts, function(z) {
      a <- suppressWarnings(as.numeric(trimws(z[1])))
      b <- suppressWarnings(as.numeric(trimws(z[2])))
      if (!is.finite(a) || !is.finite(b) || b == 0) return(NA_real_)
      a / b
    }, numeric(1))
  }
  out
}

kegg_pathway <- kegg_all %>%
  filter(
    source == "KEGG_Pathway",
    str_detect(as.character(ID), "^map[0-9]+$")
  ) %>%
  mutate(
    module = as.character(module),
    ID = as.character(ID),
    p.adjust = suppressWarnings(as.numeric(p.adjust)),
    Count = suppressWarnings(as.numeric(Count)),
    GeneRatio_numeric = parse_gene_ratio(GeneRatio)
  ) %>%
  filter(
    module %in% module_order,
    p.adjust < 0.05,
    is.finite(Count), Count > 0,
    is.finite(GeneRatio_numeric), GeneRatio_numeric >= 0
  ) %>%
  arrange(module, p.adjust, desc(Count)) %>%
  distinct(module, ID, .keep_all = TRUE)

if (!nrow(kegg_pathway)) {
  stop("No mapXXXX pathway entries satisfy p.adjust < 0.05.")
}

mapping <- data.frame(
  pathway_id = character(),
  pathway_name = character(),
  mapping_status = character(),
  stringsAsFactors = FALSE
)

mapping_input <- if (file.exists(new_mapping)) new_mapping else old_mapping

if (file.exists(mapping_input)) {
  old <- read.csv(mapping_input, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("pathway_id", "pathway_name") %in% colnames(old))) {
    mapping <- old %>%
      transmute(
        pathway_id = as.character(pathway_id),
        pathway_name = as.character(pathway_name),
        mapping_status = if ("mapping_status" %in% colnames(old)) {
          as.character(mapping_status)
        } else {
          NA_character_
        }
      )
  }
}

mapping <- mapping %>%
  filter(!is.na(pathway_id)) %>%
  distinct(pathway_id, .keep_all = TRUE)

mapped_with_names <- mapping$pathway_id[
  !is.na(mapping$pathway_name) &
    str_trim(mapping$pathway_name) != ""
]
query_ids <- setdiff(unique(kegg_pathway$ID), mapped_with_names)

if (length(query_ids)) {
  queried <- lapply(query_ids, function(pathway_id) {
    ans <- tryCatch(KEGGREST::keggGet(pathway_id)[[1]], error = function(e) NULL)
    pathway_name <- NA_character_
    if (!is.null(ans) && !is.null(ans$NAME)) {
      pathway_name <- str_split(as.character(ans$NAME)[1], "[;\\n]", simplify = TRUE)[1]
      pathway_name <- str_trim(pathway_name)
    }
    data.frame(
      pathway_id = pathway_id,
      pathway_name = pathway_name,
      mapping_status = ifelse(
        is.na(pathway_name) || pathway_name == "",
        "KEGGREST_failed_or_empty",
        "KEGGREST_success"
      ),
      stringsAsFactors = FALSE
    )
  })
  mapping <- bind_rows(mapping, bind_rows(queried))
}

mapping <- mapping %>%
  mutate(
    has_name = !is.na(pathway_name) & str_trim(pathway_name) != ""
  ) %>%
  arrange(pathway_id, desc(has_name)) %>%
  distinct(pathway_id, .keep_all = TRUE) %>%
  dplyr::select(-has_name)

# Fixed local labels override retrieved names without changing enrichment statistics.
local_pathway_names <- c(
  map04075 = "Plant hormone signal transduction",
  map04933 = "AGE-RAGE signaling pathway in diabetic complications",
  map04670 = "Leukocyte transendothelial migration",
  map04120 = "Ubiquitin mediated proteolysis",
  map04140 = "Autophagy - animal",
  map04137 = "Mitophagy - animal",
  map03013 = "Nucleocytoplasmic transport",
  map00970 = "Aminoacyl-tRNA biosynthesis",
  map04110 = "Cell cycle",
  map00062 = "Fatty acid elongation",
  map00073 = "Cutin, suberine and wax biosynthesis",
  map00561 = "Glycerolipid metabolism",
  map04016 = "MAPK signaling pathway - plant",
  map04626 = "Plant-pathogen interaction",
  map00360 = "Phenylalanine metabolism"
)

local_idx <- match(names(local_pathway_names), mapping$pathway_id)
missing_local <- is.na(local_idx)
if (any(missing_local)) {
  mapping <- bind_rows(
    mapping,
    data.frame(
      pathway_id = names(local_pathway_names)[missing_local],
      pathway_name = unname(local_pathway_names[missing_local]),
      mapping_status = "local_fallback",
      stringsAsFactors = FALSE
    )
  )
  local_idx <- match(names(local_pathway_names), mapping$pathway_id)
}

mapping$pathway_name[local_idx] <- unname(local_pathway_names)
mapping$mapping_status[local_idx] <- "local_fallback"
write.csv(mapping, new_mapping, row.names = FALSE)

kegg_pathway <- kegg_pathway %>%
  left_join(mapping, by = c("ID" = "pathway_id")) %>%
  mutate(
    pathway_name = str_trim(pathway_name),
    display_name = case_when(
      !is.na(pathway_name) & pathway_name != "" ~ pathway_name,
      TRUE ~ ID
    )
  )

# The explicit list and keyword rule filter displayed candidates; upstream results remain complete.
excluded_kegg <- c(
  "Autophagy - animal",
  "Mitophagy - animal",
  "AGE-RAGE signaling pathway in diabetic complications",
  "Leukocyte transendothelial migration",
  "Viral myocarditis",
  "MAPK signaling pathway - fly",
  "MAPK signaling pathway",
  "Hedgehog signaling pathway - fly",
  "Osteoclast differentiation",
  "Fc epsilon RI signaling pathway",
  "Natural killer cell mediated cytotoxicity",
  "B cell receptor signaling pathway",
  "Toll-like receptor signaling pathway",
  "Toll and Imd signaling pathway",
  "NF-kappa B signaling pathway",
  "NOD-like receptor signaling pathway",
  "Axon guidance",
  "Neurotrophin signaling pathway",
  "Inflammatory mediator regulation of TRP channels",
  "Adipocytokine signaling pathway",
  "PPAR signaling pathway",
  "Ras signaling pathway",
  "Rap1 signaling pathway",
  "VEGF signaling pathway",
  "Wnt signaling pathway",
  "Bacterial invasion of epithelial cells",
  "Adherens junction",
  "Focal adhesion",
  "Salmonella infection",
  "Epithelial cell signaling in Helicobacter pylori infection",
  "Shigellosis",
  "Pertussis",
  "Toxoplasmosis",
  "Measles",
  "Leishmaniasis",
  "Chagas disease",
  "Tuberculosis"
)

animal_human_pattern <- paste(
  c(
    "- animal$",
    "- fly$",
    "- yeast$",
    "diabetic",
    "leukocyte",
    "cancer",
    "carcinoma",
    "carcinogenesis",
    "cardiomyopathy",
    "alzheimer",
    "parkinson",
    "huntington",
    "human immunodeficiency",
    "hepatitis",
    "influenza",
    "myocarditis",
    "chagas disease"
  ),
  collapse = "|"
)

kegg_pathway <- kegg_pathway %>%
  mutate(
    display_name_lower = str_to_lower(str_trim(display_name)),
    explicit_exclusion = display_name %in% excluded_kegg,
    regex_exclusion = str_detect(display_name_lower, animal_human_pattern),
    exclusion_reason = case_when(
      explicit_exclusion ~ "explicit non-plant pathway",
      regex_exclusion ~ "animal/human disease keyword",
      TRUE ~ NA_character_
    ),
    excluded_nonplant = explicit_exclusion | regex_exclusion
  )

excluded_pathways <- kegg_pathway %>%
  filter(excluded_nonplant) %>%
  dplyr::select(
    ID, Description, display_name, module,
    p.adjust, GeneRatio, Count, exclusion_reason
  )

write.csv(excluded_pathways, excluded_table, row.names = FALSE)

kegg_pathway <- kegg_pathway %>%
  filter(!excluded_nonplant) %>%
  dplyr::select(-display_name_lower, -explicit_exclusion,
                -regex_exclusion, -exclusion_reason,
                -excluded_nonplant)

kegg_pathway <- kegg_pathway %>%
  group_by(module) %>%
  arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
  mutate(
    module_rank = row_number(),
    selection_status = if_else(module_rank <= 3, "selected", "excluded_top3_limit")
  ) %>%
  ungroup()

write.csv(kegg_pathway, candidate_table, row.names = FALSE)

selected <- kegg_pathway %>%
  filter(selection_status == "selected") %>%
  mutate(module = factor(module, levels = module_order))

if (!nrow(selected)) stop("No pathways remain to plot after top-three selection.")

write.csv(selected, selected_table, row.names = FALSE)

audit_lines <- c(
  paste0("analysis_date: ", Sys.time()),
  paste0("R_version: ", R.version.string),
  paste0("input_table: ", kegg_table),
  paste0("mapping_output: ", new_mapping),
  paste0("candidate_table: ", candidate_table),
  paste0("selected_table: ", selected_table),
  paste0("excluded_table: ", excluded_table),
  "",
  "Filtering:",
  "  source == KEGG_Pathway",
  "  ID matches ^map[0-9]+$",
  "  ko-level entries excluded",
  "  p.adjust < 0.05",
  "  non-plant pathway filtering after enrichment and before representative selection",
  paste0("  excluded_pathway_rows: ", nrow(excluded_pathways)),
  "  maximum 3 pathways per module",
  "  ranking: p.adjust ascending, Count descending",
  "",
  paste0("candidate_rows: ", nrow(kegg_pathway)),
  paste0("selected_rows: ", nrow(selected)),
  paste0("selected_modules: ", paste(unique(as.character(selected$module)), collapse = ", ")),
  "",
  "sessionInfo:",
  capture.output(sessionInfo())
)
writeLines(audit_lines, selection_audit)

message("KEGG pathway mapping, filtering, and selection completed.")
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

selected_table <- file.path(
  publication_dir,
  "kegg/MainFig7B_kegg_terms.csv"
)

output_pdf <- file.path(
  publication_dir,
  "kegg/MainFig7B_kegg.pdf"
)

output_png <- file.path(
  publication_dir,
  "kegg/MainFig7B_kegg.png"
)

module_order <- c(
  "MEgrey60", "MEpurple", "MEturquoise", "MEred",
  "MEgreen", "MEpink", "MEmagenta", "MEbrown"
)

if (!file.exists(selected_table)) {
  stop("Selected KEGG table not found: ", selected_table)
}

plot_data <- read.csv(
  selected_table,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "module", "ID", "display_name", "p.adjust",
  "GeneRatio_numeric", "Count"
)

missing_columns <- setdiff(required_columns, colnames(plot_data))
if (length(missing_columns)) {
  stop("The input table is missing required columns: ", paste(missing_columns, collapse = ", "))
}

plot_data <- plot_data %>%
  mutate(
    module = as.character(module),
    display_name = as.character(display_name),
    p.adjust = suppressWarnings(as.numeric(p.adjust)),
    GeneRatio_numeric = suppressWarnings(as.numeric(GeneRatio_numeric)),
    Count = suppressWarnings(as.numeric(Count))
  ) %>%
  filter(
    module %in% module_order,
    !is.na(display_name), display_name != "",
    is.finite(p.adjust),
    is.finite(GeneRatio_numeric),
    is.finite(Count), Count > 0
  )

if (!nrow(plot_data)) stop("No KEGG pathway entries are available to plot.")

# Order pathways by best adjusted P, then by decreasing module count.
term_order <- plot_data %>%
  group_by(ID, display_name) %>%
  summarise(
    best_p = min(p.adjust, na.rm = TRUE),
    n_modules = n_distinct(module),
    .groups = "drop"
  ) %>%
  arrange(best_p, desc(n_modules), display_name) %>%
  pull(display_name)

plot_data <- plot_data %>%
  mutate(
    module = factor(module, levels = module_order),
    term = factor(display_name, levels = rev(term_order))
  )

ratio_breaks <- pretty(range(plot_data$GeneRatio_numeric, na.rm = TRUE), n = 4)
ratio_breaks <- ratio_breaks[
  ratio_breaks >= min(plot_data$GeneRatio_numeric, na.rm = TRUE) &
    ratio_breaks <= max(plot_data$GeneRatio_numeric, na.rm = TRUE)
]
if (length(ratio_breaks) < 2) ratio_breaks <- waiver()

p <- ggplot(plot_data, aes(x = module, y = term)) +
  geom_point(
    aes(size = GeneRatio_numeric, fill = p.adjust),
    shape = 21,
    colour = "#222222",
    stroke = 0.55,
    alpha = 0.96
  ) +
  scale_fill_gradient(
    low = "#E96B6B",
    high = "#3F86C3",
    name = "p.adjust",
    guide = guide_colorbar(
      order = 1,
      direction = "horizontal",
      reverse = TRUE,
      title.position = "left",
      title.hjust = 1,
      barwidth = grid::unit(42, "mm"),
      barheight = grid::unit(9.0, "mm"),
      frame.colour = NA,
      ticks.colour = "white"
    )
  ) +
  scale_size_continuous(
    name = "GeneRatio",
    range = c(3.2, 10.5),
    breaks = ratio_breaks,
    labels = label_number(accuracy = 0.01, trim = TRUE),
    guide = guide_legend(
      order = 2,
      direction = "horizontal",
      title.position = "left",
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
    limits = rev(term_order),
    labels = function(x) str_wrap(x, width = 38),
    expand = expansion(add = c(0.45, 0.45))
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_family = "Arial", base_size = 11) +
  theme(
    axis.text.x = element_text(
      family = "Arial", face = "plain", size = 16,
      angle = 45, hjust = 1, vjust = 1,
      colour = "black", margin = margin(t = 5)
    ),
    axis.text.y = element_text(
      family = "Arial", face = "plain", size = 16,
      hjust = 1, colour = "black", lineheight = 0.75,
      margin = margin(r = 6)
    ),
    axis.title.x = element_text(family = "Arial", face = "plain", size = 14),
    axis.title.y = element_text(family = "Arial", face = "plain", size = 14),
    panel.background = element_rect(fill = "#FAFAFA", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.42),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "#3F3F3F", fill = NA, linewidth = 0.65
    ),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.45),
    axis.ticks.length.y = grid::unit(2.2, "mm"),
    legend.position = "none",
    plot.margin = margin(8, 3, 6, 4)
  )

plot_height <- max(7.5, min(11.5, 1.5 + length(term_order) * 0.40))
plot_width <- 7.3

ggsave(
  output_pdf, p,
  width = plot_width, height = plot_height,
  units = "in", device = grDevices::cairo_pdf,
  limitsize = FALSE
)

ggsave(
  output_png, p,
  width = plot_width, height = plot_height,
  units = "in", dpi = 600,
  limitsize = FALSE
)

p_for_legend <- p +
  scale_fill_gradient(
    low = "#E96B6B",
    high = "#3F86C3",
    name = "p.adjust",
    breaks = scales::breaks_pretty(n = 3),
    labels = scales::label_number(accuracy = 0.001, trim = TRUE),
    guide = guide_colorbar(
      order = 1,
      direction = "horizontal",
      reverse = TRUE,
      title.position = "left",
      title.hjust = 1,
      barwidth = grid::unit(42, "mm"),
      barheight = grid::unit(9.0, "mm"),
      frame.colour = NA,
      ticks = FALSE
    )
  ) +
  theme(
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.title = element_text(
    family = "Arial",
    face = "plain",
    size = 17,
    colour = "black"
  ),
  legend.text = element_text(
    family = "Arial",
    face = "plain",
    size = 15,
    colour = "black"
  ),
  legend.key = element_blank(),
  legend.margin = margin(t = 1, r = 0, b = 1, l = 0),
  legend.spacing.x = grid::unit(4, "mm"),
  legend.box.spacing = grid::unit(2, "mm")
)

g <- ggplotGrob(p_for_legend)
legend_index <- which(g$layout$name == "guide-box-bottom")

if (length(legend_index) != 1) {
  stop("Cannot uniquely find bottom legend.")
}

legend_grob <- g$grobs[[legend_index]]

legend_pdf <- file.path(
  publication_dir,
  "kegg/MainFig7B_kegg_legend.pdf"
)

grDevices::cairo_pdf(
  legend_pdf,
  width = 12,
  height = 2.7
)
grid::grid.newpage()
grid::grid.draw(legend_grob)
dev.off()

legend_png <- file.path(
  publication_dir,
  "kegg/MainFig7B_kegg_legend.png"
)

png(
  legend_png,
  width = 12,
  height = 2.7,
  units = "in",
  res = 600
)
grid::grid.newpage()
grid::grid.draw(legend_grob)
dev.off()

message("KEGG plot-only outputs written to: ", publication_dir)
}))

