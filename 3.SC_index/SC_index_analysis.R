#!/usr/bin/env Rscript

# Analyze subgenome contribution and response amplitude from normalized Cm scaledTPM
# and curated OG annotations. Write Tables S7/S8, Fig. 4, and Supplementary Figs. 7/8.
# cHEB uses (B-S)/(B+S); amplitudes summarize the 0/24/48 hpi expression trajectory.
# These statistics use within-species expression, without VST or cross-species factors.

# Keep analysis and final-figure variables separate while preserving the working directory.
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
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(ggrepel)
})

project_dir <- RELEASE_ROOT
normalization_dir <- file.path(project_dir, "1.normalization")
normalization_result_dir <- file.path(
  normalization_dir, "results", "01_tximport", "Cm"
)

og_path <- file.path(project_dir, "OG_info", "OG_info_260906_OK.xlsx")
normalized_counts_path <- file.path(
  normalization_result_dir,
  "Cm_scaledTPM_withinSpecies_normalized.tsv"
)

output_dir <- file.path(project_dir, "3.SC_index")
table_dir <- file.path(output_dir, "tables")
main_figure_dir <- file.path(output_dir, "Main_Fig_4")
supp_figure_7_dir <- file.path(output_dir, "Supp_Fig_7")
supp_figure_8_dir <- file.path(output_dir, "Supp_Fig_8")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_figure_7_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_figure_8_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(og_path, normalized_counts_path)
if (any(!file.exists(required_files))) {
  stop(
    "Missing input files: ",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
}

species_prefix <- "Cm"
timepoint_levels <- c("0", "24", "48")
module_order <- c(
  "PRRs", "RLCKs", "NLRs", "MAPKs", "WRKYs",
  "Ca2+", "SAs", "JAs & ETs", "Cucurbitacins"
)

# Eight final Cm libraries are retained; Cm_48h_D2 is excluded.
expected_samples <- c(
  "Cm_0h_D1", "Cm_0h_D2", "Cm_0h_D3",
  "Cm_24h_D1", "Cm_24h_D2", "Cm_24h_D3",
  "Cm_48h_D1", "Cm_48h_D3"
)

write_audit <- function(name, data) {
  write_tsv(data, file.path(output_dir, name), na = "NA")
}

read_expression_matrix <- function(path, expected_sample_names) {
  x <- read_tsv(path, show_col_types = FALSE)
  if (!"GID" %in% names(x)) {
    stop("Expression matrix must contain a GID column: ", path)
  }
  if (anyDuplicated(x$GID)) stop("Expression matrix contains duplicate GIDs: ", path)

  sample_names <- setdiff(names(x), "GID")
  if (!setequal(sample_names, expected_sample_names)) {
    stop(
      "Expression samples do not match the final Cm sample set. Observed: ",
      paste(sample_names, collapse = ", "),
      "; expected: ", paste(expected_sample_names, collapse = ", ")
    )
  }
  if (any(!is.finite(as.matrix(x[, sample_names])))) {
    stop("Expression matrix contains NA/NaN/Inf: ", path)
  }
  if (any(as.matrix(x[, sample_names]) < 0)) {
    stop("Expression matrix contains negative values: ", path)
  }

  mat <- as.matrix(x[, expected_sample_names])
  rownames(mat) <- x$GID
  mat
}

aggregate_replicates <- function(mat, prefix, timepoints) {
  sample_info <- tibble(sample = colnames(mat)) |>
    mutate(
      tp = str_match(sample, paste0("^", prefix, "_(\\d+)h_D\\d+$"))[, 2],
      rep = str_match(sample, paste0("^", prefix, "_\\d+h_D(\\d+)$"))[, 2]
    )

  if (anyNA(sample_info$tp) || anyNA(sample_info$rep)) {
    stop("Expression sample names cannot be parsed into timepoint and replicate.")
  }
  if (!setequal(unique(sample_info$tp), timepoints)) {
    stop("Expression matrix is missing expected timepoints: ", paste(timepoints, collapse = ", "))
  }

  split_columns <- split(
    sample_info$sample,
    factor(sample_info$tp, levels = timepoints)
  )
  result <- vapply(split_columns, function(cols) {
    x <- mat[, cols, drop = FALSE]
    values <- rowMeans(x, na.rm = TRUE)
    # A gene with no positive replicate at this timepoint does not count toward n_B/n_S.
    values[rowSums(x > 0, na.rm = TRUE) < 1] <- NA_real_
    values
  }, FUN.VALUE = numeric(nrow(mat)))
  rownames(result) <- rownames(mat)
  result
}

make_og_time_table <- function(aggregated, og_map) {
  result <- as_tibble(aggregated, rownames = "GID") |>
    pivot_longer(-GID, names_to = "tp", values_to = "expr") |>
    inner_join(og_map, by = "GID") |>
    group_by(OG_name, Module, Type, tp, subg) |>
    summarise(
      E_sum = sum(expr, na.rm = TRUE),
      n = sum(is.finite(expr)),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = subg,
      values_from = c(E_sum, n),
      values_fill = 0,
      names_sep = "_"
    )

  for (nm in c("E_sum_B", "E_sum_S", "n_B", "n_S")) {
    if (!nm %in% names(result)) result[[nm]] <- 0
  }

  result |>
    mutate(tp = factor(tp, levels = timepoint_levels)) |>
    arrange(OG_name, tp)
}

# Absolute amplitude sums absolute successive changes; directional amplitude is final minus initial.
path_length <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < 2L) return(NA_real_)
  sum(abs(diff(values)))
}

net_change <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < 2L) return(NA_real_)
  tail(values, 1) - head(values, 1)
}

classify_origin <- function(type) {
  type <- toupper(as.character(type))
  case_when(
    str_detect(type, "BS") ~ "BS_pair",
    type %in% c("B", "BB") ~ "B_only",
    type %in% c("S", "SS") ~ "S_only",
    TRUE ~ "Other"
  )
}

og_info <- read_excel(og_path, sheet = 1)
required_og_columns <- c(
  "Cmo", "GeneName_base", "Module", "subgenome", "Type"
)
if (!all(required_og_columns %in% names(og_info))) {
  stop(
    "OG_info_260906_OK.xlsx is missing required columns: ",
    paste(setdiff(required_og_columns, names(og_info)), collapse = ", ")
  )
}

og_map <- og_info |>
  transmute(
    GID = as.character(Cmo),
    OG_name = as.character(GeneName_base),
    Module = if_else(as.character(Module) == "JAs", "JAs & ETs", as.character(Module)),
    Type = as.character(Type),
    subg = case_when(
      as.character(subgenome) %in% c("CmB", "B") ~ "B",
      as.character(subgenome) %in% c("CmS", "S") ~ "S",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(GID), GID != "", !is.na(OG_name), !is.na(subg)) |>
  distinct(GID, OG_name, Module, Type, subg)

if (anyDuplicated(og_map$GID)) {
  stop("A Cm GID is assigned to multiple OG/subgenome records in OG_info.")
}

write_audit(
  "input_og_mapping_audit.tsv",
  tibble(
    input_file = basename(og_path),
    n_rows_after_mapping = nrow(og_map),
    n_unique_gene = n_distinct(og_map$GID),
    n_unique_og = n_distinct(og_map$OG_name),
    n_B = sum(og_map$subg == "B"),
    n_S = sum(og_map$subg == "S"),
    duplicate_gene_check = ifelse(anyDuplicated(og_map$GID), "FAIL", "PASS")
  )
)

normalized_mat <- read_expression_matrix(
  normalized_counts_path,
  expected_samples
)

sample_audit <- tibble(sample = colnames(normalized_mat)) |>
  mutate(
    species = str_match(sample, "^([A-Za-z]+)_")[, 2],
    timepoint = str_match(sample, "^[A-Za-z]+_(\\d+h)_")[, 2],
    replicate = str_match(sample, "_D(\\d+)$")[, 2]
  )
write_audit("final_Cm_sample_audit.tsv", sample_audit)

write_audit(
  "SC_index_input_audit.tsv",
  tibble(
    input_file = basename(normalized_counts_path),
    input_path = normalized_counts_path,
    n_gene = nrow(normalized_mat),
    n_library = ncol(normalized_mat),
    sample_set = paste(colnames(normalized_mat), collapse = ";"),
    excluded_samples = "Cm_48h_D2; Ls_24h_D1 (Ls is not included in this Cm-only SC-index matrix)",
    scaled_TPM_within_species_normalized = "PASS",
    VST_used_for_amplitude = "NO",
    finite_check = "PASS",
    nonnegative_check = "PASS"
  )
)

og_map_normalized <- og_map |>
  filter(GID %in% rownames(normalized_mat))
normalized_time <- aggregate_replicates(
  normalized_mat, species_prefix, timepoint_levels
)

sc_analysis <- make_og_time_table(normalized_time, og_map_normalized) |>
  mutate(
    SC_index = if_else(
      E_sum_B + E_sum_S > 0,
      (E_sum_B - E_sum_S) / (E_sum_B + E_sum_S),
      NA_real_
    ),
    tp = as.character(tp)
  )

if (any(!is.finite(sc_analysis$E_sum_B)) || any(!is.finite(sc_analysis$E_sum_S))) {
  stop("OG-level cHEB input contains non-finite expression values.")
}

supp_table_s7 <- sc_analysis |>
  transmute(
    OG = OG_name,
    Timepoint = tp,
    E_sum_B,
    E_sum_S,
    n_B,
    n_S,
    cHEB = SC_index
  )
write_csv(supp_table_s7, file.path(table_dir, "SC_index_by_OG_time_CalMean.csv"))

amplitude_time <- sc_analysis |>
  mutate(
    tp_num = as.integer(tp),
    E_total = E_sum_B + E_sum_S
  ) |>
  arrange(OG_name, tp_num)

og_type <- og_map_normalized |>
  distinct(OG_name, Module, Type)
if (anyDuplicated(og_type[c("OG_name", "Module")])) {
  stop("An OG/module has multiple Type assignments.")
}

amplitude_analysis <- amplitude_time |>
  group_by(OG_name, Module) |>
  summarise(
    A_pathlen_OG = path_length(E_total),
    A_pathlen_B = path_length(E_sum_B),
    A_pathlen_S = path_length(E_sum_S),
    A_dir_OG = net_change(E_total),
    A_dir_B = net_change(E_sum_B),
    A_dir_S = net_change(E_sum_S),
    .groups = "drop"
  ) |>
  left_join(og_type, by = c("OG_name", "Module")) |>
  mutate(
    transient_B = abs(A_dir_B) < A_pathlen_B,
    transient_S = abs(A_dir_S) < A_pathlen_S,
    Transient = case_when(
      transient_B & transient_S ~ "Both",
      transient_B & !transient_S ~ "B",
      !transient_B & transient_S ~ "S",
      TRUE ~ "None"
    )
  )

supp_table_s8 <- amplitude_analysis |>
  transmute(
    OG = OG_name,
    Aabs_OG = A_pathlen_OG,
    Aabs_B = A_pathlen_B,
    Aabs_S = A_pathlen_S,
    Adir_OG = A_dir_OG,
    Adir_B = A_dir_B,
    Adir_S = A_dir_S,
    Transient
  )
write_csv(
  supp_table_s8,
  file.path(table_dir, "OG_summary_A_CalMean_scaledTPM_normalized.csv")
)

write_audit(
  "SC_index_amplitude_audit.tsv",
  tibble(
    analysis = c("cHEB/SC-index", "response amplitude"),
    input = c("Cm_scaledTPM_withinSpecies_normalized.tsv",
              "Cm_scaledTPM_withinSpecies_normalized.tsv"),
    cross_species_scaling_applied = c("NO", "NO"),
    VST_applied = c("NO", "NO"),
    method = c("centered HEB: (B-S)/(B+S)", "absolute path length across 0/24/48 hpi")
  )
)

sc_main <- sc_analysis |>
  mutate(
    tp = if_else(tp == "00", "0", tp),
    Module = factor(Module, levels = module_order),
    Type = as.character(Type)
  )

amplitude_main <- amplitude_analysis |>
  mutate(
    Module = factor(Module, levels = module_order),
    Type = as.character(Type)
  )

module_response <- amplitude_main |>
  mutate(origin_class = classify_origin(Type)) |>
  group_by(Module, origin_class) |>
  summarise(
    Aabs_B_sum = sum(A_pathlen_B, na.rm = TRUE),
    Aabs_S_sum = sum(A_pathlen_S, na.rm = TRUE),
    Aabs_mean = mean(A_pathlen_OG, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    B_share = Aabs_B_sum / (Aabs_B_sum + Aabs_S_sum + 1e-6),
    B_share_centered = 2 * B_share - 1
  )

size_limits <- range(module_response$Aabs_mean, na.rm = TRUE)
if (!all(is.finite(size_limits))) size_limits <- c(0, 1)
if (size_limits[1] == size_limits[2]) size_limits <- size_limits + c(-1, 1) * 1e-6

theme_card_strip <- theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, vjust = 3),
    plot.margin = margin(t = 20, r = 4, b = 2, l = 4)
  )

theme_card_mid <- theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 2, r = 15, b = 8, l = 15)
  )

theme_card_bottom <- theme_grey(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 4, b = 2, l = 4)
  )

make_composition_panel <- function(module) {
  data <- amplitude_main |>
    filter(Module == module) |>
    mutate(origin_class = classify_origin(Type)) |>
    filter(origin_class != "Other") |>
    count(origin_class, name = "n") |>
    mutate(
      prop = n / sum(n),
      origin_class = factor(origin_class, levels = c("B_only", "BS_pair", "S_only"))
    )

  ggplot(data, aes(x = 1, y = prop, fill = origin_class)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = c(
      S_only = "#ee6055", BS_pair = "#f1a93b", B_only = "#1e88e5"
    )) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = as.character(module)) +
    theme_card_strip
}

make_density_panel <- function(module) {
  data <- sc_main |>
    filter(Module == module, Type %in% c("BS", "BS_pair")) |>
    filter(!is.na(SC_index), between(SC_index, -1, 1))

  ggplot(data, aes(x = SC_index, color = tp, group = tp)) +
    geom_density(linewidth = 0.85, adjust = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    scale_color_manual(values = c("0" = "#8ef0b4", "24" = "#4abf8a", "48" = "#3b7e9b")) +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(x = "cHEB (BS pairs)", y = "Density") +
    theme_card_mid
}

make_response_panel <- function(module) {
  data <- module_response |> filter(Module == module)

  ggplot(data, aes(x = B_share_centered, y = 0, group = origin_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_point(
      aes(size = Aabs_mean, fill = B_share_centered),
      shape = 21, color = "grey50", stroke = 0.25
    ) +
    coord_cartesian(xlim = c(-1.1, 1.1)) +
    scale_fill_gradient2(
      midpoint = 0, limits = c(-0.75, 0.75),
      low = "#e63746", high = "#287fb8", oob = scales::squish,
      breaks = c(-0.5, 0, 0.5), name = "Response bias"
    ) +
    scale_size_continuous(
      limits = size_limits, range = c(1.5, 8.5), name = "Aabs mean"
    ) +
    labs(x = "Response bias", y = NULL) +
    theme_card_bottom
}

make_module_card <- function(module) {
  density_panel <- make_density_panel(module)
  response_panel <- make_response_panel(module)
  (make_composition_panel(module) /
    density_panel /
    response_panel /
    plot_spacer()) +
    plot_layout(heights = c(0.18, 1, 0.35, 0.05))
}

modules_main <- module_order[module_order %in% as.character(unique(sc_main$Module))]
main_figure <- wrap_plots(
  lapply(modules_main, make_module_card), ncol = 3
) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

ggsave(
  file.path(main_figure_dir, "SC_summary_new_addSingleA_meanVersion.pdf"),
  main_figure, width = 8.5, height = 11, units = "in"
)

# Auxiliary module panels are retained; the final Fig. S7 composite is generated below.
pal_B <- brewer.pal(3, "Blues")
pal_S <- brewer.pal(3, "Reds")
pal_BS <- brewer.pal(3, "Greens")
color_levels <- c(
  "single_B_0", "single_B_24", "single_B_48",
  "single_S_0", "single_S_24", "single_S_48",
  "BS_0", "BS_24", "BS_48"
)
color_map <- setNames(c(pal_B, pal_S, pal_BS), color_levels)

sc_supp7 <- sc_analysis |>
  mutate(
    tp = factor(tp, levels = timepoint_levels),
    Type_group = case_when(
      Type %in% c("B", "BB") ~ "single_B",
      Type %in% c("S", "SS") ~ "single_S",
      TRUE ~ "BS"
    ),
    expr_total = E_sum_B + E_sum_S,
    expr_total = if_else(expr_total > 0, expr_total, NA_real_),
    color_key = paste0(Type_group, "_", as.character(tp))
  ) |>
  filter(!is.na(expr_total))

baseline_epsilon <- 0.5
baseline_table <- sc_supp7 |>
  group_by(Module, OG_name) |>
  summarise(
    baseline_raw = exp(mean(log(expr_total + baseline_epsilon), na.rm = TRUE)),
    .groups = "drop"
  )

sc_supp7 <- sc_supp7 |>
  left_join(baseline_table, by = c("Module", "OG_name")) |>
  mutate(
    ratio = (expr_total + baseline_epsilon) / (baseline_raw + baseline_epsilon),
    clipped_change = pmin(pmax(log2(ratio), -1.5), 1.5),
    size_final = pmin(pmax(3 * 2^(0.60 * clipped_change), 1.8), 6)
  )

adir_long <- amplitude_analysis |>
  select(OG_name, Module, Adir_B = A_dir_B, Adir_S = A_dir_S) |>
  pivot_longer(c(Adir_B, Adir_S), names_to = "Subgenome", values_to = "Adir") |>
  mutate(Subgenome = recode(Subgenome, Adir_B = "B", Adir_S = "S")) |>
  left_join(
    sc_supp7 |> distinct(Module, OG_name, Type_group),
    by = c("Module", "OG_name")
  )

fill_colors <- c(B = "#a4c9df", S = "#f49fbe")
modules_supp7 <- sort(intersect(unique(adir_long$Module), unique(sc_supp7$Module)))

for (module in modules_supp7) {
  module_data <- sc_supp7 |> filter(Module == module)

  left_plot <- ggplot(
    adir_long |> filter(Module == module),
    aes(x = OG_name, y = Adir, fill = Subgenome)
  ) +
    geom_col(width = 0.8, position = "stack") +
    scale_fill_manual(values = fill_colors, name = "Subgenome") +
    coord_flip() +
    labs(x = NULL, y = "Adir", title = module) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")

  right_plot <- ggplot(module_data, aes(x = SC_index, y = OG_name)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_point(aes(size = size_final, colour = color_key), alpha = 0.92, stroke = 0) +
    scale_size_identity() +
    scale_colour_manual(values = color_map, guide = "none") +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(x = "cHEB", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          panel.grid.minor = element_blank())

  combined_plot <- left_plot + right_plot + plot_layout(widths = c(1.2, 1))
  safe_module_name <- gsub("[\\/:*?\"<>|]", "_", as.character(module))
  ggsave(
    file.path(supp_figure_7_dir, paste0(safe_module_name, "_combo.pdf")),
    combined_plot, width = 8.5, height = max(2.5, 0.22 * nrow(module_data) + 1.5),
    units = "in", dpi = 300
  )
}

amplitude_supp8 <- amplitude_analysis |>
  mutate(
    Transient = factor(Transient, levels = c("Both", "B", "S", "None")),
    shape = Type == "BS"
  )

label_data <- amplitude_supp8 |>
  group_by(Transient) |>
  mutate(distance_from_origin = sqrt(A_dir_B^2 + A_dir_S^2)) |>
  slice_max(distance_from_origin, n = 5, with_ties = FALSE) |>
  ungroup()

supp_figure_8 <- ggplot(amplitude_supp8, aes(x = A_dir_B, y = A_dir_S)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  geom_point(
    aes(size = A_pathlen_OG, color = Transient, shape = shape),
    stroke = 1, alpha = 0.65
  ) +
  geom_text_repel(
    data = label_data, aes(label = OG_name),
    size = 3.2, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.2
  ) +
  scale_color_manual(values = c(
    Both = "orange", B = "#1e88e5", S = "#ee6055", None = "#5fd8a4"
  )) +
  scale_size_continuous(name = expression(A^{abs}[OG]), range = c(1, 12)) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  labs(x = "Adir(B)", y = "Adir(S)") +
  theme_classic() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 11)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(~Transient)

ggsave(
  file.path(supp_figure_8_dir, "A_level_xyplot_facet_label.pdf"),
  supp_figure_8, width = 12, height = 8, units = "in"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "SC_index_sessionInfo.txt")
)

message("SC-index/cHEB and response-amplitude analysis completed.")
}))

# Render final Fig. S7 from the exported S7/S8 tables without recomputing statistics.
.run_release_stage(quote({

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

project_dir <- RELEASE_ROOT
if (!dir.exists(project_dir)) stop("Project directory does not exist: ", project_dir)

sc_dir <- file.path(project_dir, "3.SC_index")
table_dir <- file.path(sc_dir, "tables")
og_path <- file.path(project_dir, "OG_info", "OG_info_260906_OK.xlsx")

s7_path <- file.path(table_dir, "SC_index_by_OG_time_CalMean.csv")
s8_path <- file.path(
  table_dir,
  "OG_summary_A_CalMean_scaledTPM_normalized.csv"
)

output_dir <- file.path(
  sc_dir,
  "Supp_Fig_7",
  "results"
)
individual_output_dir <- file.path(output_dir, "individual_modules")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(og_path, s7_path, s8_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Missing Supplementary Fig. 7 plotting inputs: \n",
    paste(missing_files, collapse = "\n")
  )
}

module_order <- c(
  "PRRs", "RLCKs", "NLRs", "MAPKs", "WRKYs",
  "Ca2+", "SAs", "JAs & ETs", "Cucurbitacins"
)

timepoint_levels <- c("0", "24", "48")

pal_B <- brewer.pal(3, "Blues")
pal_S <- brewer.pal(3, "Reds")
pal_BS <- brewer.pal(3, "Greens")
color_levels <- c(
  "single_B_0", "single_B_24", "single_B_48",
  "single_S_0", "single_S_24", "single_S_48",
  "BS_0", "BS_24", "BS_48"
)
color_map <- setNames(c(pal_B, pal_S, pal_BS), color_levels)
fill_colors <- c(B = "#a4c9df", S = "#f49fbe")

# These within-module baseline and winsorization settings control point size only.
winsor_lo <- 0.05
winsor_hi <- 0.95
size_min <- 1.8
size_max <- 6
baseline_epsilon <- 0.5
clip_change <- 1.5
change_weight <- 0.60
row_height <- 0.22
module_fixed_height <- 0.8
left_width <- 3.6
middle_width <- 1
right_width <- 3.2
plot_base_size <- 11
plot_base_family <- "Helvetica"
og_label_size <- 3

s7 <- read_csv(s7_path, show_col_types = FALSE) |>
  transmute(
    OG_name = as.character(OG),
    tp = as.character(as.integer(Timepoint)),
    E_sum_B = as.numeric(E_sum_B),
    E_sum_S = as.numeric(E_sum_S),
    SC_index = as.numeric(cHEB)
  ) |>
  mutate(tp = factor(tp, levels = timepoint_levels))

s8 <- read_csv(s8_path, show_col_types = FALSE) |>
  transmute(
    OG_name = as.character(OG),
    Adir_B = as.numeric(Adir_B),
    Adir_S = as.numeric(Adir_S)
  )

og_info <- read_excel(og_path, sheet = 1)
required_og_columns <- c(
  "GeneName_base", "Module", "Type"
)
if (!all(required_og_columns %in% names(og_info))) {
  stop(
    "OG_info is missing required columns: ",
    paste(setdiff(required_og_columns, names(og_info)), collapse = ", ")
  )
}

og_annotation <- og_info |>
  transmute(
    OG_name = as.character(GeneName_base),
    Module = if_else(
      as.character(Module) == "JAs",
      "JAs & ETs",
      as.character(Module)
    ),
    Type = as.character(Type)
  ) |>
  filter(!is.na(OG_name), OG_name != "", !is.na(Module), Module != "") |>
  distinct()

if (anyDuplicated(og_annotation[c("OG_name", "Module")])) {
  stop("OG_info contains duplicate OG annotations within a module.")
}

sc_supp7 <- s7 |>
  left_join(og_annotation, by = "OG_name") |>
  mutate(
    Type_group = case_when(
      Type %in% c("B", "BB") ~ "single_B",
      Type %in% c("S", "SS") ~ "single_S",
      TRUE ~ "BS"
    ),
    expr_total = pmax(E_sum_B, 0) + pmax(E_sum_S, 0),
    expr_total = if_else(expr_total > 0, expr_total, NA_real_),
    color_key = paste0(Type_group, "_", as.character(tp))
  ) |>
  filter(!is.na(Module), !is.na(expr_total))

baseline_table <- sc_supp7 |>
  group_by(Module, OG_name) |>
  summarise(
    baseline_raw = exp(mean(log(expr_total + baseline_epsilon), na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(baseline_comp = asinh(baseline_raw))

baseline_parameters <- baseline_table |>
  group_by(Module) |>
  summarise(
    qlo = quantile(baseline_comp, winsor_lo, na.rm = TRUE, names = FALSE),
    qhi = quantile(baseline_comp, winsor_hi, na.rm = TRUE, names = FALSE),
    r_min = min(baseline_comp, na.rm = TRUE),
    r_max = max(baseline_comp, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    from_lo = if_else(is.finite(qlo), qlo, r_min),
    from_hi = if_else(is.finite(qhi), qhi, r_max),
    from_lo = if_else(from_hi > from_lo, from_lo, r_min),
    from_hi = if_else(from_hi > from_lo, from_hi, r_max)
  ) |>
  dplyr::select(Module, from_lo, from_hi)

baseline_table <- baseline_table |>
  left_join(baseline_parameters, by = "Module") |>
  mutate(
    denominator = pmax(from_hi - from_lo, .Machine$double.eps),
    scaled = pmin(pmax((baseline_comp - from_lo) / denominator, 0), 1),
    R_base = size_min + scaled * (size_max - size_min),
    R_base = if_else(
      !is.finite(denominator) | denominator <= .Machine$double.eps,
      (size_min + size_max) / 2,
      R_base
    )
  ) |>
  dplyr::select(Module, OG_name, baseline_raw, R_base)

sc_supp7 <- sc_supp7 |>
  left_join(baseline_table, by = c("Module", "OG_name")) |>
  mutate(
    ratio = (expr_total + baseline_epsilon) / (baseline_raw + baseline_epsilon),
    clipped_change = pmin(pmax(log2(ratio), -clip_change), clip_change),
    size_final = R_base * 2^(change_weight * clipped_change),
    size_final = pmin(pmax(size_final, size_min), size_max)
  )

adir_long <- s8 |>
  left_join(og_annotation, by = "OG_name") |>
  pivot_longer(
    c(Adir_B, Adir_S),
    names_to = "Subgenome",
    values_to = "Adir"
  ) |>
  mutate(Subgenome = recode(Subgenome, Adir_B = "B", Adir_S = "S")) |>
  left_join(
    sc_supp7 |> distinct(Module, OG_name, Type_group),
    by = c("Module", "OG_name")
  ) |>
  filter(!is.na(Module), !is.na(Adir))

order_ogs_for_module <- function(data) {
  data |>
    group_by(OG_name) |>
    summarise(
      Type_major = dplyr::first(Type_group),
      mean_sc = mean(SC_index, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      group_rank = match(Type_major, c("single_S", "BS", "single_B")),
      within_value = case_when(
        Type_major == "single_B" ~ -mean_sc,
        TRUE ~ mean_sc
      )
    ) |>
    arrange(group_rank, within_value) |>
    pull(OG_name)
}

make_module_plot <- function(module) {
  module_sc <- sc_supp7 |> filter(Module == module)
  module_adir <- adir_long |> filter(Module == module)
  og_levels <- order_ogs_for_module(module_sc)
  left_data <- module_adir |>
    mutate(OG = factor(OG_name, levels = og_levels))

  middle_data <- tibble(
    OG_name = og_levels,
    OG = factor(og_levels, levels = og_levels)
  )

  right_data <- module_sc |>
    mutate(OG = factor(OG_name, levels = og_levels)) |>
    arrange(OG, desc(size_final))

  plot_height <- length(og_levels) * row_height + module_fixed_height

  left_plot <- ggplot(
    left_data,
    aes(x = OG, y = Adir, fill = Subgenome)
  ) +
    geom_col(width = 0.8, position = "stack") +
    scale_fill_manual(values = fill_colors, name = "Subgenome") +
    coord_flip() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = "Adir", title = module) +
    theme_minimal(base_size = plot_base_size, base_family = plot_base_family) +
    theme(
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      plot.margin = margin(5.5, 2, 5.5, 4)
    )

  middle_plot <- ggplot(
    middle_data,
    aes(x = 0, y = OG, label = OG_name)
  ) +
    geom_text(
      size = og_label_size, family = plot_base_family,
      hjust = 0.5, check_overlap = TRUE, na.rm = TRUE
    ) +
    coord_cartesian(xlim = c(-1, 1), clip = "off") +
    scale_y_discrete(drop = FALSE) +
    theme_void(base_size = plot_base_size, base_family = plot_base_family) +
    theme(plot.margin = margin(5.5, 0, 5.5, 0))

  right_plot <- ggplot(right_data, aes(x = SC_index, y = OG)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_point(
      aes(size = size_final, colour = color_key),
      alpha = 0.92,
      stroke = 0
    ) +
    scale_size_identity() +
    scale_colour_manual(values = color_map, guide = "none") +
    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.5)
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "cHEB", y = NULL) +
    theme_minimal(base_size = plot_base_size, base_family = plot_base_family) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(5.5, 4, 5.5, 2)
    )

  combined_plot <- left_plot + middle_plot + right_plot +
    plot_layout(widths = c(left_width, middle_width, right_width)) +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

  attr(combined_plot, "supp7_height") <- plot_height
  combined_plot
}

modules_supp7 <- module_order[
  module_order %in% unique(adir_long$Module) &
    module_order %in% unique(sc_supp7$Module)
]

paired_rows <- list(
  c("PRRs", "Ca2+"),
  c("RLCKs", "JAs & ETs"),
  c("NLRs", "SAs")
)
paired_rows <- lapply(paired_rows, function(x) x[x %in% modules_supp7])

# Keep nested patchwork panels intact within the outer layout.
as_layout_panel <- function(plot) {
  patchwork::wrap_elements(full = plot)
}

for (module in modules_supp7) {
  module_plot <- make_module_plot(module)
  module_height <- attr(module_plot, "supp7_height")
  safe_module_name <- gsub("[\\/:*?\"<>|]", "_", module)
  ggsave(
    file.path(individual_output_dir, paste0(safe_module_name, "_combo.pdf")),
    module_plot,
    width = left_width + middle_width + right_width,
    height = module_height,
    units = "in",
    dpi = 300
  )

}

# Preserve each module's physical height when assembling the two columns.
left_modules <- c("PRRs", "RLCKs", "NLRs", "MAPKs", "WRKYs")
right_modules <- c("Ca2+", "JAs & ETs", "SAs", "Cucurbitacins")
composite_modules <- c(left_modules, right_modules)
composite_modules <- composite_modules[composite_modules %in% modules_supp7]
composite_plots <- setNames(
  lapply(composite_modules, make_module_plot), composite_modules
)
composite_heights <- vapply(
  composite_plots, function(p) attr(p, "supp7_height"), numeric(1)
)
left_total_height <- sum(composite_heights[left_modules])

left_bottom <- left_total_height
composite_layers <- list()
for (module in left_modules[left_modules %in% composite_modules]) {
  plot_height <- composite_heights[[module]]
  left_bottom <- left_bottom - plot_height
  composite_layers[[length(composite_layers) + 1L]] <- cowplot::draw_plot(
    composite_plots[[module]],
    x = 0, y = left_bottom / left_total_height,
    width = 0.5, height = plot_height / left_total_height
  )
}

ca_height <- composite_heights[["Ca2+"]]
ja_height <- composite_heights[["JAs & ETs"]]
sa_height <- composite_heights[["SAs"]]
cucurbitacin_height <- composite_heights[["Cucurbitacins"]]
ca_to_ja_gap <- max(0, composite_heights[["PRRs"]] - ca_height)
# Overlap blank outer margins only; do not overlap data regions.
sa_cuci_margin_overlap <- 0.15
ja_to_sa_gap <- left_total_height - ca_height - ca_to_ja_gap - ja_height -
  sa_height - cucurbitacin_height + sa_cuci_margin_overlap
if (ja_to_sa_gap < 0) {
  stop("The right column exceeds the left-column height; bottom alignment would overlap data regions.")
}

right_top <- left_total_height
right_top <- right_top - ca_height
composite_layers[[length(composite_layers) + 1L]] <- cowplot::draw_plot(
  composite_plots[["Ca2+"]],
  x = 0.5, y = right_top / left_total_height,
  width = 0.5, height = ca_height / left_total_height
)
right_top <- right_top - ca_to_ja_gap - ja_height
composite_layers[[length(composite_layers) + 1L]] <- cowplot::draw_plot(
  composite_plots[["JAs & ETs"]],
  x = 0.5, y = right_top / left_total_height,
  width = 0.5, height = ja_height / left_total_height
)
right_top <- right_top - ja_to_sa_gap - sa_height
composite_layers[[length(composite_layers) + 1L]] <- cowplot::draw_plot(
  composite_plots[["SAs"]],
  x = 0.5, y = right_top / left_total_height,
  width = 0.5, height = sa_height / left_total_height
)
right_top <- right_top + sa_cuci_margin_overlap - cucurbitacin_height
composite_layers[[length(composite_layers) + 1L]] <- cowplot::draw_plot(
  composite_plots[["Cucurbitacins"]],
  x = 0.5, y = right_top / left_total_height,
  width = 0.5, height = cucurbitacin_height / left_total_height
)

supp_fig7 <- cowplot::ggdraw()
for (layer in composite_layers) supp_fig7 <- supp_fig7 + layer

pdf_path <- file.path(output_dir, "Supp_Fig_7.pdf")
png_path <- file.path(output_dir, "Supp_Fig_7.png")
combined_height <- left_total_height

ggsave(
  pdf_path,
  supp_fig7,
  width = 11,
  height = combined_height,
  units = "in",
  dpi = 300
)

ggsave(
  png_path,
  supp_fig7,
  width = 11,
  height = combined_height,
  units = "in",
  dpi = 300
)

message("Supplementary Fig. 7 plotting completed.")
message("Output directory: ", output_dir)
}))
