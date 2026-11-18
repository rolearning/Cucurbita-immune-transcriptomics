#!/usr/bin/env Rscript

# Analyze final model 0.5_22_2_50 from the Stage 1 workspace and manual-label table.
# Write module/trait statistics, hub tables, DR enrichment, role-space plots and formal_results.rds.
# Use hybrid bicor for trait associations and bicor for module membership and connectivity.
# Select modules at nominal P < 0.05 and report FDR; hubs require MM >= 0.80 and
# within-module kWithin percentile >= 0.70, without a GS cutoff.

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

cfg2 <- list(
  stage1_workspace = file.path(RELEASE_ROOT, "6.WGCNA/stage1/stage1_workspace.rds"),
  candidate_file = file.path(RELEASE_ROOT, "6.WGCNA/stage1/candidate_models.csv"),
  output_dir = file.path(RELEASE_ROOT, "6.WGCNA/stage2_only_0.5_22_2_50"),

  nominal_p_cut = 0.05,
  fdr_cut = 0.05,
  trait_correlation = "hybrid_bicor",

  hub_mm_cut = 0.80,
  hub_kwithin_rank_cut = 0.70,

  # Manual labels match GeneName; Shift_x/Shift_y set the initial ggrepel displacement.
  label_mode = "manual",
  manual_label_file = file.path(RELEASE_ROOT, "6.WGCNA/role_space_manual_labels.tsv"),
  manual_label_gene_column = "GeneName",

  top_dr_hubs_to_label_per_module = 0,
  manual_label_ogs = character(),

  # Enrichment is disabled here; formal eggNOG ORA is a shared downstream step.
  run_functional_enrichment = FALSE
)

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

get_colors <- function(net, genes) {
  colors <- net$colors
  if (is.null(colors)) stop("Network object has no colors.")
  if (is.null(names(colors))) names(colors) <- genes
  colors[genes]
}

get_ordered_mes <- function(net, dat_expr, colors) {
  mes <- net$MEs
  if (is.null(mes)) mes <- WGCNA::moduleEigengenes(dat_expr, colors)$eigengenes
  WGCNA::orderMEs(mes)
}

cor_and_p <- function(x, y, method, stage1_cfg) {
  if (method == "pearson") {
    r <- WGCNA::cor(x, y, use = "pairwise.complete.obs")
    p <- WGCNA::corPvalueStudent(r, nSamples = nrow(x))
    return(list(cor = r, p = p, method = "Pearson"))
  }
  if (method == "hybrid_bicor") {
    z <- WGCNA::bicorAndPvalue(
      x, y, use = "pairwise.complete.obs",
      robustX = TRUE, robustY = FALSE,
      maxPOutliers = stage1_cfg$max_p_outliers,
      pearsonFallback = stage1_cfg$pearson_fallback
    )
    return(list(cor = z$bicor, p = z$p, method = "hybrid bicor"))
  }
  stop("trait_correlation must be 'pearson' or 'hybrid_bicor'.")
}

classify_modules <- function(cor_mat, p_mat, fdr_mat, p_cut) {
  priority <- c(
    mid_24_vs_mean_0_48 = "mid_transient",
    late_48_vs_rest = "late_response",
    early_0_vs_rest = "early_response",
    t_trend = "trend"
  )

  purrr::map_dfr(rownames(cor_mat), function(me) {
    selected_trait <- NA_character_
    selected_class <- "none"
    for (tr in names(priority)) {
      if (tr %in% colnames(p_mat) && is.finite(p_mat[me, tr]) &&
          p_mat[me, tr] < p_cut) {
        selected_trait <- tr
        selected_class <- unname(priority[tr])
        break
      }
    }
    if (is.na(selected_trait)) {
      return(tibble(
        module = me, class = selected_class, direction = NA_character_,
        top_trait = NA_character_, correlation = NA_real_, p_value = NA_real_,
        FDR = NA_real_
      ))
    }
    r <- cor_mat[me, selected_trait]
    tibble(
      module = me,
      class = selected_class,
      direction = ifelse(r > 0, "up", "down"),
      top_trait = selected_trait,
      correlation = unname(r),
      p_value = unname(p_mat[me, selected_trait]),
      FDR = unname(fdr_mat[me, selected_trait])
    )
  })
}

module_connectivity_same_network <- function(dat_expr, colors, power, stage1_cfg) {
  genes <- colnames(dat_expr)
  modules <- setdiff(unique(as.character(colors)), "grey")
  purrr::map_dfr(modules, function(mod) {
    mod_genes <- genes[colors == mod]
    n <- length(mod_genes)
    if (n == 1L) {
      return(tibble(Gene = mod_genes, Module = mod,
                    kWithin = 0, kWithin_rank = 0, kWithin_scaled = 0))
    }
    adj <- WGCNA::adjacency(
      dat_expr[, mod_genes, drop = FALSE],
      type = stage1_cfg$network_type,
      power = power,
      corFnc = "bicor",
      corOptions = list(
        use = "pairwise.complete.obs",
        maxPOutliers = stage1_cfg$max_p_outliers,
        pearsonFallback = stage1_cfg$pearson_fallback
      )
    )
    k <- rowSums(adj) - diag(adj)
    tibble(
      Gene = mod_genes,
      Module = mod,
      kWithin = as.numeric(k),
      kWithin_rank = dplyr::percent_rank(k),
      kWithin_scaled = if (max(k) > 0) k / max(k) else 0
    )
  })
}

trait_display_name <- function(x) {
  dplyr::recode(x,
    early_0_vs_rest = "0 hpi",
    mid_24_vs_mean_0_48 = "24 hpi",
    late_48_vs_rest = "48 hpi",
    t_trend = "Trend",
    .default = x
  )
}

module_order_from_summary <- function(summary_df) {
  class_levels <- c("early_response", "mid_transient", "late_response", "trend")
  summary_df %>%
    filter(!is.na(top_trait)) %>%
    mutate(class = factor(class, levels = class_levels)) %>%
    arrange(class, desc(direction == "up"), desc(abs(correlation))) %>%
    pull(module) %>% unique()
}

main_fig7a_module_order <- c(
  "MEgrey60", "MEpurple", "MEturquoise", "MEred",
  "MEgreen", "MEpink", "MEmagenta", "MEbrown"
)

make_heatmap_glyph <- function(cor_mat, p_mat, mes, sample_info,
                               module_summary, cor_method) {
  missing_modules <- setdiff(
    main_fig7a_module_order,
    intersect(rownames(cor_mat), colnames(mes))
  )
  if (length(missing_modules) > 0) {
    stop(
      "Main Fig. 7A is missing fixed-order module(s): ",
      paste(missing_modules, collapse = ", ")
    )
  }
  modules <- main_fig7a_module_order

  hm <- as.data.frame(cor_mat[modules, , drop = FALSE]) %>%
    tibble::rownames_to_column("module") %>%
    pivot_longer(-module, names_to = "trait", values_to = "correlation") %>%
    left_join(
      as.data.frame(p_mat[modules, , drop = FALSE]) %>%
        tibble::rownames_to_column("module") %>%
        pivot_longer(-module, names_to = "trait", values_to = "p_value"),
      by = c("module", "trait")
    ) %>%
    mutate(
      module = factor(module, levels = rev(modules)),
      trait_label = factor(trait_display_name(trait),
                           levels = c("0 hpi", "24 hpi", "48 hpi", "Trend")),
      star = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      label = sprintf("%.2f%s", correlation, star)
    )

  p_heat <- ggplot(hm, aes(trait_label, module, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = cor_method) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1,
                                 color = "#4D4D4D"),
      axis.text.y = element_text(size = 12, color = "#4D4D4D"),
      legend.position = "left",
      plot.margin = margin(5.5, 0, 5.5, 5.5)
    )

  me_long <- mes[, modules, drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "module", values_to = "ME") %>%
    left_join(as.data.frame(sample_info), by = "sample") %>%
    group_by(module, time_hpi) %>%
    summarise(median_ME = median(ME), .groups = "drop") %>%
    mutate(module = factor(module, levels = modules),
           time_hpi = factor(time_hpi, levels = sort(unique(time_hpi))))

  p_glyph <- ggplot(me_long, aes(time_hpi, median_ME, group = module)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.6) +
    scale_y_continuous(position = "right",
                       expand = expansion(mult = c(0.10, 0.10))) +
    facet_grid(module ~ ., scales = "fixed") +
    labs(x = NULL, y = "ME (median)") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.spacing.y = grid::unit(0, "pt"),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y.right = element_text(),
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    )

  p_heat + p_glyph + patchwork::plot_layout(widths = c(4, 1.35))
}

dr_enrichment <- function(colors, dr_genes) {
  genes <- names(colors)
  modules <- setdiff(unique(as.character(colors)), "grey")
  purrr::map_dfr(modules, function(mod) {
    in_mod <- genes[colors == mod]
    a <- sum(in_mod %in% dr_genes)
    b <- length(in_mod) - a
    c <- sum(genes %in% dr_genes) - a
    d <- length(genes) - a - b - c
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                      alternative = "greater")
    tibble(module = paste0("ME", mod), module_size = length(in_mod),
           DR_in_module = a, nonDR_in_module = b,
           DR_outside_module = c, nonDR_outside_module = d,
           odds_ratio = unname(ft$estimate),
           log2_odds_ratio_plot = log2(((a + 0.5) * (d + 0.5)) /
                                         ((b + 0.5) * (c + 0.5))),
           p_value = ft$p.value)
  }) %>% mutate(FDR = p.adjust(p_value, method = "BH"))
}

add_annotation <- function(gene_tbl, oginfo, dr_gene_column) {
  if (is.null(oginfo)) return(gene_tbl)
  ann <- oginfo
  colnames(ann)[match(dr_gene_column, colnames(ann))] <- "Gene"
  # Drop annotation Module before joining to avoid .x/.y name collisions.
  if ("Module" %in% colnames(ann) && "Module" %in% colnames(gene_tbl)) {
    colnames(ann)[colnames(ann) == "Module"] <- "Functional_module"
  }
  gene_tbl %>% left_join(ann, by = "Gene")
}

format_subgenome_label <- function(x) {
  x <- as.character(x)
  ifelse(x %in% c("B", "S"), paste0("Cm", x), x)
}

read_manual_role_labels <- function(file) {
  if (is.null(file) || !nzchar(file) || !file.exists(file)) {
    stop("Manual role-space label file does not exist: ", file)
  }
  x <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("GeneName", "Shift_x", "Shift_y")
  missing <- setdiff(required, colnames(x))
  if (length(missing)) {
    stop("Manual label file is missing columns: ", paste(missing, collapse = ", "))
  }
  if (anyNA(x$GeneName) || any(!nzchar(trimws(x$GeneName)))) {
    stop("Manual label file contains an empty GeneName.")
  }
  if (anyDuplicated(x$GeneName)) {
    stop("Manual label file contains duplicated GeneName values: ",
         paste(unique(x$GeneName[duplicated(x$GeneName)]), collapse = ", "))
  }
  x$Shift_x <- suppressWarnings(as.numeric(x$Shift_x))
  x$Shift_y <- suppressWarnings(as.numeric(x$Shift_y))
  if (anyNA(x$Shift_x) || anyNA(x$Shift_y)) {
    stop("Shift_x and Shift_y must be numeric and cannot be empty.")
  }
  x %>% dplyr::select(GeneName, Shift_x, Shift_y)
}

read_cm_eggnog_annotation <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    stop("Cm eggNOG annotation file does not exist: ", path)
  }

  # Find the #query header after any eggNOG metadata lines.
  header_line <- readLines(path, n = 20, warn = FALSE)
  header_idx <- which(grepl("^#query\\t", header_line))[1]
  if (is.na(header_idx)) {
    stop("Cannot find the eggNOG annotation header beginning with #query: ", path)
  }
  col_names <- strsplit(sub("^#", "", header_line[header_idx]), "\\t")[[1]]

  ann <- read.delim(
    path,
    sep = "\t",
    header = FALSE,
    skip = header_idx,
    comment.char = "#",
    quote = "",
    fill = TRUE,
    col.names = col_names,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ann <- as_tibble(ann)
  required <- c("query", "GOs", "KEGG_Pathway", "KEGG_ko")
  missing <- setdiff(required, colnames(ann))
  if (length(missing)) {
    stop("Cm eggNOG annotation is missing column(s): ",
         paste(missing, collapse = ", "))
  }

  ann %>%
    transmute(
      Gene = sub("\\..*$", "", as.character(query)),
      GOs = as.character(GOs),
      KEGG_Pathway = as.character(KEGG_Pathway),
      KEGG_ko = as.character(KEGG_ko)
    ) %>%
    filter(!is.na(Gene), Gene != "") %>%
    distinct()
}

run_ora_from_terms <- function(foreground, universe, term2gene,
                               term2name = NULL, module, source,
                               cfg2) {
  foreground <- intersect(unique(as.character(foreground)),
                          unique(as.character(universe)))
  universe <- unique(as.character(universe))
  term2gene <- term2gene %>%
    transmute(term = as.character(term), gene = as.character(gene)) %>%
    filter(!is.na(term), term != "", !is.na(gene), gene != "") %>%
    distinct()
  if (!length(foreground) || !length(universe) || !nrow(term2gene)) {
    return(tibble())
  }

  annotated_universe <- intersect(universe, unique(term2gene$gene))
  annotated_foreground <- intersect(foreground, annotated_universe)
  if (!length(annotated_foreground) || !length(annotated_universe)) {
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
    minGSSize = cfg2$enrichment_min_gs_size,
    maxGSSize = cfg2$enrichment_max_gs_size
  )
  result_df <- as.data.frame(result)
  if (!nrow(result_df)) return(tibble())
  as_tibble(result_df) %>%
    mutate(
      module = module,
      source = source,
      n_foreground_annotated = length(annotated_foreground),
      n_universe_annotated = length(annotated_universe),
      .before = 1
    )
}

run_functional_enrichment <- function(run_dir, significant_modules,
                                      network_genes, colors, stage1, cfg2) {
  if (!cfg2$run_functional_enrichment) return(invisible(NULL))

  # Restrict annotation and enrichment to genes in the actual network.
  ann <- read_cm_eggnog_annotation(cfg2$eggnog_annotation_file)
  network_genes <- unique(as.character(network_genes))
  annotated_network_genes <- intersect(network_genes, ann$Gene)
  enrich_dir <- file.path(run_dir, "functional_enrichment")
  dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

  coverage <- tibble(
    universe_network_gene_n = length(network_genes),
    universe_with_eggNOG_annotation_n = length(annotated_network_genes),
    universe_annotation_coverage = ifelse(
      length(network_genes),
      length(annotated_network_genes) / length(network_genes),
      NA_real_
    )
  )
  write.csv(coverage, file.path(enrich_dir, "mapping_coverage.csv"),
            row.names = FALSE)

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

  kegg_terms <- ann %>%
    dplyr::select(Gene, KEGG_Pathway) %>%
    filter(!is.na(KEGG_Pathway), KEGG_Pathway != "", KEGG_Pathway != "-") %>%
    separate_rows(KEGG_Pathway, sep = ",") %>%
    transmute(term = trimws(KEGG_Pathway), gene = Gene) %>%
    filter(term != "", term != "-") %>%
    distinct()
  kegg_names <- kegg_terms %>%
    distinct(term) %>%
    mutate(name = term)

  # Keep KO terms separate from KEGG pathways.
  ko_terms <- ann %>%
    dplyr::select(Gene, KEGG_ko) %>%
    filter(!is.na(KEGG_ko), KEGG_ko != "", KEGG_ko != "-") %>%
    separate_rows(KEGG_ko, sep = ",") %>%
    transmute(term = trimws(KEGG_ko), gene = Gene) %>%
    filter(grepl("^ko:[Kk][0-9]{5}$|^[Kk][0-9]{5}$", term)) %>%
    mutate(term = ifelse(grepl("^ko:", term), term, paste0("ko:", term))) %>%
    distinct()
  ko_names <- ko_terms %>%
    distinct(term) %>%
    mutate(name = term)

  module_genes <- setNames(
    lapply(significant_modules, function(me) {
      mod <- sub("^ME", "", me)
      network_genes[as.character(colors) == mod]
    }),
    significant_modules
  )

  enrichment_rows <- list()
  coverage_rows <- list()
  row_index <- 1L
  coverage_index <- 1L
  for (me in names(module_genes)) {
    fg <- module_genes[[me]]
    coverage_rows[[coverage_index]] <- tibble(
      module = me,
      n_module_gene = length(fg),
      n_module_gene_with_eggNOG = length(intersect(fg, annotated_network_genes)),
      module_annotation_coverage = ifelse(
        length(fg), length(intersect(fg, annotated_network_genes)) / length(fg), NA_real_
      ),
      n_universe_gene = length(network_genes),
      n_universe_gene_with_eggNOG = length(annotated_network_genes)
    )
    coverage_index <- coverage_index + 1L

    if (cfg2$run_go) {
      for (ont in c("BP", "MF", "CC")) {
        terms <- go_terms %>%
          filter(ONTOLOGY == ont) %>%
          transmute(term = GO, gene = Gene, name = TERM)
        term2name <- terms %>% dplyr::select(term, name) %>% distinct()
        out <- run_ora_from_terms(
          fg, network_genes, terms %>% dplyr::select(term, gene), term2name,
          me, paste0("GO_", ont), cfg2
        )
        if (nrow(out)) {
          out$ONTOLOGY <- ont
          enrichment_rows[[row_index]] <- out
          row_index <- row_index + 1L
        }
      }
    }

    if (cfg2$run_kegg) {
      out <- run_ora_from_terms(
        fg, network_genes, kegg_terms, kegg_names,
        me, "KEGG_Pathway", cfg2
      )
      if (nrow(out)) {
        enrichment_rows[[row_index]] <- out
          row_index <- row_index + 1L
        }

      out_ko <- run_ora_from_terms(
        fg, network_genes, ko_terms, ko_names,
        me, "KEGG_KO", cfg2
      )
      if (nrow(out_ko)) {
        enrichment_rows[[row_index]] <- out_ko
        row_index <- row_index + 1L
      }
      }
  }

  enrichment <- bind_rows(enrichment_rows)
  coverage_by_module <- bind_rows(coverage_rows)
  write.csv(coverage_by_module,
            file.path(enrich_dir, "module_annotation_coverage.csv"),
            row.names = FALSE)
  write.csv(enrichment,
            file.path(enrich_dir, "GO_KEGG_enrichment_all_terms.csv"),
            row.names = FALSE)
  write.csv(
    enrichment %>% filter(source == "KEGG_Pathway"),
    file.path(enrich_dir, "KEGG_Pathway_enrichment_all_terms.csv"),
    row.names = FALSE
  )
  write.csv(
    enrichment %>% filter(source == "KEGG_KO"),
    file.path(enrich_dir, "KEGG_KO_enrichment_all_terms.csv"),
    row.names = FALSE
  )

  if (nrow(enrichment)) {
    write.csv(
      enrichment %>% filter(!is.na(p.adjust), p.adjust <= cfg2$enrichment_p_cut),
      file.path(enrich_dir, "GO_KEGG_enrichment_padj_significant.csv"),
      row.names = FALSE
    )
    write.csv(
      enrichment %>% filter(source == "GO_BP",
                            !is.na(p.adjust), p.adjust <= cfg2$enrichment_p_cut),
      file.path(enrich_dir, "GO_BP_enrichment_padj_significant.csv"),
      row.names = FALSE
    )
    write.csv(
      enrichment %>% filter(source == "KEGG_Pathway",
                            !is.na(p.adjust), p.adjust <= cfg2$enrichment_p_cut),
      file.path(enrich_dir, "KEGG_enrichment_padj_significant.csv"),
      row.names = FALSE
    )
    write.csv(
      enrichment %>% filter(source == "KEGG_KO",
                            !is.na(p.adjust), p.adjust <= cfg2$enrichment_p_cut),
      file.path(enrich_dir, "KEGG_KO_enrichment_padj_significant.csv"),
      row.names = FALSE
    )
  }
  invisible(enrichment)
}

stage1 <- readRDS(cfg2$stage1_workspace)
candidates <- read.csv(cfg2$candidate_file, check.names = FALSE)
if (!"use_in_stage2" %in% colnames(candidates)) {
  stop("candidate_file must contain use_in_stage2.")
}
use <- tolower(trimws(as.character(candidates$use_in_stage2))) %in%
  c("true", "t", "1", "yes", "y")
selected_ids <- "0.5_22_2_50"
if (!length(selected_ids)) {
  stop("No models selected. Set use_in_stage2=TRUE in candidate_models.csv.")
}
missing_ids <- setdiff(selected_ids, names(stage1$runs))
if (length(missing_ids)) stop("Unknown run_id(s): ", paste(missing_ids, collapse = ", "))

dir.create(cfg2$output_dir, recursive = TRUE, showWarnings = FALSE)

candidate_comparison <- list()

for (run_id in selected_ids) {
  message("Writing formal outputs for ", run_id)
  run <- stage1$runs[[run_id]]
  pre <- stage1$preprocess[[run$keep_key]]
  dat_expr <- pre$datExpr
  traits <- pre$traits
  sample_info <- pre$sample_info
  net <- run$net
  colors <- get_colors(net, colnames(dat_expr))
  mes <- get_ordered_mes(net, dat_expr, colors)

  run_dir <- file.path(cfg2$output_dir, safe_name(run_id))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  mt <- cor_and_p(mes, traits, cfg2$trait_correlation, stage1$config)
  mt_fdr <- matrix(
    p.adjust(as.vector(mt$p), method = "BH"),
    nrow = nrow(mt$p), dimnames = dimnames(mt$p)
  )
  module_summary <- classify_modules(mt$cor, mt$p, mt_fdr, cfg2$nominal_p_cut)
  size_tab <- table(colors)
  module_summary <- module_summary %>%
    mutate(module_color = sub("^ME", "", module),
           module_size = as.numeric(size_tab[module_color])) %>%
    dplyr::select(-module_color) %>%
    arrange(factor(class, levels = c("early_response", "mid_transient",
                                     "late_response", "trend", "none")),
            desc(abs(correlation)))
  significant_modules <- module_summary$module[!is.na(module_summary$top_trait)]

  write.csv(mt$cor, file.path(run_dir, "module_trait_correlation.csv"))
  write.csv(mt$p, file.path(run_dir, "module_trait_raw_P.csv"))
  write.csv(mt_fdr, file.path(run_dir, "module_trait_FDR.csv"))
  write.csv(module_summary, file.path(run_dir, "module_summary.csv"), row.names = FALSE)

  combined <- make_heatmap_glyph(mt$cor, mt$p, mes, sample_info,
                                 module_summary, mt$method)
  if (!is.null(combined)) {
    ggsave(file.path(run_dir, "module_trait_heatmap_and_ME_glyph.pdf"),
           combined, width = 6,
           height = max(5, 0.55 * length(significant_modules) + 1.5))
  }

  me_long <- mes %>% as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "ME", values_to = "value") %>%
    left_join(as.data.frame(sample_info), by = "sample") %>%
    filter(ME %in% significant_modules)
  if (nrow(me_long)) {
    p_me <- ggplot(me_long, aes(factor(time_hpi), value, color = factor(time_hpi))) +
      geom_boxplot(outlier.shape = NA, width = 0.55) +
      geom_jitter(width = 0.08, size = 1.8) +
      facet_wrap(~ME, scales = "free_y", ncol = 4, nrow = 2) +
      labs(x = "Time (hpi)", y = "Module eigengene", color = "hpi") +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(size = 11))
    ggsave(file.path(run_dir, "module_eigengene_profiles.pdf"), p_me,
           width = 12, height = 7.5)
  }

  # MM and its P values use bicor, matching the network adjacency.
  mm_result <- WGCNA::bicorAndPvalue(
    dat_expr,
    mes,
    use = "pairwise.complete.obs",
    robustX = TRUE,
    robustY = TRUE,
    maxPOutliers = stage1$config$max_p_outliers,
    pearsonFallback = stage1$config$pearson_fallback
  )
  mm_all <- mm_result$bicor
  mm_p_all <- mm_result$p

  # Hybrid bicor robustifies expression but leaves experimental traits non-robust.
  gs_result <- cor_and_p(
    dat_expr,
    traits,
    cfg2$trait_correlation,
    stage1$config
  )
  gs_all <- gs_result$cor
  gs_p_all <- gs_result$p
  connectivity <- module_connectivity_same_network(
    dat_expr, colors, run$power, stage1$config
  )

  best_trait <- module_summary %>%
    filter(!is.na(top_trait)) %>%
    transmute(ME = module, Module = sub("^ME", "", module), BestTrait = top_trait,
              BestTrait_corME = correlation, BestTrait_pME = p_value,
              BestTrait_FDR_ME = FDR, class, direction)

  gene_tbl <- tibble(Gene = colnames(dat_expr), Module = as.character(colors)) %>%
    filter(Module != "grey") %>%
    mutate(ME = paste0("ME", Module)) %>%
    filter(ME %in% significant_modules) %>%
    left_join(connectivity, by = c("Gene", "Module")) %>%
    left_join(best_trait, by = c("ME", "Module")) %>%
    rowwise() %>%
    mutate(
      MM = mm_all[Gene, ME],
      MM_p = mm_p_all[Gene, ME],
      GS = gs_all[Gene, BestTrait],
      GS_p = gs_p_all[Gene, BestTrait]
    ) %>%
    ungroup() %>%
    mutate(
      isDR = Gene %in% stage1$dr_genes,
      isHUB = MM >= cfg2$hub_mm_cut &
        kWithin_rank >= cfg2$hub_kwithin_rank_cut,
      BestTrait_label = trait_display_name(BestTrait)
    )
  gene_tbl <- add_annotation(gene_tbl, stage1$oginfo,
                             stage1$config$dr_gene_column)

  # Divergent DR orthogroups span multiple modules across all nongrey network genes.
  all_module_membership <- tibble(
    Gene = colnames(dat_expr),
    Module = as.character(colors),
    ME = paste0("ME", as.character(colors)),
    isDR = colnames(dat_expr) %in% stage1$dr_genes
  ) %>%
    filter(Module != "grey")
  all_module_membership <- add_annotation(
    all_module_membership, stage1$oginfo, stage1$config$dr_gene_column
  )

  og_col <- intersect(c("GeneName_base", "OG"),
                      colnames(all_module_membership))[1]
  if (length(og_col) && !is.na(og_col)) {
    div_ogs <- all_module_membership %>%
      filter(isDR, !is.na(.data[[og_col]])) %>%
      group_by(.data[[og_col]]) %>%
      summarise(n_ME_modules = n_distinct(ME), .groups = "drop") %>%
      filter(n_ME_modules > 1) %>%
      pull(1)
    gene_tbl$isDiv <- gene_tbl$isDR & gene_tbl[[og_col]] %in% div_ogs
  } else {
    gene_tbl$isDiv <- FALSE
  }

  subgenome_col <- intersect(c("subgenome", "Subgenome"),
                             colnames(gene_tbl))[1]
  if (length(og_col) && !is.na(og_col) &&
      length(subgenome_col) && !is.na(subgenome_col)) {
    gene_tbl$plot_label <- ifelse(
      !is.na(gene_tbl[[og_col]]) & !is.na(gene_tbl[[subgenome_col]]),
      paste0(gene_tbl[[og_col]], "_",
             format_subgenome_label(gene_tbl[[subgenome_col]])),
      gene_tbl$Gene
    )
  } else {
    gene_tbl$plot_label <- gene_tbl$Gene
  }

  write.csv(gene_tbl, file.path(run_dir, "WGCNA_gene_metrics_ALL.csv"), row.names = FALSE)
  write.csv(filter(gene_tbl, isDR),
            file.path(run_dir, "WGCNA_gene_metrics_DR.csv"), row.names = FALSE)
  write.csv(filter(gene_tbl, isHUB),
            file.path(run_dir, "WGCNA_hubs.csv"), row.names = FALSE)
  write.csv(filter(gene_tbl, isHUB & isDR),
            file.path(run_dir, "WGCNA_DR_hubs.csv"), row.names = FALSE)

  dr_tbl <- dr_enrichment(colors, stage1$dr_genes) %>%
    left_join(module_summary %>% dplyr::select(module, class, direction, top_trait),
              by = c("module" = "module")) %>%
    arrange(FDR, p_value)
  write.csv(dr_tbl, file.path(run_dir, "DR_module_enrichment.csv"), row.names = FALSE)
  if (nrow(dr_tbl)) {
    preferred_dr_order <- module_order_from_summary(module_summary)
    preferred_dr_order <- c(preferred_dr_order,
                            setdiff(dr_tbl$module, preferred_dr_order))
    dr_plot <- dr_tbl %>%
      mutate(neg_log10_fdr = -log10(pmax(FDR, .Machine$double.xmin)),
             module = factor(module, levels = rev(preferred_dr_order))) %>%
      ggplot(aes(log2_odds_ratio_plot, module)) +
      geom_point(aes(size = DR_in_module, color = neg_log10_fdr)) +
      scale_color_gradient(low = "#abdcff", high = "darkblue") +
      labs(x = "Enrichment strength (continuity-corrected log2 odds ratio)",
           y = "WGCNA module",
           color = expression(-log[10](FDR)), size = "DR genes in module") +
      theme_bw(base_size = 10)
    ggsave(file.path(run_dir, "DR_module_enrichment.pdf"), dr_plot,
           width = 5.5, height = max(4, 0.35 * nrow(dr_tbl) + 1.5))
  }

  plot_tbl <- gene_tbl %>%
    mutate(
      absMM = abs(MM),
      point_type = case_when(
        isHUB & isDR ~ "Central IC genes",
        isDR ~ "IC genes",
        TRUE ~ "Other genes"
      )
    )

  label_mode <- match.arg(cfg2$label_mode,
                          c("automatic", "manual", "both"))
  plot_tbl <- plot_tbl %>%
    mutate(label = NA_character_, Shift_x = 0, Shift_y = 0)

  if (label_mode %in% c("automatic", "both")) {
    label_gene_ids <- plot_tbl %>%
      filter(isHUB & isDR) %>% group_by(ME) %>%
      arrange(desc(MM), desc(kWithin_rank), .by_group = TRUE) %>%
      slice_head(n = cfg2$top_dr_hubs_to_label_per_module) %>%
      ungroup() %>% pull(Gene)
    if (length(cfg2$manual_label_ogs) && length(og_col) && !is.na(og_col)) {
      label_gene_ids <- union(
        label_gene_ids,
        plot_tbl$Gene[plot_tbl[[og_col]] %in% cfg2$manual_label_ogs]
      )
    }
    plot_tbl$label <- ifelse(
      plot_tbl$Gene %in% label_gene_ids,
      plot_tbl$plot_label,
      NA_character_
    )
  }

  if (label_mode %in% c("manual", "both")) {
    manual_labels <- read_manual_role_labels(cfg2$manual_label_file) %>%
      transmute(.manual_key = as.character(GeneName),
                .manual_label = as.character(GeneName),
                .manual_shift_x = Shift_x,
                .manual_shift_y = Shift_y)
    key_col <- cfg2$manual_label_gene_column
    if (!key_col %in% colnames(plot_tbl)) {
      stop("manual_label_gene_column is absent from gene annotations: ", key_col)
    }
    plot_tbl$.manual_key <- as.character(plot_tbl[[key_col]])
    matched_keys <- plot_tbl$.manual_key[
      plot_tbl$.manual_key %in% manual_labels$.manual_key
    ]
    if (anyDuplicated(matched_keys)) {
      stop("A manual GeneName matches more than one plotted gene: ",
           paste(unique(matched_keys[duplicated(matched_keys)]), collapse = ", "))
    }
    plot_tbl <- plot_tbl %>%
      left_join(manual_labels, by = ".manual_key") %>%
      mutate(
        label = if_else(!is.na(.manual_label), .manual_label, label),
        Shift_x = coalesce(.manual_shift_x, Shift_x),
        Shift_y = coalesce(.manual_shift_y, Shift_y)
      ) %>%
      dplyr::select(-.manual_key, -.manual_label,
             -.manual_shift_x, -.manual_shift_y)

    manual_label_status <- manual_labels %>%
      transmute(GeneName = .manual_key,
                Shift_x = .manual_shift_x,
                Shift_y = .manual_shift_y,
                found_in_candidate = .manual_key %in% matched_keys)
    write.csv(manual_label_status,
              file.path(run_dir, "role_space_manual_label_status.csv"),
              row.names = FALSE)
  }

  role_label_tbl <- filter(plot_tbl, !is.na(label))

  point_colours <- c(
    "Other genes" = "#BDBDBD",
    "IC genes" = "#F2B6CF",
    "Central IC genes" = "#8E63B6"
  )

  # Overlay outlined divergent genes after unoutlined genes, with central divergent genes last.
  # Use the same data, labels and scales for the main and supplementary role-space layouts.
  make_role_space_plot <- function(module_ids, ncol, nrow) {
    plot_tbl_subset <- plot_tbl %>%
      filter(ME %in% module_ids) %>%
      mutate(ME = factor(ME, levels = module_ids))
    role_label_tbl_subset <- role_label_tbl %>%
      filter(ME %in% module_ids) %>%
      mutate(ME = factor(ME, levels = module_ids))

    ggplot(plot_tbl_subset, aes(absMM, GS, size = kWithin)) +
      geom_point(
        data = filter(plot_tbl_subset, point_type == "Other genes"),
        aes(color = point_type), alpha = 0.95
      ) +
      geom_point(
        data = filter(plot_tbl_subset, point_type == "IC genes", !isDiv),
        aes(color = point_type), alpha = 0.95
      ) +
      geom_point(
        data = filter(plot_tbl_subset, point_type == "Central IC genes", !isDiv),
        aes(color = point_type), alpha = 0.9
      ) +
      geom_point(
        data = filter(plot_tbl_subset, point_type == "IC genes", isDiv),
        aes(fill = point_type, shape = "Divergent IC genes"),
        color = "black", alpha = 0.90, stroke = 0.6
      ) +
      geom_point(
        data = filter(plot_tbl_subset, point_type == "Central IC genes", isDiv),
        aes(fill = point_type, shape = "Divergent IC genes"),
        color = "black", alpha = 0.9, stroke = 0.6
      ) +
      ggrepel::geom_text_repel(
        data = role_label_tbl_subset,
        aes(label = label),
        size = 5.5,
        nudge_x = role_label_tbl_subset$Shift_x,
        nudge_y = role_label_tbl_subset$Shift_y,
        point.padding = 0.25,
        box.padding = 0.35,
        force = 3.5,
        force_pull = 0.05,
        max.time = 10,
        max.iter = 10000,
        direction = "both",
        segment.color = "grey35",
        segment.size = 0.45,
        segment.alpha = 0.85,
        max.overlaps = Inf,
        min.segment.length = 0,
        seed = 1
      ) +
      scale_color_manual(values = point_colours) +
      scale_fill_manual(values = point_colours, guide = "none") +
      scale_shape_manual(name = NULL, values = c("Divergent IC genes" = 21)) +
      scale_size_continuous(
        name = "kWithin",
        range = c(3, 10)
      ) +
      guides(
        size = guide_legend(
          title = "kWithin",
          ncol = 1,
          keyheight = grid::unit(10, "mm")
        )
      ) +
      facet_wrap(~ME, scales = "fixed", ncol = ncol, nrow = nrow) +
      labs(x = "|MM|", y = "GS", color = NULL, size = "kWithin") +
  theme_bw(base_family = "Arial", base_size = 18) +
    theme(
      axis.text.x = element_text(
        family = "Arial", face = "plain", size = 18,
        colour = "black"
      ),
      axis.text.y = element_text(
        family = "Arial", face = "plain", size = 18,
        colour = "black"
      ),
      axis.title.x = element_text(
        family = "Arial", face = "plain", size = 20
      ),
      axis.title.y = element_text(
        family = "Arial", face = "plain", size = 20
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        family = "Arial", face = "plain", size = 20,
        colour = "black"
      ),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
      legend.box.just = "left",
      legend.title = element_text(
        family = "Arial", face = "plain", size = 20
      ),
      legend.text = element_text(
        family = "Arial", face = "plain", size = 18
      ),
        legend.key = element_blank(),
        legend.spacing.x = grid::unit(4, "mm"),
        panel.grid.major = element_line(
          colour = "#E6E6E6", linewidth = 0.35
        ),
        panel.grid.minor = element_blank(),
        plot.margin = margin(12, 90, 12, 12)
      )
  }

  # Split the six main modules and two supplementary modules into independent figures.
  main_modules <- c(
    "MEpurple",
    "MEturquoise",
    "MEred",
    "MEgreen",
    "MEpink",
    "MEbrown"
  )

  separate_modules <- c(
    "MEgrey60",
    "MEmagenta"
  )
  main_modules <- intersect(main_modules, as.character(significant_modules))
  separate_modules <- intersect(separate_modules, as.character(significant_modules))
  if (length(main_modules)) {
    p_role_main <- make_role_space_plot(main_modules, ncol = 3, nrow = 2)
    ggsave(file.path(run_dir, "2D_role_space_IC_6modules.pdf"), p_role_main,
           width = 17, height = 10,device = grDevices::cairo_pdf)
  }

  if (length(separate_modules)) {
    p_role_separate <- make_role_space_plot(separate_modules, ncol = 2, nrow = 1)
    ggsave(file.path(run_dir, "2D_role_space_IC_2suppleME.pdf"),
           p_role_separate,  width = 19 * 2 / 3, height = 10.5 / 2,
           device = grDevices::cairo_pdf)
  }

  # If enabled, use module genes as foreground and all network genes as the universe.
  run_functional_enrichment(
    run_dir = run_dir,
    significant_modules = significant_modules,
    network_genes = colnames(dat_expr),
    colors = colors,
    stage1 = stage1,
    cfg2 = cfg2
  )

  candidate_comparison[[run_id]] <- tibble(
    run_id = run_id,
    n_nominal_modules = length(significant_modules),
    n_nominal_trait_cells = sum(mt$p < cfg2$nominal_p_cut),
    n_fdr_trait_cells = sum(mt_fdr < cfg2$fdr_cut),
    n_hubs = sum(gene_tbl$isHUB),
    n_dr_hubs = sum(gene_tbl$isHUB & gene_tbl$isDR),
    n_nominal_dr_enriched_modules = sum(dr_tbl$p_value < cfg2$nominal_p_cut),
    n_fdr_dr_enriched_modules = sum(dr_tbl$FDR < cfg2$fdr_cut)
  )

  saveRDS(
    list(run_id = run_id, module_summary = module_summary,
         gene_metrics = gene_tbl, dr_enrichment = dr_tbl),
    file.path(run_dir, "formal_results.rds")
  )
  gc(verbose = FALSE)
}

candidate_comparison <- bind_rows(candidate_comparison) %>%
  left_join(stage1$structural_summary, by = "run_id")
write.csv(candidate_comparison,
          file.path(cfg2$output_dir, "candidate_formal_output_comparison.csv"),
          row.names = FALSE)
capture.output(sessionInfo(), file = file.path(cfg2$output_dir, "sessionInfo.txt"))
message("Stage 2 complete for ", length(selected_ids), " candidate model(s).")
