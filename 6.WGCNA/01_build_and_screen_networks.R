#!/usr/bin/env Rscript

# Build and screen signed bicor networks from the Cm VST matrix and curated OG table.
# Retain the top half by MAD plus DR genes, using the upstream sample set and VST values.
# Write network workspaces, candidate rankings, stability summaries and checkpoints to stage1/.
# Structural metrics guide model selection; trait and DR summaries remain descriptive.

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)

# Bound OpenMP/BLAS threads separately from WGCNA threads.
Sys.setenv(
  OMP_NUM_THREADS = "8",
  OPENBLAS_NUM_THREADS = "8",
  MKL_NUM_THREADS = "8",
  BLIS_NUM_THREADS = "8"
)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(WGCNA)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 26 )

cfg <- list(
  expression_file = file.path(RELEASE_ROOT, "1.normalization/results/01_tximport/Cm/Cm_VST_blindFALSE.tsv"),
  oginfo_file = file.path(RELEASE_ROOT, "OG_info/OG_info_260906_OK.xlsx"),
  output_dir = file.path(RELEASE_ROOT, "6.WGCNA/stage1"),

  gene_id_column = "GID",
  dr_gene_column = "Cmo",

  # The first regex capture group is time in hours.
  sample_regex = "^Cm_([0-9]+)h_D[0-9]+$",
  expected_times = c(0, 24, 48),

  mad_keep_proportions = c(0.50),
  # DR genes bypass the MAD cutoff; goodSamplesGenes may still remove zero-variance genes.
  always_retain_dr_genes = TRUE,

  # Use the same correlation and network settings in both stages.
  network_type = "signed",
  tom_type = "signed",
  cor_type = "bicor",
  max_p_outliers = 0.10,
  pearson_fallback = "individual",

  soft_threshold_powers = c(1:10, seq(12, 30, 2)),
  scale_free_r2_target = 0.80,
  network_powers = c(8, 10, 12, 14, 16, 18, 20, 22),
  deep_splits = 0:2,
  min_module_sizes = c(30, 50),
  merge_cut_height = 0.25,
  random_seed = 54321,
  # Retain a single block when the gene count permits.
  max_block_size = 20000,

  # Exact kWithin uses the same bicor adjacency as the network.
  calculate_connectivity_metrics = TRUE,
  write_per_run_diagnostics = TRUE,

  # These downstream thresholds do not define structural quality.
  nominal_p_cut = 0.05,
  fdr_cut = 0.05,
  hub_mm_cut = 0.80,
  hub_kwithin_rank_cut = 0.70,

  # These flags rank candidates without excluding networks.
  max_grey_fraction = 0.20,
  max_largest_module_fraction = 0.60,
  max_floor_module_fraction = 0.50,
  min_median_abs_kme = 0.50,
  min_neighbor_ari = 0.50
)

assert_columns <- function(x, cols, object_name) {
  missing_cols <- setdiff(cols, colnames(x))
  if (length(missing_cols)) {
    stop(object_name, " is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

make_traits <- function(sample_names, sample_regex, expected_times) {
  time_text <- sub(sample_regex, "\\1", sample_names)
  matched <- grepl(sample_regex, sample_names)
  if (!all(matched)) {
    stop("Sample names not matching sample_regex: ",
         paste(sample_names[!matched], collapse = ", "))
  }
  time_num <- as.numeric(time_text)
  if (anyNA(time_num)) stop("Failed to parse time from sample names.")
  if (!all(time_num %in% expected_times)) {
    stop("Unexpected time values: ",
         paste(setdiff(unique(time_num), expected_times), collapse = ", "))
  }

  sample_info <- data.frame(
    sample = sample_names,
    time_hpi = time_num,
    replicate = ave(time_num, time_num, FUN = seq_along),
    row.names = sample_names,
    check.names = FALSE
  )

  # Encode the 0/24/48 h contrasts used for downstream trait associations.
  traits <- data.frame(
    t_trend = match(time_num, sort(expected_times)),
    early_0_vs_rest = as.numeric(time_num == expected_times[1]),
    mid_24_vs_mean_0_48 = as.numeric(time_num == expected_times[2]),
    late_48_vs_rest = as.numeric(time_num == expected_times[3]),
    row.names = sample_names,
    check.names = FALSE
  )
  list(sample_info = sample_info, traits = traits)
}

pearson_and_p <- function(x, y) {
  r <- WGCNA::cor(x, y, use = "pairwise.complete.obs")
  p <- WGCNA::corPvalueStudent(r, nSamples = nrow(x))
  list(cor = r, p = p)
}

hybrid_bicor_and_p <- function(x, y, max_p_outliers) {
  z <- WGCNA::bicorAndPvalue(
    x, y,
    use = "pairwise.complete.obs",
    robustX = TRUE,
    robustY = FALSE,       # Traits are binary/ordinal; do not robustify them.
    maxPOutliers = max_p_outliers,
    pearsonFallback = "individual"
  )
  list(cor = z$bicor, p = z$p)
}

get_colors <- function(net, genes) {
  colors <- net$colors
  if (is.null(colors)) stop("Network object has no colors.")
  if (is.null(names(colors))) names(colors) <- genes
  colors[genes]
}

get_ordered_mes <- function(net, dat_expr, colors) {
  mes <- net$MEs
  if (is.null(mes)) {
    mes <- WGCNA::moduleEigengenes(dat_expr, colors = colors)$eigengenes
  }
  WGCNA::orderMEs(mes)
}

module_connectivity_same_network <- function(dat_expr, colors, power, cfg) {
  genes <- colnames(dat_expr)
  modules <- setdiff(unique(as.character(colors)), "grey")

  purrr::map_dfr(modules, function(mod) {
    mod_genes <- genes[colors == mod]
    n <- length(mod_genes)
    if (n == 1L) {
      return(tibble(gene = mod_genes, module = mod,
                    kWithin = 0, kWithin_rank = 0, kWithin_scaled = 0))
    }

    adj <- WGCNA::adjacency(
      dat_expr[, mod_genes, drop = FALSE],
      type = cfg$network_type,
      power = power,
      corFnc = "bicor",
      corOptions = list(
        use = "pairwise.complete.obs",
        maxPOutliers = cfg$max_p_outliers,
        pearsonFallback = cfg$pearson_fallback
      )
    )
    k <- rowSums(adj) - diag(adj)
    tibble(
      gene = mod_genes,
      module = mod,
      kWithin = as.numeric(k),
      kWithin_rank = dplyr::percent_rank(k),
      kWithin_scaled = if (max(k) > 0) k / max(k) else 0
    )
  })
}

choose2 <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  sum_ij <- sum(choose2(tab))
  sum_a <- sum(choose2(rowSums(tab)))
  sum_b <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- sum_a * sum_b / total
  denom <- 0.5 * (sum_a + sum_b) - expected
  if (denom == 0) return(NA_real_)
  (sum_ij - expected) / denom
}

coclustering_jaccard <- function(a, b) {
  tab <- table(a, b)
  intersection <- sum(choose2(tab))
  pairs_a <- sum(choose2(rowSums(tab)))
  pairs_b <- sum(choose2(colSums(tab)))
  union <- pairs_a + pairs_b - intersection
  if (union == 0) NA_real_ else intersection / union
}

are_neighbor_runs <- function(a, b, grid_values) {
  if (a$keep_prop != b$keep_prop) return(FALSE)
  vals_a <- c(a$power, a$deep_split, a$min_module_size)
  vals_b <- c(b$power, b$deep_split, b$min_module_size)
  changed <- which(vals_a != vals_b)
  if (length(changed) != 1L) return(FALSE)
  grids <- list(grid_values$power, grid_values$deep_split,
                grid_values$min_module_size)
  g <- sort(unique(grids[[changed]]))
  abs(match(vals_a[changed], g) - match(vals_b[changed], g)) == 1L
}

plot_soft_threshold <- function(fit, keep_label, r2_target, file) {
  grDevices::pdf(file, width = 8, height = 4)
  old <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
  on.exit({par(old); grDevices::dev.off()}, add = TRUE)
  plot(fit$Power, fit$signed_SFT_R2, type = "n", xlab = "Soft-threshold power",
       ylab = "Scale-free topology fit (signed R2)", main = keep_label)
  text(fit$Power, fit$signed_SFT_R2, labels = fit$Power, col = "blue")
  abline(h = r2_target, col = "red", lty = 2)
  plot(fit$Power, fit$mean.k., type = "n", xlab = "Soft-threshold power",
       ylab = "Mean connectivity", main = keep_label)
  text(fit$Power, fit$mean.k., labels = fit$Power, col = "blue")
}

write_run_diagnostics <- function(run, dat_expr, traits, out_dir, cfg) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  net <- run$net
  colors <- get_colors(net, colnames(dat_expr))
  mes <- get_ordered_mes(net, dat_expr, colors)
  mt <- pearson_and_p(mes, traits)
  fdr <- matrix(p.adjust(mt$p, method = "BH"), nrow = nrow(mt$p),
                dimnames = dimnames(mt$p))
  stars <- matrix("", nrow(mt$p), ncol(mt$p), dimnames = dimnames(mt$p))
  stars[mt$p < cfg$nominal_p_cut] <- "*"
  stars[fdr < cfg$fdr_cut] <- "**"

  grDevices::pdf(file.path(out_dir, "dendrograms.pdf"), width = 14, height = 6)
  for (block in seq_along(net$dendrograms)) {
    block_genes <- net$blockGenes[[block]]
    WGCNA::plotDendroAndColors(
      net$dendrograms[[block]], colors[block_genes],
      main = paste(run$run_id, "block", block),
      cex.dendro = 0.5, cex.colorLabels = 0.5
    )
  }
  grDevices::dev.off()

  grDevices::pdf(file.path(out_dir, "module_trait_heatmap_Pearson.pdf"),
                 width = 6, height = max(5, 0.22 * nrow(mt$cor) + 2))
  WGCNA::labeledHeatmap(
    Matrix = mt$cor,
    xLabels = colnames(mt$cor), yLabels = rownames(mt$cor),
    ySymbols = rownames(mt$cor), colorLabels = FALSE,
    colors = colorRampPalette(c("blue", "white", "red"))(100),
    textMatrix = stars, zlim = c(-1, 1),
    main = paste(run$run_id, "module-trait relationships")
  )
  grDevices::dev.off()

  if (ncol(mes) > 1) {
    me_cor <- WGCNA::bicor(
      mes,
      use = "pairwise.complete.obs",
      maxPOutliers = cfg$max_p_outliers
    )

    grDevices::pdf(
      file.path(out_dir, "module_eigengene_clustering.pdf"),
      width = 10, height = 7
    )
    plot(
      hclust(as.dist(1 - me_cor), method = "average"),
      main = run$run_id
    )
    abline(h = cfg$merge_cut_height, col = "red", lty = 2)
    grDevices::dev.off()
  }
}

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$output_dir, "soft_threshold"), showWarnings = FALSE)
dir.create(file.path(cfg$output_dir, "run_diagnostics"), showWarnings = FALSE)

vst_tbl <- data.table::fread(cfg$expression_file, data.table = FALSE)
assert_columns(vst_tbl, cfg$gene_id_column, "VST expression table")
if (anyDuplicated(vst_tbl[[cfg$gene_id_column]])) {
  stop("Duplicated gene IDs in VST expression table.")
}

rownames(vst_tbl) <- vst_tbl[[cfg$gene_id_column]]
vst_tbl[[cfg$gene_id_column]] <- NULL
if (any(!vapply(vst_tbl, is.numeric, logical(1)))) {
  stop("All VST expression columns must be numeric.")
}

expr <- as.matrix(vst_tbl)
storage.mode(expr) <- "double"
if (any(!is.finite(expr))) {
  stop("VST expression matrix contains NA, NaN, or infinite values.")
}
if (ncol(expr) < 15) {
  warning(
    "Only ", ncol(expr),
    " samples are available. WGCNA results should be treated as exploratory; ",
    "sample QC and parameter sensitivity are especially important."
  )
}

parsed <- make_traits(colnames(expr), cfg$sample_regex, cfg$expected_times)
sample_info <- parsed$sample_info
traits <- parsed$traits

oginfo <- NULL
dr_genes <- character()
if (!is.null(cfg$oginfo_file) && file.exists(cfg$oginfo_file)) {
  oginfo <- as.data.frame(
    readxl::read_excel(cfg$oginfo_file, sheet = 1),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  assert_columns(oginfo, cfg$dr_gene_column, "oginfo table")
  dr_genes <- unique(na.omit(as.character(oginfo[[cfg$dr_gene_column]])))
}

message(
  "Upstream VST matrix loaded directly: ",
  nrow(expr), " genes x ", ncol(expr), " samples."
)

# Sample QC is descriptive; sample exclusions must be made upstream.
sample_cor_pearson <- stats::cor(expr, method = "pearson", use = "pairwise.complete.obs")
sample_cor_bicor <- WGCNA::bicor(
  expr,
  use = "pairwise.complete.obs",
  maxPOutliers = cfg$max_p_outliers,
  pearsonFallback = cfg$pearson_fallback
)
diag(sample_cor_pearson) <- NA_real_
diag(sample_cor_bicor) <- NA_real_

pca <- stats::prcomp(t(expr), center = TRUE, scale. = FALSE)
pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
sample_qc <- sample_info %>%
  mutate(
    median_Pearson_to_other_samples = apply(sample_cor_pearson, 2, median, na.rm = TRUE),
    median_bicor_to_other_samples = apply(sample_cor_bicor, 2, median, na.rm = TRUE),
    PC1 = pca$x[rownames(sample_info), 1],
    PC2 = pca$x[rownames(sample_info), 2]
  )
write.csv(
  sample_qc,
  file.path(cfg$output_dir, "sample_QC_metrics.csv"),
  row.names = FALSE
)

grDevices::pdf(file.path(cfg$output_dir, "sample_QC_plots.pdf"), width = 9, height = 7)
sample_tree <- hclust(dist(t(expr)), method = "average")
plot(sample_tree, main = "Sample clustering from upstream VST matrix",
     xlab = "", sub = "")
plot(
  pca$x[, 1], pca$x[, 2],
  pch = 19,
  col = as.integer(factor(sample_info$time_hpi)),
  xlab = paste0("PC1 (", round(pca_var[1], 1), "%)"),
  ylab = paste0("PC2 (", round(pca_var[2], 1), "%)"),
  main = "PCA of VST expression"
)
text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.75)
legend(
  "topright", legend = sort(unique(sample_info$time_hpi)),
  col = seq_along(sort(unique(sample_info$time_hpi))), pch = 19,
  title = "hpi", bty = "n"
)
boxplot(as.data.frame(expr), las = 2, outline = FALSE,
        main = "VST expression distributions", ylab = "VST expression")
grDevices::dev.off()

mad_values <- apply(expr, 1, stats::mad)

preprocess <- list()
soft_threshold_tables <- list()
runs <- list()
structural_rows <- list()
biological_rows <- list()
module_observation_rows <- list()

for (keep_prop in cfg$mad_keep_proportions) {
  keep_key <- paste0("keep_", format(keep_prop, trim = TRUE))
  mad_cut <- unname(quantile(mad_values, probs = 1 - keep_prop, na.rm = TRUE))
  retain <- mad_values >= mad_cut
  if (cfg$always_retain_dr_genes && length(dr_genes)) {
    retain[names(mad_values) %in% dr_genes] <- TRUE
  }
  expr_var <- expr[retain, , drop = FALSE]
  dat_expr <- t(expr_var)
  gsg <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  bad_samples <- rownames(dat_expr)[!gsg$goodSamples]
  if (length(bad_samples) > 0) {
    stop(
      "goodSamplesGenes() flagged input sample(s); automatic sample removal is not allowed: ",
      paste(bad_samples, collapse = ", ")
    )
  }
  dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  traits_use <- traits[rownames(dat_expr), , drop = FALSE]
  sample_info_use <- sample_info[rownames(dat_expr), , drop = FALSE]

  preprocess[[keep_key]] <- list(
    keep_prop = keep_prop,
    mad_cut = mad_cut,
    datExpr = dat_expr,
    traits = traits_use,
    sample_info = sample_info_use
  )

  sft <- WGCNA::pickSoftThreshold(
    dat_expr,
    powerVector = cfg$soft_threshold_powers,
    networkType = cfg$network_type,
    corFnc = "bicor",
    corOptions = list(
      use = "pairwise.complete.obs",
      maxPOutliers = cfg$max_p_outliers,
      pearsonFallback = cfg$pearson_fallback
    ),
    moreNetworkConcepts = TRUE,
    verbose = 5
  )
  fit <- as.data.frame(sft$fitIndices)
  # A scale-free fit requires a negative log(k)-log(p(k)) slope.
  fit$signed_SFT_R2 <- -sign(fit$slope) * fit$SFT.R.sq
  fit$keep_prop <- keep_prop
  soft_threshold_tables[[keep_key]] <- fit
  plot_soft_threshold(
    fit, keep_key, cfg$scale_free_r2_target,
    file.path(cfg$output_dir, "soft_threshold",
              paste0(safe_name(keep_key), "_soft_threshold.pdf"))
  )

  for (power in cfg$network_powers) {
    for (deep_split in cfg$deep_splits) {
      for (min_module_size in cfg$min_module_sizes) {
        run_id <- paste(
          keep_prop,
          power,
          deep_split,
          min_module_size,
          sep = "_"
        )

        run_dir <- file.path(
          cfg$output_dir,
          "run_diagnostics",
          safe_name(run_id)
        )
        checkpoint_file <- file.path(run_dir, "checkpoint.rds")
        success_file <- file.path(run_dir, "_SUCCESS")

        # Resume only from a complete saved network and its diagnostics.
        if (file.exists(checkpoint_file)) {
          checkpoint <- tryCatch(
            readRDS(checkpoint_file),
            error = function(e) {
              warning("Could not read checkpoint for ", run_id, ": ", conditionMessage(e))
              NULL
            }
          )

          required_checkpoint_items <- c(
            "run", "structural_row", "biological_row", "module_obs"
          )

          if (!is.null(checkpoint) &&
              all(required_checkpoint_items %in% names(checkpoint))) {
            runs[[run_id]] <- checkpoint$run
            structural_rows[[run_id]] <- checkpoint$structural_row
            biological_rows[[run_id]] <- checkpoint$biological_row
            module_observation_rows[[run_id]] <- checkpoint$module_obs

            # Diagnostics precede the atomic checkpoint save; restore a missing success marker.
            if (!file.exists(success_file)) {
              writeLines(
                paste("Completed:", Sys.time()),
                success_file
              )
            }

            message("Loaded completed ", run_id, " from checkpoint")
            rm(checkpoint)
            next
          }

          warning("Checkpoint for ", run_id, " is incomplete; rebuilding this run.")
          rm(checkpoint)
        }

        # A success marker without a checkpoint cannot restore the network.
        if (file.exists(success_file) && !file.exists(checkpoint_file)) {
          message("Ignoring legacy _SUCCESS without checkpoint for ", run_id)
          unlink(success_file)
        }

        message("Building ", run_id)
        set.seed(cfg$random_seed)

        net <- WGCNA::blockwiseModules(
          dat_expr,
          power = power,
          corType = cfg$cor_type,
          maxPOutliers = cfg$max_p_outliers,
          pearsonFallback = cfg$pearson_fallback,
          networkType = cfg$network_type,
          TOMType = cfg$tom_type,
          minModuleSize = min_module_size,
          deepSplit = deep_split,
          mergeCutHeight = cfg$merge_cut_height,
          numericLabels = FALSE,
          pamRespectsDendro = TRUE,
          saveTOMs = FALSE,
          maxBlockSize = cfg$max_block_size,
          randomSeed = cfg$random_seed,
          useCorOptionsThroughout = TRUE,
          verbose = 3
        )

        colors <- get_colors(net, colnames(dat_expr))
        mes <- get_ordered_mes(net, dat_expr, colors)
        tab <- table(colors)
        assigned_sizes <- as.numeric(tab[names(tab) != "grey"])
        assigned_n <- sum(assigned_sizes)
        n_modules <- length(assigned_sizes)

        # Assess gene-to-eigengene coherence with robust correlations.
        kme_bicor <- WGCNA::bicor(
          dat_expr, mes,
          use = "pairwise.complete.obs",
          maxPOutliers = cfg$max_p_outliers,
          pearsonFallback = cfg$pearson_fallback
        )
        own_me <- paste0("ME", colors)
        own_kme_bicor <- vapply(seq_along(colors), function(j) {
          if (colors[j] == "grey" || !own_me[j] %in% colnames(kme_bicor)) return(NA_real_)
          kme_bicor[j, own_me[j]]
        }, numeric(1))

        connectivity <- NULL
        if (cfg$calculate_connectivity_metrics) {
          connectivity <- module_connectivity_same_network(
            dat_expr, colors, power, cfg
          )
        }

        me_max_cor <- NA_real_
        me_pairs_over_merge <- 0L
        if (ncol(mes) > 1) {
          me_cor <- WGCNA::bicor(
            mes, use = "pairwise.complete.obs",
            maxPOutliers = cfg$max_p_outliers,
            pearsonFallback = cfg$pearson_fallback
          )
          vals <- abs(me_cor[upper.tri(me_cor)])
          me_max_cor <- max(vals, na.rm = TRUE)
          me_pairs_over_merge <- sum(vals > (1 - cfg$merge_cut_height), na.rm = TRUE)
        }

        sft_row <- fit[match(power, fit$Power), , drop = FALSE]
        structural_rows[[run_id]] <- tibble(
          run_id = run_id,
          keep_prop = keep_prop,
          power = power,
          deep_split = deep_split,
          min_module_size = min_module_size,
          n_genes = ncol(dat_expr),
          n_modules = n_modules,
          grey_fraction = if ("grey" %in% names(tab)) unname(tab["grey"]) / length(colors) else 0,
          median_module_size = if (n_modules) median(assigned_sizes) else NA_real_,
          min_module_size_observed = if (n_modules) min(assigned_sizes) else NA_real_,
          max_module_size_observed = if (n_modules) max(assigned_sizes) else NA_real_,
          largest_module_fraction = if (assigned_n) max(assigned_sizes) / assigned_n else NA_real_,
          floor_module_fraction = if (n_modules) mean(assigned_sizes <= 1.25 * min_module_size) else NA_real_,
          median_abs_kme_bicor = median(abs(own_kme_bicor), na.rm = TRUE),
          prop_abs_kme_bicor_ge_0.5 = mean(abs(own_kme_bicor) >= 0.5, na.rm = TRUE),
          median_kwithin_scaled = if (!is.null(connectivity)) median(connectivity$kWithin_scaled, na.rm = TRUE) else NA_real_,
          median_module_density = if (!is.null(connectivity)) {
            connectivity %>% group_by(module) %>%
              summarise(density = mean(kWithin) / max(n() - 1, 1), .groups = "drop") %>%
              summarise(x = median(density, na.rm = TRUE)) %>% pull(x)
          } else NA_real_,
          max_abs_me_correlation = me_max_cor,
          n_me_pairs_above_merge_correlation = me_pairs_over_merge,
          scale_free_r2 = if (nrow(sft_row)) sft_row$signed_SFT_R2 else NA_real_,
          mean_connectivity_sft = if (nrow(sft_row)) sft_row$mean.k. else NA_real_
        )

        # Trait and DR summaries do not enter the structural score.
        mt_pearson <- pearson_and_p(mes, traits_use)
        mt_bicor <- hybrid_bicor_and_p(mes, traits_use, cfg$max_p_outliers)
        fdr_pearson <- p.adjust(as.vector(mt_pearson$p), method = "BH")
        fdr_bicor <- p.adjust(as.vector(mt_bicor$p), method = "BH")
        nominal_modules <- rownames(mt_pearson$p)[apply(mt_pearson$p < cfg$nominal_p_cut, 1, any)]
        nominal_colors <- sub("^ME", "", nominal_modules)

        dr_in_network <- intersect(names(colors), dr_genes)
        dr_in_nominal <- dr_in_network[colors[dr_in_network] %in% nominal_colors]

        dr_tests <- tibble(
          module = setdiff(names(tab), "grey"),
          module_size = as.numeric(tab[setdiff(names(tab), "grey")]),
          DR_in_module = vapply(setdiff(names(tab), "grey"), function(mod) {
            sum(names(colors)[colors == mod] %in% dr_genes)
          }, numeric(1)),
          odds_ratio = NA_real_, p = NA_real_, fdr = NA_real_
        )
        best_dr_p <- best_dr_fdr <- best_dr_or <- NA_real_
        n_dr_enriched_nominal <- n_dr_enriched_fdr <- 0L
        if (length(dr_genes) && n_modules) {
          dr_tests <- purrr::map_dfr(setdiff(names(tab), "grey"), function(mod) {
            in_mod <- names(colors)[colors == mod]
            a <- sum(in_mod %in% dr_genes)
            b <- length(in_mod) - a
            c <- sum(names(colors) %in% dr_genes) - a
            d <- length(colors) - a - b - c
            ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                              alternative = "greater")
            tibble(module = mod, module_size = length(in_mod),
                   DR_in_module = a,
                   odds_ratio = unname(ft$estimate), p = ft$p.value)
          }) %>% mutate(fdr = p.adjust(p, method = "BH"))
          best_i <- which.min(dr_tests$p)
          best_dr_p <- dr_tests$p[best_i]
          best_dr_fdr <- dr_tests$fdr[best_i]
          best_dr_or <- dr_tests$odds_ratio[best_i]
          n_dr_enriched_nominal <- sum(dr_tests$p < cfg$nominal_p_cut)
          n_dr_enriched_fdr <- sum(dr_tests$fdr < cfg$fdr_cut)
        }

        # Use bicor MM and the same adjacency for kWithin so hub rankings remain consistent.
        n_hubs <- n_dr_hubs <- NA_integer_
        if (!is.null(connectivity)) {
          kme_hub <- WGCNA::bicor(
            dat_expr, mes,
            use = "pairwise.complete.obs",
            maxPOutliers = cfg$max_p_outliers,
            pearsonFallback = cfg$pearson_fallback
          )
          own_kme_hub <- vapply(seq_along(colors), function(j) {
            if (colors[j] == "grey" || !own_me[j] %in% colnames(kme_hub)) return(NA_real_)
            kme_hub[j, own_me[j]]
          }, numeric(1))
          names(own_kme_hub) <- names(colors)
          hub_tbl <- connectivity %>%
            mutate(MM = own_kme_hub[gene],
                   is_hub = MM >= cfg$hub_mm_cut &
                     kWithin_rank >= cfg$hub_kwithin_rank_cut)
          n_hubs <- sum(hub_tbl$is_hub, na.rm = TRUE)
          n_dr_hubs <- sum(hub_tbl$is_hub & hub_tbl$gene %in% dr_genes, na.rm = TRUE)
        }

        biological_rows[[run_id]] <- tibble(
          run_id = run_id,
          n_nominal_trait_cells_pearson = sum(mt_pearson$p < cfg$nominal_p_cut),
          n_fdr_trait_cells_pearson = sum(fdr_pearson < cfg$fdr_cut),
          n_nominal_trait_modules_pearson = length(nominal_modules),
          n_nominal_trait_cells_hybrid_bicor = sum(mt_bicor$p < cfg$nominal_p_cut),
          n_fdr_trait_cells_hybrid_bicor = sum(fdr_bicor < cfg$fdr_cut),
          n_dr_genes_in_network = length(dr_in_network),
          n_dr_genes_in_nominal_modules = length(dr_in_nominal),
          best_dr_enrichment_odds_ratio = best_dr_or,
          best_dr_enrichment_p = best_dr_p,
          best_dr_enrichment_fdr = best_dr_fdr,
          n_dr_enriched_modules_nominal = n_dr_enriched_nominal,
          n_dr_enriched_modules_fdr = n_dr_enriched_fdr,
          n_hubs = n_hubs,
          n_dr_hubs = n_dr_hubs
        )

        trait_obs <- purrr::map_dfr(rownames(mt_pearson$p), function(me) {
          pvals <- mt_pearson$p[me, ]
          rs <- mt_pearson$cor[me, ]
          best <- which.min(pvals)
          tibble(
            module = sub("^ME", "", me),
            any_nominal_trait = any(pvals < cfg$nominal_p_cut, na.rm = TRUE),
            strongest_trait = colnames(mt_pearson$p)[best],
            strongest_trait_correlation = rs[best],
            minimum_trait_p = pvals[best]
          )
        })
        module_obs <- trait_obs %>%
          left_join(dr_tests, by = "module") %>%
          mutate(run_id = run_id, .before = 1)
        module_observation_rows[[run_id]] <- module_obs

        if (cfg$write_per_run_diagnostics) {
          diag_dir <- file.path(cfg$output_dir, "run_diagnostics", safe_name(run_id))
          dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
          write.csv(module_obs, file.path(diag_dir, "module_observations.csv"),
                    row.names = FALSE)
          nominal_dr_counts <- module_obs$DR_in_module[module_obs$any_nominal_trait]
          if (length(nominal_dr_counts) && any(is.finite(nominal_dr_counts))) {
            grDevices::pdf(file.path(diag_dir,
                                     "DR_counts_in_nominal_trait_modules.pdf"),
                           width = 6, height = 4)
            hist(nominal_dr_counts, breaks = "FD",
                 main = run_id, xlab = "DR genes per nominal trait-associated module")
            grDevices::dev.off()
          }
        }

        runs[[run_id]] <- list(
          run_id = run_id,
          keep_key = keep_key,
          keep_prop = keep_prop,
          power = power,
          deep_split = deep_split,
          min_module_size = min_module_size,
          net = net
        )

        if (cfg$write_per_run_diagnostics) {
          write_run_diagnostics(
            runs[[run_id]], dat_expr, traits_use,
            run_dir, cfg
          )
        }

        # Rename the completed temporary checkpoint atomically to avoid partial saves.
        dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
        checkpoint <- list(
          run = runs[[run_id]],
          structural_row = structural_rows[[run_id]],
          biological_row = biological_rows[[run_id]],
          module_obs = module_observation_rows[[run_id]]
        )

        tmp_checkpoint <- file.path(run_dir, "checkpoint.rds.tmp")
        checkpoint_file <- file.path(run_dir, "checkpoint.rds")
        success_file <- file.path(run_dir, "_SUCCESS")

        if (file.exists(tmp_checkpoint)) unlink(tmp_checkpoint)
        saveRDS(checkpoint, tmp_checkpoint, compress = FALSE)

        if (file.exists(checkpoint_file)) unlink(checkpoint_file)
        if (!file.rename(tmp_checkpoint, checkpoint_file)) {
          stop("Failed to finalize checkpoint for ", run_id)
        }

        writeLines(
          paste("Completed:", Sys.time()),
          success_file
        )

        rm(checkpoint)
        rm(net, connectivity, kme_bicor)
        gc(verbose = FALSE)
      }
    }
  }
}

# Compare neighboring models using label-invariant module assignments.

run_meta <- bind_rows(lapply(runs, function(x) {
  tibble(run_id = x$run_id, keep_prop = x$keep_prop, power = x$power,
         deep_split = x$deep_split, min_module_size = x$min_module_size)
}))

pairs <- combn(names(runs), 2, simplify = FALSE)
pairwise_stability <- purrr::map_dfr(pairs, function(ids) {
  a <- runs[[ids[1]]]
  b <- runs[[ids[2]]]
  da <- preprocess[[a$keep_key]]$datExpr
  db <- preprocess[[b$keep_key]]$datExpr
  ca <- get_colors(a$net, colnames(da))
  cb <- get_colors(b$net, colnames(db))
  common <- intersect(names(ca), names(cb))
  common <- common[ca[common] != "grey" & cb[common] != "grey"]
  if (length(common) < 2) return(tibble())

  meta_a <- run_meta[match(ids[1], run_meta$run_id), ]
  meta_b <- run_meta[match(ids[2], run_meta$run_id), ]
  neighbor <- are_neighbor_runs(
    meta_a, meta_b,
    list(power = cfg$network_powers, deep_split = cfg$deep_splits,
         min_module_size = cfg$min_module_sizes)
  )
  tibble(
    run_a = ids[1], run_b = ids[2], n_common_assigned = length(common),
    adjusted_rand_index = adjusted_rand_index(ca[common], cb[common]),
    coclustering_jaccard = coclustering_jaccard(ca[common], cb[common]),
    same_keep_prop = a$keep_prop == b$keep_prop,
    is_neighbor = neighbor
  )
})

neighbor_summary <- bind_rows(
  pairwise_stability %>% transmute(run_id = run_a, other = run_b,
                                   same_keep_prop, is_neighbor,
                                   adjusted_rand_index, coclustering_jaccard),
  pairwise_stability %>% transmute(run_id = run_b, other = run_a,
                                   same_keep_prop, is_neighbor,
                                   adjusted_rand_index, coclustering_jaccard)
) %>%
  group_by(run_id) %>%
  summarise(
    n_neighbor_models = sum(is_neighbor),
    median_neighbor_ari = median(adjusted_rand_index[is_neighbor], na.rm = TRUE),
    min_neighbor_ari = min(adjusted_rand_index[is_neighbor], na.rm = TRUE),
    median_neighbor_coclustering_jaccard =
      median(coclustering_jaccard[is_neighbor], na.rm = TRUE),
    mean_ari_all_same_keep_prop =
      mean(adjusted_rand_index[same_keep_prop], na.rm = TRUE),
    .groups = "drop"
  )

structural_summary <- bind_rows(structural_rows) %>%
  left_join(neighbor_summary, by = "run_id") %>%
  mutate(
    scale_free_flag = !is.na(scale_free_r2) & scale_free_r2 >= cfg$scale_free_r2_target,
    grey_flag = grey_fraction <= cfg$max_grey_fraction,
    largest_module_flag = largest_module_fraction <= cfg$max_largest_module_fraction,
    fragmentation_flag = floor_module_fraction <= cfg$max_floor_module_fraction,
    coherence_flag = median_abs_kme_bicor >= cfg$min_median_abs_kme,
    stability_flag = median_neighbor_ari >= cfg$min_neighbor_ari,
    n_structural_flags_passed = rowSums(across(ends_with("_flag")), na.rm = TRUE),
    all_structural_flags_passed = scale_free_flag & grey_flag &
      largest_module_flag & fragmentation_flag & coherence_flag & stability_flag
  ) %>%
  arrange(desc(all_structural_flags_passed), desc(n_structural_flags_passed),
          desc(median_neighbor_ari), desc(median_abs_kme_bicor), grey_fraction,
          floor_module_fraction)

biological_observations <- bind_rows(biological_rows)
module_biological_observations <- bind_rows(module_observation_rows)
overview <- structural_summary %>%
  left_join(biological_observations, by = "run_id")

candidate_models <- structural_summary %>%
  transmute(
    use_in_stage2 = FALSE,
    run_id, keep_prop, power, deep_split, min_module_size,
    all_structural_flags_passed, n_structural_flags_passed,
    median_neighbor_ari, median_neighbor_coclustering_jaccard,
    grey_fraction, n_modules, median_abs_kme_bicor,
    note = "Set use_in_stage2 to TRUE for the shortlisted models"
  )

soft_threshold_summary <- bind_rows(soft_threshold_tables)
write.csv(soft_threshold_summary,
          file.path(cfg$output_dir, "soft_threshold_summary.csv"), row.names = FALSE)
write.csv(structural_summary,
          file.path(cfg$output_dir, "model_structural_summary.csv"), row.names = FALSE)
write.csv(biological_observations,
          file.path(cfg$output_dir, "model_biological_observations.csv"), row.names = FALSE)
write.csv(module_biological_observations,
          file.path(cfg$output_dir, "module_biological_observations.csv"),
          row.names = FALSE)
write.csv(pairwise_stability,
          file.path(cfg$output_dir, "pairwise_model_stability.csv"), row.names = FALSE)
write.csv(overview,
          file.path(cfg$output_dir, "model_selection_overview.csv"), row.names = FALSE)
write.csv(candidate_models,
          file.path(cfg$output_dir, "candidate_models.csv"), row.names = FALSE)
write.csv(sample_info,
          file.path(cfg$output_dir, "sample_information.csv"), row.names = FALSE)

workspace <- list(
  config = cfg,
  sample_info = sample_info,
  traits = traits,
  oginfo = oginfo,
  dr_genes = dr_genes,
  preprocess = preprocess,
  soft_threshold = soft_threshold_tables,
  runs = runs,
  structural_summary = structural_summary,
  biological_observations = biological_observations,
  module_biological_observations = module_biological_observations,
  pairwise_stability = pairwise_stability
)
saveRDS(workspace, file.path(cfg$output_dir, "stage1_workspace.rds"), compress = FALSE)

capture.output(sessionInfo(), file = file.path(cfg$output_dir, "sessionInfo.txt"))
message("Stage 1 complete. Edit candidate_models.csv, then run stage 2.")
