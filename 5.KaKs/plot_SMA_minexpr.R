#!/usr/bin/env Rscript
# Plot the raw bootstrap statistics produced by run_SMA_minexpr.R.
# Required packages: data.table, ggplot2.
# The lineage correlation plot shows B_Ls, Pooled, and S_Ta in that order.
# Pooled values must be calculated upstream using all eligible pairs together.

# Command-line arguments and input validation
args <- commandArgs(trailingOnly = TRUE)

if ('--help' %in% args) {
  cat('Usage: Rscript plot_SMA_minexpr.R --results RESULTS [--out PLOTS] [--min-expr 0.1,0.3,1]\n')
  quit(save = 'no')
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

flags <- args[startsWith(args, '--')]
valid_flags <- c('--results', '--out', '--min-expr')

if (length(args) != 2L * length(flags) ||
    anyDuplicated(flags) ||
    any(!flags %in% valid_flags)) {
  stop('Invalid arguments. Run with --help.')
}

arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) {
    return(default)
  }
  if (i == length(args) || startsWith(args[i + 1L], '--')) {
    stop('Missing value for ', flag)
  }
  args[i + 1L]
}

root <- normalizePath(arg('--results', '.'), mustWork = TRUE)
plot_root <- arg('--out', file.path(root, 'plots'))

st <- fread(file.path(root, 'final_raw_bootstrap_statistics.csv'))
rate <- fread(file.path(root, 'MA_rates_exact_immune_pairs.csv'))

need_stats <- c(
  'minexpr', 'grouping', 'module', 'stratum', 'metric', 'rho_raw',
  'raw_CI_low', 'raw_CI_high', 'raw_q_BH'
)
need_rate <- c('pair_id', 'OG', 'Analysis_module', 'Subgenome', 'Ka', 'Ks', 'omega')

if (!all(need_stats %in% names(st))) {
  stop('Incomplete final_raw_bootstrap_statistics.csv')
}
if (!all(need_rate %in% names(rate))) {
  stop('Incomplete MA_rates_exact_immune_pairs.csv')
}
if (!nrow(st) || anyDuplicated(st, by = c('minexpr', 'grouping', 'module', 'stratum', 'metric'))) {
  stop('Empty or duplicate statistics.')
}
if (!nrow(rate) || anyDuplicated(rate$pair_id)) {
  stop('Empty or duplicate MA pairs.')
}

for (column in c('rho_raw', 'raw_CI_low', 'raw_CI_high', 'raw_q_BH')) {
  set(st, j = column, value = as.numeric(st[[column]]))
}
for (column in c('Ka', 'Ks', 'omega')) {
  set(rate, j = column, value = as.numeric(rate[[column]]))
}

if (any(!is.finite(st$minexpr) | st$minexpr <= 0)) {
  stop('Invalid threshold in statistics.')
}

thresholds <- sort(unique(st$minexpr))
requested <- arg('--min-expr')

if (!is.null(requested)) {
  selected <- suppressWarnings(as.numeric(strsplit(requested, ',', fixed = TRUE)[[1]]))
  if (!length(selected) || any(!is.finite(selected)) || any(!selected %in% thresholds)) {
    stop('Choose thresholds present in the CSV: ', paste(thresholds, collapse = ', '))
  }
  thresholds <- sort(unique(selected))
}

# Shared labels, colors, plotting style, and output helpers
labels <- c(
  A_SMA_sq = 'Baseline SMA divergence',
  B_response24 = '0-24 Response divergence',
  B_response24_48 = '24-48 Response divergence'
)

module_levels <- sort(unique(rate$Analysis_module))
rate[, Analysis_module := factor(Analysis_module, levels = module_levels)]
rate[, Subgenome := factor(Subgenome, levels = c('CmB', 'CmS'))]

dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
out <- plot_root

cols <- c(CmB = '#1e88e5', CmS = '#e31a1c')
linecols <- c(B_Ls = '#1e88e5', S_Ta = '#e31a1c')
correlation_palette <- c(
  B_Ls = linecols[['B_Ls']],
  Pooled = '#e69f00',
  S_Ta = linecols[['S_Ta']]
)

base <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = 'bottom'
  )

savefig <- function(p, name, w, h) {
  ggsave(
    filename = file.path(out, paste0(name, '.pdf')),
    plot = p,
    width = w,
    height = h,
    limitsize = FALSE
  )
}

rate_name <- function(r) {
  if (r == 'omega') {
    return('omega')
  }
  r
}

# Ka, Ks, and omega distribution plots
rate_long <- melt(
  rate,
  id.vars = c('pair_id', 'OG', 'Analysis_module', 'Subgenome'),
  measure.vars = c('Ka', 'Ks', 'omega'),
  variable.name = 'rate',
  value.name = 'value'
)
rate_long <- rate_long[is.finite(value) & value >= 0]
limits <- c(Ka = 0.4, Ks = 1.5, omega = 0.6)

rate_module <- function(rate_name_id) {
  d <- copy(rate_long[rate == rate_name_id])
  d[, module_x := match(as.character(Analysis_module), module_levels)]

  # Width-scaled halves use full data values before plot limits are applied.
  halves <- d[, {
    if (.N >= 2L && min(value) < max(value)) {
      dens <- density(value, from = min(value), to = max(value))
      direction <- if (Subgenome[1L] == 'CmB') -1 else 1
      .(
        x = c(
          module_x[1L],
          module_x[1L] + direction * 0.425 * dens$y / max(dens$y),
          module_x[1L]
        ),
        y = c(dens$x[1L], dens$x, tail(dens$x, 1L))
      )
    } else {
      .(x = numeric(), y = numeric())
    }
  }, by = .(Analysis_module, Subgenome)]

  medians <- d[, .(
    x = module_x[1L] + ifelse(Subgenome[1L] == 'CmB', -0.18, 0.18),
    value = median(value)
  ), by = .(Analysis_module, Subgenome)]

  ggplot(d, aes(x = module_x, y = value, group = Analysis_module)) +
    geom_polygon(
      data = halves,
      aes(
        x = x,
        y = y,
        fill = Subgenome,
        group = interaction(Analysis_module, Subgenome)
      ),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    geom_boxplot(
      width = 0.20,
      outlier.shape = NA,
      fill = 'white',
      color = 'black',
      alpha = 0.85,
      linewidth = 0.45
    ) +
    geom_point(
      data = medians,
      aes(x = x, y = value, color = Subgenome),
      inherit.aes = FALSE,
      shape = 95,
      size = 7,
      show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(module_levels),
      labels = module_levels,
      position = 'top'
    ) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_color_manual(values = cols, drop = FALSE) +
    base +
    labs(x = NULL, y = rate_name(rate_name_id), fill = NULL) +
    coord_cartesian(
      xlim = c(0.5, length(module_levels) + 0.5),
      ylim = c(0, limits[[rate_name_id]])
    ) +
    theme(
      axis.text.x.top = element_text(hjust = 0, vjust = 0),
      axis.text = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
}

rate_sub <- function(r) {
  d <- rate_long[rate == r]

  ggplot(d, aes(x = Subgenome, y = value, fill = Subgenome)) +
    geom_violin(width = 0.75, alpha = 0.20, color = NA) +
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.5, color = 'black') +
    scale_x_discrete(position = 'top') +
    scale_fill_manual(values = cols, drop = FALSE) +
    base +
    theme(
      axis.text.x = element_text(angle = 0),
      axis.text.x.top = element_text(vjust = 0.5, hjust = 0.5),
      legend.position = 'none'
    ) +
    labs(x = NULL, y = rate_name(r)) +
    coord_cartesian(ylim = c(0, limits[[r]]))
}

stack_plots <- function(top, bottom, name, w, h) {
  grDevices::pdf(file.path(out, paste0(name, '.pdf')), width = w, height = h)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  lay <- grid::grid.layout(2, 1, heights = grid::unit(c(0.50, 0.50), 'null'))
  grid::pushViewport(grid::viewport(layout = lay))

  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(ggplotGrob(top))
  grid::popViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.draw(ggplotGrob(bottom))
  grid::popViewport(2)
}

for (r in c('Ka', 'Ks', 'omega')) {
  if (!nrow(rate_long[rate == r])) {
    next
  }
  savefig(rate_module(r), paste0('rates_', r, '_modules_half_violin'), 7, 3.5)
  savefig(rate_sub(r), paste0('rates_', r, '_subgenomes_boxplot'), 2, 2)
}

# Bootstrap correlation plots (estimates and inference come from input CSV)
make_corr <- function(d, metric_name, title, categories, palette, threshold) {
  d[, category := factor(category, levels = categories)]
  d[, stratum := factor(stratum, levels = names(palette))]
  d[, star := fifelse(is.finite(raw_q_BH) & raw_q_BH < 0.05, '*', '')]
  d[, star_y := pmin(
    1.1,
    fifelse(is.finite(raw_CI_high), raw_CI_high, rho_raw) + 0.10
  )]

  pd <- position_dodge(width = 0.58)

  ggplot(d, aes(x = category, y = rho_raw, color = stratum, group = stratum)) +
    geom_hline(yintercept = 0, color = 'grey65', linetype = 'dashed') +
    geom_errorbar(
      aes(ymin = raw_CI_low, ymax = raw_CI_high),
      position = pd,
      width = 0.14,
      linewidth = 0.65,
      na.rm = TRUE
    ) +
    geom_point(position = pd, size = 2.7, na.rm = TRUE) +
    geom_text(
      aes(y = star_y, label = star),
      position = pd,
      size = 6,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = c(-1, 1.2)) +
    labs(x = NULL, y = 'Correlation (rho)', color = NULL) +
    base +
    theme(
      legend.position = 'none',
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
}

# The overall/Pooled row is a separate correlation calculated on combined pairs.
# It is NOT an average of the B_Ls and S_Ta correlations.
# Some upstream versions name the overall grouping 'overall_descriptive'.
prepare_lineage_corr <- function(d) {
  lineage_data <- copy(d[grouping == 'lineage'])
  pooled_data <- copy(d[grouping %in% c('overall', 'overall_descriptive', 'pooled')])

  if (nrow(pooled_data) > 1L) {
    stop('Multiple pooled correlation rows for one threshold and metric; check input statistics.')
  }

  if (nrow(pooled_data)) {
    pooled_data[, stratum := 'Pooled']
    lineage_data <- rbindlist(list(lineage_data, pooled_data), use.names = TRUE, fill = TRUE)
  } else {
    warning('No overall/Pooled correlation row in input CSV; pooled estimate cannot be plotted.')
  }

  lineage_data <- lineage_data[stratum %in% names(correlation_palette)]
  lineage_data[, category := stratum]
  lineage_data
}

for (threshold in thresholds) {
  out <- file.path(
    plot_root,
    paste0('minexpr_', format(threshold, scientific = FALSE, trim = TRUE))
  )
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  for (m in names(labels)) {
    d <- st[minexpr == threshold & metric == m]

    # Subgenome correlation: B_Ls / Pooled / S_Ta.
    lineage_data <- prepare_lineage_corr(d)
    if (nrow(lineage_data)) {
      lineage_plot <- make_corr(
        lineage_data, m, 'by lineage',
        names(correlation_palette), correlation_palette, threshold
      )
      savefig(lineage_plot, paste0('correlation_lineage_', m), 2, 2)

      if (nrow(rate_long[rate == 'Ka'])) {
        stack_plots(rate_sub('Ka'), lineage_plot, paste0('Ka_and_lineage_', m), 6, 3)
      }
    }

    # Functional modules: same B_Ls / Pooled / S_Ta palette as before.
    module_data <- copy(d[grouping %in% c('module_pooled', 'module_by_lineage')])
    if (nrow(module_data)) {
      module_data[, category := module]
      palette <- if (any(module_data$grouping == 'module_by_lineage')) {
        correlation_palette
      } else {
        correlation_palette['Pooled']
      }

      module_plot <- make_corr(
        module_data, m, '', module_levels, palette, threshold
      )
      savefig(module_plot, paste0('Ka_correlation_modules_B_pooled_S_', m), 6, 2.5)
    }
  }
}
