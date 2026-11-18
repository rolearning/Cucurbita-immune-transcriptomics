#!/usr/bin/env Rscript
# Ka–expression associations with whole-orthogroup bootstrap.
# Required packages: data.table, readxl. Run with --help for usage.
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], '--')) stop('Missing value for ', flag)
  args[i + 1L]
}
if ('--help' %in% args) {
  cat('Usage: Rscript run_SMA_minexpr.R --root PROJECT --out RESULTS --min-expr 0.1,0.3,1 --background-mode multi\n',
  'Optional: --boot 1000 --min-bg 50 --min-og-test 5 --contrast-boot 1000 --broad-modules no --module-subgenome yes\n',
  'Inputs can be overridden with --synteny, --subgenome-map, --module-map, --expr-cm, --expr-ls, --expr-ta and --factors.\n',
  sep = '')
  quit(save = 'no')
}
known <- c('--root', '--out', '--min-expr', '--background-mode', '--boot', '--min-bg',
  '--min-og-test', '--contrast-boot', '--broad-modules', '--module-subgenome',
  '--synteny', '--subgenome-map', '--module-map', '--expr-cm', '--expr-ls', '--expr-ta', '--factors')
flags <- args[startsWith(args, '--')]
if (length(args) != 2L * length(flags) || anyDuplicated(flags) || any(!flags %in% known))
  stop('Invalid arguments. Run with --help. Unknown options: ', paste(setdiff(flags, known), collapse = ', '))
minexpr_grid <- suppressWarnings(as.numeric(strsplit(arg('--min-expr', '0.1,0.3,1'), ',', fixed = TRUE)[[1]]))
if (!length(minexpr_grid) || any(!is.finite(minexpr_grid) | minexpr_grid <= 0))
  stop('--min-expr must contain positive thresholds, e.g. 0.1,0.3,1')
minexpr_grid <- sort(unique(minexpr_grid))
int_arg <- function(flag, default, minimum) {
  value <- suppressWarnings(as.numeric(arg(flag, default)))
  if (!is.finite(value) || value < minimum || value != floor(value) || value > .Machine$integer.max)
    stop(flag, ' must be an integer >= ', minimum)
  as.integer(value)
}
minbg <- int_arg('--min-bg', '50', 2L)
bootN <- int_arg('--boot', '1000', 100L)
minOG <- int_arg('--min-og-test', '5', 5L)
contrastN <- int_arg('--contrast-boot', '0', 0L)
if (contrastN > 0L && contrastN < 200L) stop('--contrast-boot must be 0 or >= 200')

root <- normalizePath(arg('--root', '.'), mustWork = TRUE)
out <- arg('--out', file.path(root, 'SMA_raw_minexpr'))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
getfile <- function(n) {
  p <- c(file.path(root, n), file.path(root, 'data', n))
  m <- p[file.exists(p)]
  if (length(m)) m[1] else NA_character_
}
syfile <- arg('--synteny', NA_character_)
if (is.na(syfile)) {
  for (n in c(
  'All_single_syngroup_w_notes.xlsx',
  'All_single_syngroup_w_notes.csv',
  'All_single_syngroup_w_notes.tsv'
  )) {
    if (is.na(syfile) && !is.na(getfile(n))) syfile <- getfile(n)
  }
}
paths <- c(
  ka = getfile('KaKs_table_revision1.xlsx'),
  og = getfile('OG_table_revision1_1025.xlsx'),
  synteny = syfile,
  factors = getfile('all_candidate_sample_factors.tsv')
)
for (sp in c('Cm', 'Ls', 'Ta')) {
  paths[sp] <- arg(paste0('--expr-', tolower(sp)), getfile(paste0(sp, '_scaledTPM_counts.tsv')))
}
paths['factors'] <- arg('--factors', paths['factors'])
if (anyNA(paths) || any(!file.exists(paths))) {
  stop('Missing inputs: ', paste(names(paths)[is.na(paths) | !file.exists(paths)], collapse = ', '))
}
read_any <- function(p) {
  if (grepl('\\.xlsx?$', p, ignore.case = TRUE)) as.data.table(read_excel(p)) else fread(p)
}
save <- function(d, n) {
  fwrite(as.data.table(d), file.path(out, paste0(n, '.csv')), sep = ',', na = 'NA')
}
clean <- function(z) {
  z <- trimws(as.character(z))
  z[z == '' | toupper(z) %in% c('NA', 'NULL')] <- NA_character_
  sub('^TAN', 'Tan', z)
}
subcanon <- function(s) {
  s <- toupper(gsub('[ _-]', '', clean(s)))
  fifelse(s %in% c('B', 'CMB', 'BENI', 'BENINCASEAE'),
  'CmB',
  fifelse(s %in% c('S', 'CMS', 'SICY', 'SICYOEAE'), 'CmS', NA_character_))
}
optname <- function(d, choices) {
  n <- names(d)
  hit <- which(tolower(n) %in% tolower(choices))
  if (length(hit)) n[hit[1]] else NA_character_
}
# Scale all count-like expression using the prior validated all_BS/TMM multipliers
fac <- fread(paths['factors'])
stopifnot(all(c('reference_set', 'method', 'sample', 'multiplier') %in% names(fac)))
fac <- fac[reference_set == 'all_BS' & method == 'TMM']
if (!nrow(fac) || anyDuplicated(fac$sample)) stop('Invalid all_BS / TMM factors')
read_expression <- function(sp) {
  d <- fread(paths[sp])
  if (!'GID' %in% names(d)) stop('Missing GID in ', sp)
  setnames(d, 'GID', 'gene')
  d[, gene := clean(gene)]
  if (anyDuplicated(d$gene)) stop('Duplicate expression gene IDs ', sp)
  cols <- setdiff(names(d), 'gene')
  v <- fac$multiplier[match(cols, fac$sample)]
  if (any(!is.finite(v)) || any(v <= 0)) stop('Missing or nonpositive multipliers for ', sp)
  for (i in seq_along(cols)) set(d, j = cols[i], value = as.numeric(d[[cols[i]]]) * v[i])
  tt <- regmatches(cols, regexec(paste0('^', sp, '_(0|24|48)h_D[0-9]+$'), cols))
  if (any(lengths(tt) != 2L)) stop('Unrecognized sample IDs for ', sp)
  tp <- vapply(tt, '[', character(1), 2)
  o <- d[, .(gene)]
  for (t in c('0', '24', '48')) {
    cc <- cols[tp == t]
    if (!length(cc)) stop('No ', sp, ' at ', t)
    o[, (paste0('t', t)) := rowMeans(as.matrix(d[, ..cc]))]
  }
  o
}
means <- lapply(c('Cm', 'Ls', 'Ta'), read_expression)
names(means) <- c('Cm', 'Ls', 'Ta')
# Immune Ka-pair table, exact pair matching to approved OGs
ka <- read_any(paths['ka'])
req <- c('Sequence', 'Method', 'Ka')
if (!all(req %in% names(ka))) stop('Ka table missing ', paste(setdiff(req, names(ka)), collapse = ','))
ka <- ka[Method == 'MA']
ka[, c('Cm_gene', 'proxy_gene') := tstrsplit(as.character(Sequence), '-', fixed = TRUE)]
ka[, `:=`(Cm_gene = clean(Cm_gene), proxy_gene = clean(proxy_gene), Ka = as.numeric(Ka))]
ka <- ka[grepl('^Cmo', Cm_gene) & (grepl('^Lsi', proxy_gene) | grepl('^Tan', proxy_gene))]
ka[, lineage := fifelse(grepl('^Lsi', proxy_gene), 'B_Ls', 'S_Ta')]
og <- read_any(paths['og'])
need <- c('Cmo_Gene_ID', 'OG', 'Module', 'Subgenome', 'Progenitor_Lineage_Proxy_Gene_ID')
if (!all(need %in% names(og))) stop('OG table missing columns: ',
  paste(setdiff(need, names(og)), collapse = ','))
og[, `:=`(Cmo_Gene_ID = clean(Cmo_Gene_ID),
  OG = as.character(OG),
  Module = trimws(as.character(Module)),
  Subgenome = subcanon(Subgenome))]
og[, proxy_gene := strsplit(gsub('[;， ]+', ',', as.character(Progenitor_Lineage_Proxy_Gene_ID)),
  ',',
  fixed = TRUE)]
og <- unique(og[, .(proxy_gene = clean(unlist(proxy_gene))), by = .(Cm_gene = Cmo_Gene_ID,
  OG,
  Module,
  Subgenome)][!is.na(proxy_gene)])
ka[, ka_row := .I]
j <- merge(ka, og, by = c('Cm_gene', 'proxy_gene'), all.x = TRUE, allow.cartesian = TRUE)
save(j[is.na(OG) | is.na(Subgenome) | !(lineage == 'B_Ls' & Subgenome == 'CmB' | lineage == 'S_Ta' & Subgenome == 'CmS')], 'QC_immune_pairs_excluded')
q <- j[!is.na(OG) & ((lineage == 'B_Ls' & Subgenome == 'CmB') | (lineage == 'S_Ta' & Subgenome == 'CmS')) & is.finite(Ka) & Ka >= 0]
if (!nrow(q)) stop('No valid Ka–OG matches')
if (anyDuplicated(q$ka_row)) stop('Ka pair maps to >1 OG; check inputs')
q[, pair_id := paste(Cm_gene, proxy_gene, sep = '__')]
if (anyDuplicated(q$pair_id)) stop('Duplicate Ka pairs')
# Background: SIMPLE groups, explicit ancestry; strict one-to-one OR multi-copy weighted background
# --background-mode strict|multi: multi retains repeated proxies and weights each proxy equally within lineage.
bgmode <- tolower(arg('--background-mode', 'strict'))
if (!bgmode %in% c('strict', 'multi')) stop('--background-mode must be strict or multi')
sy <- read_any(paths['synteny'])
typ <- optname(sy, c('type'))
if (is.na(typ)) stop('Synteny table needs type (simple/complex)')
sy <- sy[tolower(trimws(as.character(get(typ)))) == 'simple']
if (!nrow(sy)) stop('No simple syngroups')
cb <- optname(sy, c('CmB'))
cs <- optname(sy, c('CmS'))
cm <- optname(sy, c('Cm'))
anc <- optname(sy, c('Subgenome', 'Cm_subgenome', 'ancestry', 'Cm_origin'))
ls <- optname(sy, c('Lsi', 'Ls'))
ta <- optname(sy, c('Tan', 'Ta'))
if (is.na(ls) || is.na(ta)) stop('Synteny requires Lsi and Tan columns')
if (!is.na(cb) && !is.na(cs)) {
  b <- data.table(Cm_gene = clean(sy[[cb]]), proxy_gene = clean(sy[[ls]]), lineage = 'B_Ls')
  s <- data.table(Cm_gene = clean(sy[[cs]]), proxy_gene = clean(sy[[ta]]), lineage = 'S_Ta')
  bg <- rbindlist(list(b, s))
} else if (!is.na(cm)) {
  if (is.na(anc)) {
    mp <- arg('--subgenome-map', NA_character_)
    if (is.na(mp) || !file.exists(mp)) stop('The screenshot columns alone do NOT identify CmB/CmS. Add Subgenome column to syngroup sheet OR provide --subgenome-map TSV with Cm,Subgenome; stopping rather than inferring ancestry.')
    m <- read_any(mp)
    mc <- optname(m, c('Cm', 'Cm_gene', 'Cmo_Gene_ID'))
    ma <- optname(m, c('Subgenome', 'ancestry'))
    if (is.na(mc) || is.na(ma)) stop('Map must have Cm and Subgenome columns')
    m <- unique(data.table(Cm_gene = clean(m[[mc]]), anc = subcanon(m[[ma]])))
    if (anyDuplicated(m$Cm_gene)) stop('Ambiguous subgenome-map genes')
    sy[, Cm_gene_tmp := clean(get(cm))]
    sy <- merge(sy, m, by.x = 'Cm_gene_tmp', by.y = 'Cm_gene', all.x = TRUE, sort = FALSE)
    anc <- 'anc'
  }
  ss <- subcanon(sy[[anc]])
  bg <- rbindlist(list(data.table(Cm_gene = clean(sy[[cm]]), proxy_gene = clean(sy[[ls]]), lineage = 'B_Ls', subgenome = ss),
  data.table(Cm_gene = clean(sy[[cm]]), proxy_gene = clean(sy[[ta]]), lineage = 'S_Ta', subgenome = ss)))
  bg <- bg[(lineage == 'B_Ls' & subgenome == 'CmB') | (lineage == 'S_Ta' & subgenome == 'CmS')][, subgenome := NULL]
} else stop('Synteny requires Cm with Subgenome OR both CmB and CmS columns')
bg <- unique(bg[!is.na(Cm_gene) & !is.na(proxy_gene)])
bg[, n_proxy := uniqueN(proxy_gene), by = .(lineage, Cm_gene)]
bg[, n_cm := uniqueN(Cm_gene), by = .(lineage, proxy_gene)]
save(bg, 'QC_background_before_reciprocal_filter')
# In this source each Cm is already unique, while Lsi/Tan proxy IDs can recur.
# Strict: reciprocal one-to-one. Multi: keep every distinct edge; weights prevent
# a reused proxy (or Cm) from contributing linearly more than a single-copy gene.
if (bgmode == 'strict') bg <- bg[n_proxy == 1 & n_cm == 1]
if (!nrow(bg)) stop('No background syngroups after ', bgmode, ' selection')
bg[, pair_weight := if (bgmode == 'strict') 1 else 1 / (n_cm * n_proxy)]
save(bg, paste0('background_simple_', bgmode))

attach <- function(d) {
  d <- copy(d)
  d[, row_id__ := .I]
  m <- copy(means$Cm)
  setnames(m, c('gene', 't0', 't24', 't48'), c('Cm_gene', 'cm_0', 'cm_24', 'cm_48'))
  d <- merge(d, m, by = 'Cm_gene', all.x = TRUE, sort = FALSE)
  for (t in c('0', '24', '48')) d[, (paste0('proxy_', t)) := NA_real_]
  for (l in c('B_Ls', 'S_Ta')) {
    ids <- which(d$lineage == l)
    if (!length(ids)) next
    sp <- if (l == 'B_Ls') 'Ls' else 'Ta'
    m <- means[[sp]]
    ix <- match(d$proxy_gene[ids], m$gene)
    for (t in c('0', '24', '48')) set(d, i = ids, j = paste0('proxy_', t), value = m[[paste0('t', t)]][ix])
  }
  setorder(d, row_id__)
  d[, row_id__ := NULL]
  d
}
q <- attach(q)
bg <- attach(bg)
save(q, 'immune_pair_expression_inputs')
q[, `:=`(Module = trimws(Module), Module_group = NA_character_)]
mapfile <- arg('--module-map', NA_character_)
if (!is.na(mapfile)) {
  mp <- read_any(mapfile)
  if (!all(c('Module',
  'Module_group') %in% names(mp))) stop('module-map must contain Module and Module_group')
  if (anyDuplicated(mp$Module)) stop('Duplicate module mapping')
  q[, Module_group := as.character(mp$Module_group[match(Module, mp$Module)])]
} else {
  z <- tolower(q$Module)
  q[, Module_group := fcase(grepl('rlck|mapk|calcium|ca2|ca\\+|signaling|signalling', z),
  'Signaling',
  grepl('hormon|salicyl|jasmon|ethylene|auxin|abscisic|^sa$|^sas$|^jas\\s*&\\s*ets$|^jas|^ets$', z),
  'Hormone',
  grepl('prr|nlr|recognition|receptor', z),
  'Recognition',
  grepl('wrky|transcription', z),
  'Transcription',
  grepl('cucurbit|secondary|specialized|biosynth', z),
  'Defense_metabolism',
  default = 'Other')]
}
q[is.na(Module_group) | trimws(Module_group) == '', Module_group := 'Other']
save(unique(q[, .(Module, Module_group)]), 'QC_module_mapping')
q[, Module_fine := trimws(as.character(Module))]
q[is.na(Module_fine) | Module_fine == '', Module_fine := 'Unassigned']
use_broad <- tolower(arg('--broad-modules', 'no')) %in% c('yes', 'true', '1')
module_by_lineage <- tolower(arg('--module-subgenome', 'yes')) %in% c('yes', 'true', '1')
q[, Analysis_module := if (use_broad) Module_group else Module_fine]
if (!use_broad && uniqueN(q$Analysis_module) != 9L)
  warning('Expected nine functional modules; found ',
  uniqueN(q$Analysis_module),
  ' (see QC_module_mapping.csv).')

fit_sma <- function(d) {
  x <- log2(d$cm_0)
  y <- log2(d$proxy_0)
  w <- d$pair_weight
  if (length(x) < minbg) stop('SMA background insufficient')
  if (any(!is.finite(w) | w <= 0)) stop('Bad SMA weights')
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  vx <- sum(w * (x - mx) ^ 2)
  vy <- sum(w * (y - my) ^ 2)
  if (vx <= 0 || vy <= 0) stop('Constant background expression')
  r <- sum(w * (x - mx) * (y - my)) / sqrt(vx * vy)
  if (!is.finite(r) || abs(r) < .05) stop('Weak background correlation; inspect')
  b <- sign(r) * sqrt(vy / vx)
  list(intercept = my - b * mx,
  slope = b,
  r = r,
  n = length(x),
  n_eff = 1 / sum(w ^ 2),
  weight_sum = sum(d$pair_weight))
}

metrics <- c('A_SMA_sq', 'B_response24', 'B_response24_48')
lineages <- c('B_Ls', 'S_Ta')
modules <- sort(unique(na.omit(q$Analysis_module)))
passes <- function(x, threshold) is.finite(x) & x >= threshold

raw_spearman <- function(x, y) {
  if (length(x) < 5L || min(x) == max(x) || min(y) == max(y)) return(NA_real_)
  cor(x, y, method = 'spearman')
}

eligible_rows <- function(d, metric_name) {
  d[is.finite(Ka) & is.finite(get(metric_name)) & !is.na(OG)]
}

cluster_estimate <- function(d, metric_name) {
  d <- eligible_rows(d, metric_name)
  blocks <- split(seq_len(nrow(d)), as.character(d$OG))
  n_og <- length(blocks)
  ans <- data.table(n_pairs = nrow(d), n_OG = n_og, n_boot_valid = 0L, rho_raw = NA_real_,
  raw_CI_low = NA_real_, raw_CI_high = NA_real_, raw_boot_SE = NA_real_,
  raw_p_cluster_normal = NA_real_)
  if (nrow(d) < 5L || n_og < minOG) return(ans)
  x <- d$Ka
  y <- d[[metric_name]]
  obs <- raw_spearman(x, y)
  ans[, rho_raw := obs]
  if (!is.finite(obs)) return(ans)
  draws <- replicate(bootN, {
    ids <- unlist(blocks[sample.int(n_og, n_og, replace = TRUE)], use.names = FALSE)
    raw_spearman(x[ids], y[ids])
  })
  draws <- draws[is.finite(draws)]
  ans[, n_boot_valid := length(draws)]
  if (length(draws) < .8 * bootN) return(ans)
  ci <- quantile(draws, c(.025, .975), names = FALSE)
  se <- sd(draws)
  p <- if (is.finite(se) && se > 0) 2 * pnorm(abs(obs / se), lower.tail = FALSE) else NA_real_
  ans[, `:=`(raw_CI_low = ci[1L], raw_CI_high = ci[2L], raw_boot_SE = se, raw_p_cluster_normal = p)]
  ans
}

# Optional direct comparisons of two correlations. OG draws are shared between groups.
bootstrap_contrast <- function(d, group_column, groups, metric_name) {
  d <- eligible_rows(d[get(group_column) %in% groups], metric_name)
  in_first <- d[[group_column]] == groups[1L]
  n_first <- uniqueN(d$OG[in_first])
  n_second <- uniqueN(d$OG[!in_first])
  ans <- data.table(grouping = group_column, group_a = groups[1L], group_b = groups[2L],
  metric = metric_name, n_OG_a = n_first, n_OG_b = n_second,
  delta_rho_raw = NA_real_, raw_CI_low = NA_real_, raw_CI_high = NA_real_,
  raw_p_boot_centered = NA_real_)
  if (min(n_first, n_second) < minOG) return(ans)
  estimate <- function(idx) {
    a <- idx[in_first[idx]]
    b <- idx[!in_first[idx]]
    raw_spearman(d$Ka[a], d[[metric_name]][a]) - raw_spearman(d$Ka[b], d[[metric_name]][b])
  }
  obs <- estimate(seq_len(nrow(d)))
  ans[, delta_rho_raw := obs]
  if (!is.finite(obs)) return(ans)
  blocks <- split(seq_len(nrow(d)), as.character(d$OG))
  n_og <- length(blocks)
  draws <- replicate(contrastN, estimate(unlist(
  blocks[sample.int(n_og, n_og, replace = TRUE)], use.names = FALSE)))
  draws <- draws[is.finite(draws)]
  if (length(draws) < .8 * contrastN) return(ans)
  ci <- quantile(draws, c(.025, .975), names = FALSE)
  ans[, `:=`(raw_CI_low = ci[1L], raw_CI_high = ci[2L],
  raw_p_boot_centered = (1 + sum(abs(draws - obs) >= abs(obs))) / (1 + length(draws)))]
  ans
}

score_at_threshold <- function(threshold) {
  scored <- copy(q)
  scored[, `:=`(A_SMA_sq = NA_real_, A_SMA_signed = NA_real_,
  B_response24 = NA_real_, B_response24_48 = NA_real_)]
  fits <- lapply(lineages, function(lineage_name) {
    d <- bg[lineage == lineage_name & passes(cm_0, threshold) & passes(proxy_0, threshold)]
    fit <- tryCatch(fit_sma(d), error = function(e) {
      stop('SMA background at minexpr=', threshold, ', ', lineage_name,
      ' (n=', nrow(d), '): ', conditionMessage(e), call. = FALSE)
    })
    idx <- which(scored$lineage == lineage_name &
    passes(scored$cm_0, threshold) & passes(scored$proxy_0, threshold))
    residual <- (log2(scored$proxy_0[idx]) - fit$intercept - fit$slope * log2(scored$cm_0[idx])) /
    sqrt(1 + fit$slope ^ 2)
    if (length(idx)) {
      set(scored, i = idx, j = 'A_SMA_signed', value = residual)
      set(scored, i = idx, j = 'A_SMA_sq', value = residual ^ 2)
    }
    data.table(minexpr = threshold, lineage = lineage_name, background_mode = bgmode,
    n = fit$n, n_effective = fit$n_eff, r = fit$r,
    intercept = fit$intercept, slope = fit$slope)
  })
  # Compare signed log2 fold changes, then take the absolute difference.
  # Both species must pass minexpr at the two endpoints of each interval.
  scored[
  passes(cm_0, threshold) & passes(cm_24, threshold) &
  passes(proxy_0, threshold) & passes(proxy_24, threshold),
  B_response24 := abs((log2(cm_24) - log2(cm_0)) -
  (log2(proxy_24) - log2(proxy_0)))
  ]
  scored[
  passes(cm_24, threshold) & passes(cm_48, threshold) &
  passes(proxy_24, threshold) & passes(proxy_48, threshold),
  B_response24_48 := abs((log2(cm_48) - log2(cm_24)) -
  (log2(proxy_48) - log2(proxy_24)))
  ]
  scored[, minexpr := threshold]
  list(pairs = scored, fits = rbindlist(fits))
}

summarize_threshold <- function(scored, threshold) {
  summarize <- function(d, grouping_name, module_name, stratum_name, family) {
    rbindlist(lapply(metrics, function(metric_name) {
      row <- cluster_estimate(d, metric_name)
      row[, `:=`(minexpr = threshold, background_mode = bgmode, grouping = grouping_name,
      module = module_name, stratum = stratum_name, metric = metric_name,
      BH_family = family)]
      row
    }))
  }
  rows <- list(summarize(scored, 'overall', NA_character_, 'Pooled', 'overall'))
  for (l in lineages) {
    rows[[length(rows) + 1L]] <- summarize(scored[lineage == l],
    'lineage', NA_character_, l, paste0('lineage_', l))
  }
  for (m in modules) {
    d <- scored[Analysis_module == m]
    rows[[length(rows) + 1L]] <- summarize(d, 'module_pooled', m, 'Pooled', 'module_pooled')
    if (module_by_lineage) for (l in lineages) {
      rows[[length(rows) + 1L]] <- summarize(d[lineage == l],
      'module_by_lineage', m, l, paste0('module_', l))
    }
  }
  ans <- rbindlist(rows)
  ans[, raw_q_BH := p.adjust(raw_p_cluster_normal, method = 'BH'), by = BH_family]
  ans[, significant := is.finite(raw_q_BH) & raw_q_BH < .05]
  setcolorder(ans, c('minexpr', 'background_mode', 'grouping', 'module', 'stratum', 'metric',
  'n_pairs', 'n_OG', 'rho_raw', 'raw_CI_low', 'raw_CI_high', 'raw_boot_SE',
  'n_boot_valid', 'raw_p_cluster_normal', 'raw_q_BH', 'BH_family', 'significant'))
  ans
}

# Correlation inference is run once per threshold/stratum/metric.
fits_all <- stats_all <- vector('list', length(minexpr_grid))
for (i in seq_along(minexpr_grid)) {
  threshold <- minexpr_grid[i]
  # A threshold gets the same starting RNG state in single- and multi-threshold runs.
  set.seed(20261027)
  message('Analyzing minexpr = ', threshold)
  result <- score_at_threshold(threshold)
  fits_all[[i]] <- result$fits
  stats_all[[i]] <- summarize_threshold(result$pairs, threshold)
  save(result$pairs, paste0('pair_scores_minexpr_', format(threshold, scientific = FALSE, trim = TRUE)))
}
save(rbindlist(fits_all), 'SMA_0h_model_parameters')

# Rate distributions describe the full set of exact MA immune pairs.
if (!'Ks' %in% names(q)) q[, Ks := NA_real_]
q[, Ks := suppressWarnings(as.numeric(Ks))]
omega_col <- optname(q, c('Ka.Ks', 'Ka/Ks', 'omega', 'Omega', 'Ka_Ks', 'KaKs'))
q[, omega := if (is.na(omega_col)) fifelse(is.finite(Ks) & Ks>0, Ka / Ks, NA_real_)
  else suppressWarnings(as.numeric(get(omega_col)))]
save(q[, .(pair_id, OG, Analysis_module, Subgenome, Ka, Ks, omega)], 'MA_rates_exact_immune_pairs')
save(rbindlist(stats_all), 'final_raw_bootstrap_statistics')

# Optional contrasts run after primary inference to preserve its RNG stream.
if (contrastN > 0L) {
  comparisons <- lapply(minexpr_grid, function(threshold) {
    set.seed(20261028)
    pairs <- fread(file.path(out,
    paste0('pair_scores_minexpr_', format(threshold, scientific = FALSE, trim = TRUE), '.csv')))
    module_pairs <- if (length(modules)>1L) combn(modules, 2L, simplify = FALSE) else list()
    ans <- rbindlist(lapply(metrics, function(metric_name) {
      rbindlist(c(list(bootstrap_contrast(pairs, 'lineage', lineages, metric_name)),
      lapply(module_pairs, function(ab) bootstrap_contrast(pairs, 'Analysis_module', ab, metric_name))))
    }))
    ans[, minexpr := threshold]
    ans[, raw_q_BH := p.adjust(raw_p_boot_centered, 'BH'), by = grouping]
    ans
  })
  save(rbindlist(comparisons), 'raw_bootstrap_correlation_contrasts')
}
