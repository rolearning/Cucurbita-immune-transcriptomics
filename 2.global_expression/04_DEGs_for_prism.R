#!/usr/bin/env Rscript

# Summarize per-species DEG tables from normalization for Prism figure assembly.
# Count nonmissing padj < 0.05 with log2FoldChange >= 1 or <= -1.
# Process the 24 vs 0 h, 48 vs 24 h, and 48 vs 0 h contrasts separately.
# Write DEG_counts_*.tsv to Fig2B_GSEA/plot_prism.

RELEASE_ROOT <- normalizePath(
  Sys.getenv("ANTI_PROJECT_ROOT", unset = "."), mustWork = TRUE
)
root <- RELEASE_ROOT
in_dir <- file.path(root, "1.normalization/results/01_tximport")
out_dir <- file.path(root, "2.global_expression/Fig2B_GSEA/plot_prism")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species <- c("Cm", "Ls", "Ta")
comparisons <- c("24_vs_00", "48_vs_24","48_vs_00")

for (cmp in comparisons) {

  res <- do.call(rbind, lapply(species, function(sp) {

    f <- file.path(in_dir, sp, paste0(sp, "_DEG_", cmp, ".tsv"))
    x <- read.delim(f, check.names = FALSE)

    data.frame(
      Species = sp,
      Down = sum(!is.na(x$padj) & x$padj < 0.05 & x$log2FoldChange <= -1),
      Up   = sum(!is.na(x$padj) & x$padj < 0.05 & x$log2FoldChange >= 1)
    )
  }))

  outfile <- file.path(
    out_dir,
    paste0("DEG_counts_", cmp, ".tsv")
  )

  write.table(
    res,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Written: ", outfile)
}
