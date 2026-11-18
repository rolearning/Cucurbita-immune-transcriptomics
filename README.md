# Cucurbita immune transcriptomics

## Workflow

| Directory | Inputs and analysis |
|---|---|
| `0.salmon/` | Paired reads, decoy-aware Salmon indexing/quantification and mapping QC |
| `1.normalization/` | Corrected tximport counts, DESeq2 with length offsets, scaledTPM and cross-species scaling |
| `2.global_expression/` | Syngroup clustering, GO enrichment, GSEA and DEG summaries |
| `3.SC_index/` | Genome-wide/copy-level bias, cHEB and response amplitude |
| `4.Additivity/` | Copy-weighted MPV comparisons and parental-relative expression deviation |
| `5.KaKs/` | Coding–expression divergence association analyses |
| `6.WGCNA/` | Network screening, final-model analysis, GO/KEGG enrichment and plots |

## Requirements

- Bash, Trim Galore, Cutadapt, FastQC and Salmon for read processing.
- R and the packages listed in [reproducibility.md](reproducibility.md).
- Reference sequences, curated mapping tables and annotations listed in [reproducibility.md](reproducibility.md). These inputs are not bundled here.

## Run

Prepare inputs as described in [reproducibility.md](reproducibility.md), then run [run_analyses.txt](run_analyses.txt) from the repository root. It sets `ANTI_PROJECT_ROOT` and runs the modules in order.

## Reproducibility

Downstream analyses use 25 retained RNA-seq libraries, excluding `Cm_48h_D2` and `Ls_24h_D1`. The primary cross-species branch is `all_BS + TMM`; `median_ratio` and `strict_1to1` are sensitivity branches. Salmon quantification uses `-l ISR`.

## Data

Raw RNA-seq data are associated with BioProject **PRJNA1556260**. Large sequencing files, reference genomes, annotations and result matrices are not included. Required inputs are listed in [reproducibility.md](reproducibility.md).

If using this code, please cite the associated publication when available.
