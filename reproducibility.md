# Reproducibility

## Required inputs

Only `0.salmon/sample_info.tsv` and `6.WGCNA/role_space_manual_labels.tsv` are supplied among the metadata inputs below. Other required inputs and large datasets must be supplied separately.

| Location relative to project root | Required content |
|---|---|
| `0.salmon/Fastq/` | Paired `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz` files |
| `ref/{Cm,Ls,Ta}.mrna`, `.genome`, `.gff3` | Matching transcript, genome and annotation versions; transcript IDs must match Salmon |
| `0.salmon/sample_info.tsv` | `sample, species, quant_file, timepoint, replicate`; 25 retained libraries; quantification paths relative to the root |
| `1.normalization/gene_correction_with_edited_cds_length.tsv` | `species, keep_gene, merge_genes, edited_cds_length, status`; semicolon-separated merge IDs, positive CDS lengths, status `OK` |
| `1.normalization/All_pairs_syngroup_w_notes_gene_corrected.csv` | Curated pair map including `comp_id, complexity, Type, Cm1, Cm2, Lsi, Tan` |
| `1.normalization/All_single_syngroup_w_notes_gene_corrected.csv` | Curated single map including `comp_id, complexity, Cm, Lsi, Tan` |
| `OG_info/OG_info_260906_OK.xlsx` | First sheet with `Cmo, GeneName, GeneName_base, Module, subgenome, Type, parent`; retain curated assignments |
| `eggnog_Cm_Ls_Ta/Cmos_eggnog_anno.emapper.annotations` | Cm annotation with `#query, GOs, KEGG_Pathway, KEGG_ko` |
| `eggnog_Cm_Ls_Ta/Lsic_eggnog_anno.emapper.annotations` | Ls annotation with `#query, GOs` |
| `eggnog_Cm_Ls_Ta/Tang_eggnog_anno.emapper.annotations` | Ta annotation with `#query, GOs` |
| `5.KaKs/` or `5.KaKs/data/` | `KaKs_table_revision1.xlsx`, `OG_table_revision1_1025.xlsx`, a curated single-syngroup background table, `Cm_scaledTPM_counts.tsv`, `Ls_scaledTPM_counts.tsv`, `Ta_scaledTPM_counts.tsv` (before TMM scaling), and `all_candidate_sample_factors.tsv`; the script applies the TMM factors |
| `6.WGCNA/role_space_manual_labels.tsv` | Supplied labels and initial label displacements for the final-model plots |

Trim Galore writes to `0.salmon/cleanFastq/<species>/`; Salmon builds indexes in `0.salmon/index_decoy/<species>/` and writes estimates to `0.salmon/quant_decoy/<species>/<sample>/`. The whole genome supplies decoys, and duplicate transcripts are retained. The tx2gene helper exports `ref/<species>_tx2gene.tsv`; normalization reads the GFF3 files directly.

## Outputs and method settings

| Step | Main output location / fixed settings |
|---|---|
| Normalization | `1.normalization/results/01_tximport/`; corrected original counts with tximport length offsets for DESeq2; scaledTPM exported separately |
| Cross-species scaling | `1.normalization/results/02_cross_species_scaling/`; joint sample factors, primary `all_BS + TMM`, sensitivity `median_ratio` and `strict_1to1` |
| Clustering / GO | `2.global_expression/Fig2A_GO/`; k-means `k=6`, seed `1`, original expression/MAD filters and GO backgrounds |
| GSEA / DEG summaries | `2.global_expression/Fig2B_GSEA/`; Wald ranks from offset-aware DESeq2; existing gene-set filters and BH thresholds |
| SC-index / bias | `3.SC_index/tables/`, figure directories and `FigS6_TableS3S6/`; within-species cHEB and the original paired bias models |
| Additivity / PRED | `4.Additivity/results/Corrected_order_NEW_all_BS_TMM/`; parent normalization before MPV construction; explicit inverse-multiplier size factors |
| KaKs associations | `5.KaKs/results/` and `plots/`; existing MA estimates, exact OG/pair matching, whole-OG bootstrap |
| WGCNA | `6.WGCNA/stage1/`, `stage2_only_0.5_22_2_50/`, `stage2_GO_kegg/`; signed bicor, seed `54321`, original screening grid and final model `0.5_22_2_50` |

WGCNA includes network screening, final-model analysis and GO/KEGG enrichment and plots. Enrichment uses the annotated network background. Retain the KEGG name cache when reproducing pathway labels.

## R dependencies

- Analysis: `tximport, DESeq2, rtracklayer, S4Vectors, edgeR, apeglm, WGCNA, fgsea, clusterProfiler, AnnotationDbi, GO.db, KEGGREST`.
- Data handling: `data.table, readxl, openxlsx, tidyverse, dplyr, tidyr, tibble, purrr, readr, stringr, forcats, reshape2`.
- Graphics: `ggplot2, ggh4x, patchwork, cowplot, ggrepel, scales, RColorBrewer, ComplexHeatmap, circlize, gtable`; base R graphics/grid are also used.

Recorded core versions include R 4.5.3, DESeq2 1.50.2, tximport 1.38.2, edgeR 4.8.2, WGCNA 1.74 and apeglm 1.32.0. See the supplied `sessionInfo.txt` files for other package versions.

Complete input datasets and baseline numerical outputs are not bundled; this code-only distribution has not been verified end to end.
