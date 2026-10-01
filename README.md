# meo-pIVC

Computational workflows for the analysis of single-cell transcriptomics, untargeted metabolomics, and bulk RNA sequencing during mouse embryogenesis under in vivo and ex utero culture conditions.

This repository brings together preprocessing commands, analysis scripts, and figure-generation code for comparing **meo-pIVC1** and **meo-pIVC2** with in vivo development. The analyses cover cell population composition, developmental progression, RNA velocity, metabolic pathway dynamics, and transcriptional responses to retinoic acid (RA) and 2-deoxy-D-glucose (2-DG).

## Study design and sample labels

The main single-cell comparison uses 42 embryo samples across 14 stage/day groups, with three biological replicates per group. The principal annotated Seurat object contains 373,319 cells.

| Condition | Label in single-cell metadata | Stage/day labels | Samples |
| --- | --- | --- | --- |
| In vivo | `Vivo` | `E9`, `E10`, `E11`, `E12` | 12 |
| meo-pIVC1 | `Vitro1` | `RD1`–`RD4` | 12 |
| meo-pIVC2 | `Vitro2` | `D1`–`D6` | 18 |

Bulk RNA-seq analyses compare Control, RA, and 2-DG conditions in embryo and yolk sac samples. Metabolomics scripts use their own sample-column conventions; these should be preserved when preparing input tables.

## Repository contents

The original filenames are retained. Most R workflows are stored as `.txt` files; the file extension does not indicate the programming language.

| File | Language / format | Purpose |
| --- | --- | --- |
| [`linux code.txt`](linux%20code.txt) | Shell commands | Cell Ranger count commands for individual single-cell libraries. |
| [`moe-aggr.txt`](moe-aggr.txt) | R notebook-style code | Import aggregated 10x matrices, assign sample identities, inspect quality metrics, filter cells, and export per-sample matrices. |
| [`scrublet-python.txt`](scrublet-python.txt) | Python notebook export | Sample-wise doublet detection with Scrublet. |
| [`moe-annotation.txt`](moe-annotation.txt) | R notebook-style code | Cell population annotation, marker visualization, and examination of lineage subsets. |
| [`moe-pseudo-bulk-correlation.txt`](moe-pseudo-bulk-correlation.txt) | R | Replicate-level pseudo-bulk analysis and stage-to-stage transcriptomic Pearson correlations. |
| [`scvelo.py`](scvelo.py) | Python | scVelo RNA velocity analysis and matched culture-group plots for Figure 3d–f. |
| [`bulk_RNAseq_pipeline.txt`](bulk_RNAseq_pipeline.txt) | Bash with embedded Python | Configurable preprocessing demo from paired-end FASTQ files to gene fragment counts and QC summaries. |
| [`Figure1.txt`](Figure1.txt) | R | meo-pIVC1 glycolysis and OxPhos/TCA metabolite heatmaps. |
| [`Figure3.txt`](Figure3.txt) | R and Python | Cell population UMAPs and proportions, pseudo-bulk correlations, RNA velocity, placode and endoderm subtypes, and metabolomics comparisons. |
| [`Figure4.txt`](Figure4.txt) | R | RA-versus-Control KEGG GSEA, signalling pathway scores, and RA, HOX/Neural, TGFB, and NODAL expression heatmaps. |
| [`ExtendedData_Figure2.txt`](ExtendedData_Figure2.txt) | R | meo-pIVC1 metabolomics PCA and pathway dynamics, 2-DG KEGG GSEA, and UPR/ER-stress gene expression heatmaps. |
| [`ExtendedData_Figure6.txt`](ExtendedData_Figure6.txt) | R | Single-cell QC, stage-resolved UMAPs, developmental age modelling, cluster-level correlations, and placode/endoderm marker and composition plots. |
| [`ExtendedData_Figure7.txt`](ExtendedData_Figure7.txt) | R | meo-pIVC2 metabolomics PCA, pathway dynamics, and comparison of supplied meo-pIVC1/meo-pIVC2 pathway summary scores. |

## Software requirements

Install the dependencies required for the workflow you intend to run.

| Analysis | Main dependencies |
| --- | --- |
| Single-cell analysis in R | Seurat, Matrix, edgeR, limma, dplyr, tidyr, tibble, plyr, ggplot2, patchwork, scales, corrplot |
| Metabolomics and heatmaps | readxl, stringr, dplyr, ggplot2, ggforce, scales, ComplexHeatmap, circlize, pheatmap |
| Bulk expression and enrichment | DESeq2, clusterProfiler, org.Mm.eg.db, AnnotationDbi |
| Python doublet detection and velocity | scrublet, scvelo, numpy, pandas, scipy, matplotlib |
| Sequencing preprocessing | Cell Ranger; STAR, Subread/featureCounts, fastp, samtools, Python, and standard command-line utilities |

The original aggregation and annotation excerpts also load pagoda2, velocyto.R, monocle, cowplot, harmony, RColorBrewer, and magrittr. These are needed when executing those excerpts unchanged.

Some figure scripts document R 4.1.3. Single-cell code uses Seurat v4-style assay slots, and the velocity code uses the scVelo API shown in the scripts. Use a compatible environment and retain its package versions with analysis outputs. The bulk preprocessing demo specifies STAR 2.7.11b, featureCounts 2.1.1, Python ≥3.9, and mouse GRCm38 / Ensembl release 89.

## License

This repository is distributed under the GNU General Public License v3.0. See [`LICENSE`](LICENSE) for the full terms.
