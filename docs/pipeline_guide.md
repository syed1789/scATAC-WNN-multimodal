# Pipeline Execution Guide

## Prerequisites

### Software Requirements
```bash
# R >= 4.3.0
# Bioconductor >= 3.18
# MACS2 (for peak calling)
conda install -c bioconda macs2

# Or via pip
pip install MACS2
```

### R Package Installation
```r
# CRAN packages
install.packages(c("Seurat", "Signac", "ggplot2", "patchwork",
                   "dplyr", "viridis", "pheatmap", "cluster"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "chromVAR",
  "motifmatchr",
  "JASPAR2020",
  "TFBSTools",
  "BSgenome.Hsapiens.UCSC.hg38",
  "EnsDb.Hsapiens.v86",
  "BiocParallel"
))

# CellChat for ligand-receptor analysis
devtools::install_github("sqjin/CellChat")
```

---

## Data Download

Download the 10x Multiome PBMC dataset from 10x Genomics:

```bash
mkdir -p data && cd data

# Filtered feature barcode matrix (RNA + ATAC)
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5

# ATAC fragment file
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
```

---

## Pipeline Execution

Run scripts sequentially from the repository root:

```r
setwd("/path/to/scATAC-WNN-multimodal")

source("R/01_quality_control.R")     # ~10 min
source("R/02_rna_processing.R")      # ~15 min
source("R/03_atac_processing.R")     # ~45 min (chromVAR is slow)
source("R/04_wnn_integration.R")     # ~20 min
source("R/05_cell_annotation.R")     # ~15 min
source("R/06_tf_motif_analysis.R")   # ~20 min
```

**Estimated total runtime:** ~2–2.5 hours on a standard laptop (8GB RAM minimum, 16GB recommended)

---

## Expected Outputs

```
figures/
├── 01a_rna_qc_metrics.pdf
├── 01b_atac_qc_metrics.pdf
├── 01c_tss_enrichment.pdf
├── 01d_nucleosome_signal.pdf
├── 01_qc_summary.csv
├── 02a_rna_pca_elbow.pdf
├── 02b_rna_umap_clusters.pdf
├── 02c_rna_marker_dotplot.pdf
├── 02d_monocyte_markers_featureplot.pdf
├── 03a_lsi_depth_correlation.pdf
├── 03b_atac_umap.pdf
├── 03c_gene_activity_monocyte.pdf
├── 03d_chromvar_tf_motifs_umap.pdf
├── 03e_chromvar_tf_violin.pdf
├── 03f_ncm_enriched_motifs.pdf
├── 03_ncm_enriched_tf_motifs.csv
├── 04a_umap_comparison_rna_atac_wnn.pdf  ← Key figure
├── 04b_wnn_modality_weights.pdf
├── 04c_monocyte_markers_wnn.pdf
├── 04d_silhouette_comparison.pdf
├── 04_silhouette_comparison.csv
├── 05a_wnn_annotated_celltypes.pdf       ← Key figure
├── 05b_monocyte_subtypes_wnn.pdf         ← Key figure
├── 05c_coverage_fcgr3a.pdf
├── 05d_coverage_cx3cr1.pdf
├── 05_cell_type_proportions.csv
├── 06a_monocyte_tf_motif_heatmap.pdf     ← Key figure
├── 06b_key_tf_motifs_wnn_umap.pdf        ← Key figure
├── 06c_peak_gene_links_fcgr3a.pdf
└── 06d_peak_gene_links_cx3cr1.pdf
```

---

## Key Figures Explained

### Figure 04a — RNA vs ATAC vs WNN UMAP Comparison
The central validation figure. Shows that WNN joint embedding resolves immune cell populations better than either modality alone. Non-classical monocytes (CD16+) appear as a distinct cluster in WNN that may be merged with classical monocytes in RNA-only analysis.

### Figure 05b — Monocyte Subtype Resolution
Zoomed view of three monocyte populations: classical (CD14+), non-classical (CD16+/CX3CR1+/PECAM1+), and intermediate (CD14+CD16+). The non-classical population is the leukemia-associated subtype from Witkowski et al. (Cancer Cell, 2020).

### Figure 06a — TF Motif Heatmap
Shows chromVAR deviation scores for top differential TF motifs across monocyte subtypes. TFs enriched in non-classical monocytes (e.g., KLF4, NR4A1) versus classical monocytes (e.g., CEBPB, SPI1) represent candidate regulatory programs driving the classical→non-classical transition.

### Figure 06b — TF Motif Activity on WNN UMAP
Spatial distribution of key myeloid TF programs across the immune landscape. Enables visual identification of which clusters express which regulatory programs.

---

## Troubleshooting

**chromVAR slow or failing:**
```r
# Reduce parallelization
register(MulticoreParam(2))
# Or run on subset first
pbmc_mono <- subset(pbmc, cell_type %in% c("Classical Monocyte", "Non-classical Monocyte"))
```

**MACS2 not found:**
```r
# Specify full path
peaks <- CallPeaks(pbmc, macs2.path = "/usr/local/bin/macs2")
# Or use conda environment path
peaks <- CallPeaks(pbmc, macs2.path = "/opt/conda/bin/macs2")
```

**Memory issues:**
```r
# Process in batches or reduce feature set
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q25")  # Use top 75% of peaks only
```
