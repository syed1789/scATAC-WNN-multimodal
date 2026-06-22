# Multimodal Single-Cell Epigenomics: scATAC-seq + scRNA-seq WNN Integration

**Author:** Syed Nurul Hasan, Ph.D.  
**Affiliation:** NYU Grossman School of Medicine (Ramasamy Lab)  
**Contact:** nayanchoton1789@gmail.com | [ORCID: 0000-0002-4564-872X](https://orcid.org/0000-0002-4564-872X)

---

## Overview

This repository demonstrates a complete end-to-end pipeline for **multimodal single-cell analysis** integrating scRNA-seq and scATAC-seq data from the same cells using 10x Genomics Multiome technology. The pipeline applies **Weighted Nearest Neighbor (WNN)** integration to jointly model gene expression and chromatin accessibility, enabling higher-resolution immune cell type identification than either modality alone.

The analytical framework here directly mirrors approaches used in the characterization of immune microenvironment dynamics in disease contexts — including leukemia, cardiovascular disease, and autoimmunity — where resolving transitional cell states at chromatin resolution is critical for identifying transcription factor (TF) networks driving immune cell reprogramming.

---

## Scientific Motivation

Single-cell ATAC-seq (scATAC-seq) captures chromatin accessibility genome-wide, providing a window into the **regulatory landscape** of individual cells. When integrated with scRNA-seq via WNN:

- Cells where chromatin is more informative receive higher ATAC weight
- Cells where transcription is more informative receive higher RNA weight
- The joint embedding resolves transitional immune states invisible to either modality alone
- **TF motif deviation scores** (chromVAR) identify which regulatory programs are active per cell

This is particularly powerful for questions like:
> *What transcription factor networks drive non-classical monocyte reprogramming in disease?*  
> *How does chromatin accessibility change in immune cells across disease progression?*  
> *Which regulatory programs are pharmacologically accessible?*

---

## Pipeline Overview

```
10x Multiome Raw Data (PBMC)
        │
        ▼
1. Quality Control & Filtering
   ├── RNA: nFeature, nCount, % mitochondrial
   └── ATAC: TSS enrichment, nucleosome signal, peak counts
        │
        ▼
2. RNA Processing (Seurat)
   ├── Normalization (SCTransform)
   ├── PCA dimensionality reduction
   └── Initial RNA-based clustering
        │
        ▼
3. ATAC Processing (Signac + ArchR)
   ├── Peak calling (MACS2)
   ├── LSI dimensionality reduction
   ├── chromVAR TF motif deviation scores
   └── Gene activity score matrix
        │
        ▼
4. WNN Integration (Seurat)
   ├── Weighted Nearest Neighbor graph
   ├── Joint UMAP embedding
   └── WNN-based clustering
        │
        ▼
5. Biological Interpretation
   ├── Immune cell type annotation
   ├── Monocyte subpopulation resolution
   ├── TF motif enrichment per cluster
   ├── Peak-to-gene linkage analysis
   └── Cell-Cell Communication (CellChat)
        │
        ▼
6. Visualization & Export
   ├── UMAP plots (RNA, ATAC, WNN)
   ├── TF deviation heatmaps
   ├── Coverage plots at key loci
   └── Publication-ready figures
```

---

## Key Features

| Feature | Tool | Purpose |
|---|---|---|
| scRNA-seq processing | Seurat v5 | Normalization, PCA, clustering |
| scATAC-seq processing | Signac | Peak calling, LSI, chromatin QC |
| TF motif analysis | chromVAR | Per-cell TF deviation scores |
| WNN integration | Seurat v5 | Joint RNA+ATAC embedding |
| Gene activity | Signac | Imputed gene expression from ATAC |
| Peak-to-gene links | Signac | Regulatory element mapping |
| Cell communication | CellChat | Ligand-receptor network analysis |
| Visualization | ggplot2, Nebulosa | Publication-ready figures |

---

## Data

This pipeline uses the publicly available **10x Genomics PBMC Multiome** dataset:
- 10,412 PBMCs from a healthy human donor
- Simultaneous scRNA-seq + scATAC-seq from same cells
- Download: [10x Genomics Dataset Portal](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0)

**Why PBMC Multiome for immune biology:**  
PBMCs contain the same immune cell populations relevant to disease contexts — classical monocytes (CD14+), non-classical monocytes (CD16+), T cells, NK cells, B cells, and DCs. This dataset provides a clean baseline to validate multimodal integration approaches before applying them to disease specimens.

---

## Repository Structure

```
scATAC-WNN-multimodal/
├── README.md                          # This file
├── R/
│   ├── 01_quality_control.R           # QC filtering for RNA and ATAC modalities
│   ├── 02_rna_processing.R            # SCTransform, PCA, RNA clustering
│   ├── 03_atac_processing.R           # Peak calling, LSI, chromVAR
│   ├── 04_wnn_integration.R           # WNN graph, joint UMAP, WNN clustering
│   ├── 05_cell_annotation.R           # Immune cell type annotation
│   ├── 06_tf_motif_analysis.R         # TF deviation scores, motif enrichment
│   ├── 07_peak_gene_links.R           # Peak-to-gene linkage analysis
│   └── 08_cellchat_communication.R    # Ligand-receptor network analysis
├── docs/
│   ├── pipeline_guide.md              # Step-by-step execution guide
│   └── biological_interpretation.md   # Interpreting TF motif results
├── figures/                           # Output figures directory
└── data/                              # Data directory (see download instructions)
```

---

## Requirements

```r
# Core packages
install.packages("Seurat")          # v5.0+
install.packages("Signac")          # v1.12+
BiocManager::install("chromVAR")
BiocManager::install("motifmatchr")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("JASPAR2020")
install.packages("CellChat")

# Visualization
install.packages("ggplot2")
install.packages("patchwork")
install.packages("viridis")
BiocManager::install("Nebulosa")

# Peak calling (system)
# conda install -c bioconda macs2
```

---

## Quick Start

```r
# Clone repository
# git clone https://github.com/syed1789/scATAC-WNN-multimodal

# Run pipeline sequentially
source("R/01_quality_control.R")
source("R/02_rna_processing.R")
source("R/03_atac_processing.R")
source("R/04_wnn_integration.R")
source("R/05_cell_annotation.R")
source("R/06_tf_motif_analysis.R")
source("R/07_peak_gene_links.R")
source("R/08_cellchat_communication.R")
```

---

## Biological Context: Application to Disease Immunology

While this pipeline is demonstrated on healthy PBMCs, it was developed to address specific questions in disease immunology:

**Leukemia immune microenvironment:**  
Non-classical monocytes (CD14dim CD16+) expand at B-ALL diagnosis and relapse and predict inferior patient survival (Witkowski et al., Cancer Cell 2020). The chromatin-level determinants of this expansion — which TF networks drive NCM emergence and sustain their pro-leukemic phenotype — remain unmapped. This pipeline provides the analytical framework to address that question using paired scRNA-seq + scATAC-seq on immune fractions from matched diagnosis/remission/relapse specimens.

**Cardiovascular immunology:**  
Myeloid cell state transitions in ischemic and diabetic hearts drive immune-metabolic pathway dysregulation. WNN integration of paired RNA and ATAC data resolves transitional myeloid states invisible to transcriptomics alone, enabling identification of TF networks driving pathological immune remodeling.

**Autoimmunity:**  
Regulatory T cell (Treg) transcriptional programs — including Bcl11b-dependent chromatin accessibility — are disrupted in autoimmune disease. scATAC-seq of Treg populations enables chromatin-level characterization of tolerance versus dysfunction, with direct implications for leukemia-associated Treg immunosuppression.

---

## Key References

1. Stuart T, et al. (2021). Multimodal single-cell chromatin analysis with Signac. *Nature Methods*, 18:1333–1341.
2. Hao Y, et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184:3573–3587.
3. Schep AN, et al. (2017). chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nature Methods*, 14:975–978.
4. Witkowski MT, et al. (2020). Extensive remodeling of the immune microenvironment in B cell acute lymphoblastic leukemia. *Cancer Cell*, 37:867–882.
5. Narang S, et al. (2024). Clonal evolution of the 3D chromatin landscape in patients with relapsed pediatric B-cell acute lymphoblastic leukemia. *Nature Communications*, 15:7425.

---

## Related Work

This pipeline is part of a broader single-cell genomics portfolio developed during postdoctoral research at NYU Grossman School of Medicine:

- **scRNA-seq analysis:** Seurat-based immune cell profiling in cardiac disease
- **Spatial transcriptomics:** VisiumHD analysis of tissue immune architecture  
- **Multi-omics integration:** WNN pipeline for cardiac immune-metabolic interactions
- **Cell-cell communication:** CellChat network modeling of myeloid-stromal interactions

---

*For questions or collaboration inquiries: nayanchoton1789@gmail.com*
