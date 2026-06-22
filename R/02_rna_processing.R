####################################################################
# 02_rna_processing.R
# Multimodal Single-Cell Epigenomics Pipeline
# Author: Syed Nurul Hasan, Ph.D.
# NYU Grossman School of Medicine
#
# Purpose: scRNA-seq normalization, dimensionality reduction,
#          and initial clustering
####################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ── Load QC-filtered object ────────────────────────────────────────
message("Loading QC-filtered object...")
pbmc <- readRDS("data/pbmc_multiome_qc_filtered.rds")
DefaultAssay(pbmc) <- "RNA"

# ── 1. Normalization with SCTransform ─────────────────────────────
# SCTransform: variance-stabilizing normalization that regresses
# out sequencing depth, more robust than log-normalization for
# comparing cells with very different library sizes
message("Running SCTransform normalization...")

pbmc <- SCTransform(
  object            = pbmc,
  assay             = "RNA",
  vars.to.regress   = "percent.mt",   # Regress out mitochondrial %
  verbose           = FALSE,
  return.only.var.genes = FALSE
)

# ── 2. PCA Dimensionality Reduction ───────────────────────────────
message("Running PCA...")

pbmc <- RunPCA(
  object   = pbmc,
  assay    = "SCT",
  npcs     = 50,
  verbose  = FALSE
)

# Elbow plot to determine optimal number of PCs
p_elbow <- ElbowPlot(pbmc, ndims = 50) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
  ggtitle("PCA Elbow Plot — RNA") +
  theme_classic()

ggsave("figures/02a_rna_pca_elbow.pdf", p_elbow, width = 7, height = 5)

# ── 3. UMAP and Clustering (RNA only) ─────────────────────────────
# This provides initial RNA-based clustering for comparison
# with the final WNN joint clustering
message("Running RNA-based UMAP and clustering...")

n_dims_rna <- 30   # Selected from elbow plot

pbmc <- RunUMAP(
  object     = pbmc,
  dims       = 1:n_dims_rna,
  reduction  = "pca",
  reduction.name = "umap.rna",
  reduction.key  = "rnaUMAP_",
  verbose    = FALSE
)

pbmc <- FindNeighbors(
  object  = pbmc,
  dims    = 1:n_dims_rna,
  verbose = FALSE
)

pbmc <- FindClusters(
  object     = pbmc,
  resolution = 0.8,
  verbose    = FALSE,
  graph.name = "SCT_snn"
)

# Store RNA clusters
pbmc$rna_clusters <- pbmc$seurat_clusters

# ── 4. RNA UMAP Visualization ─────────────────────────────────────
p_rna_umap <- DimPlot(
  object     = pbmc,
  reduction  = "umap.rna",
  group.by   = "rna_clusters",
  label      = TRUE,
  repel      = TRUE,
  pt.size    = 0.5
) +
  ggtitle("RNA-only UMAP") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("figures/02b_rna_umap_clusters.pdf", p_rna_umap, width = 7, height = 6)

# ── 5. Key Immune Marker Expression ───────────────────────────────
# Monocyte markers — central to disease immunology applications
immune_markers <- c(
  # Classical monocytes (CD14+)
  "CD14", "S100A8", "S100A9", "LYZ",
  # Non-classical monocytes (CD16+) — key B-ALL associated population
  "FCGR3A", "CX3CR1", "CDKN1C", "PECAM1",
  # T cells
  "CD3D", "CD3E", "CD8A", "CD4",
  # NK cells
  "GNLY", "NKG7", "KLRD1",
  # B cells
  "CD19", "MS4A1", "CD79A",
  # Dendritic cells
  "CLEC9A", "CD1C", "FCER1A",
  # HSPCs
  "CD34", "HOXA9"
)

p_markers <- DotPlot(
  object   = pbmc,
  features = immune_markers,
  group.by = "rna_clusters"
) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Immune Cell Marker Expression by RNA Cluster")

ggsave("figures/02c_rna_marker_dotplot.pdf", p_markers, width = 14, height = 10)

# Feature plots for monocyte subtype markers
# (Critical for validating NCM identification)
p_mono_markers <- FeaturePlot(
  object    = pbmc,
  features  = c("CD14", "FCGR3A", "CX3CR1", "PECAM1",
                "S100A8", "CDKN1C"),
  reduction = "umap.rna",
  ncol      = 3,
  pt.size   = 0.3
) &
  theme_classic() &
  theme(axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank())

ggsave("figures/02d_monocyte_markers_featureplot.pdf",
       p_mono_markers, width = 12, height = 8)

# ── 6. Save RNA-processed Object ──────────────────────────────────
message("Saving RNA-processed object...")
saveRDS(pbmc, "data/pbmc_multiome_rna_processed.rds")
message("RNA processing complete. Proceed to 03_atac_processing.R")
