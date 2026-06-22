####################################################################
# 04_wnn_integration.R
# Multimodal Single-Cell Epigenomics Pipeline
# Author: Syed Nurul Hasan, Ph.D.
# NYU Grossman School of Medicine
#
# Purpose: Weighted Nearest Neighbor (WNN) integration of
#          scRNA-seq and scATAC-seq modalities
#
# WNN Key Concept:
#   For each cell, WNN learns optimal weights for RNA vs ATAC
#   modalities based on local neighborhood structure.
#   Cells where ATAC is more informative → higher ATAC weight
#   Cells where RNA is more informative → higher RNA weight
#   The joint embedding provides higher resolution than either alone.
#
# This is especially powerful for resolving transitional immune
# states — e.g., classical→non-classical monocyte transitions
# where chromatin changes precede transcriptional changes
####################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

set.seed(42)

# ── Load ATAC-processed object ─────────────────────────────────────
message("Loading ATAC-processed object...")
pbmc <- readRDS("data/pbmc_multiome_atac_processed.rds")

# ── 1. WNN Graph Construction ─────────────────────────────────────
# FindMultiModalNeighbors computes:
# 1. Within-modality kNN graphs (RNA and ATAC separately)
# 2. Cross-modality cell weights based on neighborhood overlap
# 3. Weighted combined graph for downstream clustering and UMAP

message("Building WNN graph...")

pbmc <- FindMultiModalNeighbors(
  object             = pbmc,
  reduction.list     = list("pca", "lsi"),
  dims.list          = list(1:30, 2:30),  # LSI excludes dim 1 (depth)
  modality.weight.name = "RNA.weight",
  verbose            = TRUE
)

# ── 2. WNN UMAP ────────────────────────────────────────────────────
message("Computing WNN UMAP...")

pbmc <- RunUMAP(
  object         = pbmc,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key  = "wnnUMAP_",
  verbose        = FALSE
)

# ── 3. WNN Clustering ──────────────────────────────────────────────
message("Running WNN clustering...")

pbmc <- FindClusters(
  object     = pbmc,
  graph.name = "wsnn",
  algorithm  = 3,        # SLM algorithm
  resolution = 0.8,
  verbose    = FALSE
)

pbmc$wnn_clusters <- pbmc$seurat_clusters

# ── 4. Compare RNA, ATAC, and WNN UMAPs ───────────────────────────
message("Generating comparison UMAPs...")

# Side-by-side comparison of three embeddings
p_rna <- DimPlot(
  pbmc, reduction = "umap.rna",
  group.by = "wnn_clusters",
  label = TRUE, repel = TRUE, pt.size = 0.5
) +
  ggtitle("RNA UMAP") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p_atac <- DimPlot(
  pbmc, reduction = "umap.atac",
  group.by = "wnn_clusters",
  label = TRUE, repel = TRUE, pt.size = 0.5
) +
  ggtitle("ATAC UMAP") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p_wnn <- DimPlot(
  pbmc, reduction = "wnn.umap",
  group.by = "wnn_clusters",
  label = TRUE, repel = TRUE, pt.size = 0.5
) +
  ggtitle("WNN UMAP (Joint RNA + ATAC)") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p_comparison <- p_rna | p_atac | p_wnn

ggsave("figures/04a_umap_comparison_rna_atac_wnn.pdf",
       p_comparison, width = 18, height = 6)

# ── 5. Modality Weight Visualization ──────────────────────────────
# The RNA.weight per cell shows where each modality is more informative
# High RNA weight → gene expression drives clustering
# Low RNA weight (high ATAC weight) → chromatin drives clustering
# Monocyte subtypes often show HIGH ATAC weight (chromatin-driven)

p_weights <- FeaturePlot(
  object    = pbmc,
  features  = "RNA.weight",
  reduction = "wnn.umap",
  pt.size   = 0.3
) +
  scale_color_viridis_c(
    name   = "RNA weight",
    option = "plasma"
  ) +
  ggtitle("WNN Modality Weights\n(Low = ATAC-driven, High = RNA-driven)") +
  theme_classic()

ggsave("figures/04b_wnn_modality_weights.pdf",
       p_weights, width = 7, height = 6)

# ── 6. Monocyte Subtype Resolution ────────────────────────────────
# Key validation: WNN should resolve classical vs non-classical
# monocytes better than RNA or ATAC alone

# Classical monocyte markers
p_cd14 <- FeaturePlot(
  pbmc, features = "CD14",
  reduction = "wnn.umap", pt.size = 0.3
) + scale_color_viridis_c() + theme_classic()

# Non-classical monocyte markers
p_fcgr3a <- FeaturePlot(
  pbmc, features = "FCGR3A",
  reduction = "wnn.umap", pt.size = 0.3
) + scale_color_viridis_c() + theme_classic()

p_cx3cr1 <- FeaturePlot(
  pbmc, features = "CX3CR1",
  reduction = "wnn.umap", pt.size = 0.3
) + scale_color_viridis_c() + theme_classic()

# PECAM1 — key B-ALL associated NCM marker from Witkowski 2020
p_pecam1 <- FeaturePlot(
  pbmc, features = "PECAM1",
  reduction = "wnn.umap", pt.size = 0.3
) + scale_color_viridis_c() + theme_classic()

p_mono_wnn <- (p_cd14 | p_fcgr3a) / (p_cx3cr1 | p_pecam1) +
  plot_annotation(
    title = "Monocyte Subtype Markers on WNN UMAP",
    subtitle = "CD14+ = Classical | FCGR3A/CX3CR1/PECAM1+ = Non-classical"
  )

ggsave("figures/04c_monocyte_markers_wnn.pdf",
       p_mono_wnn, width = 12, height = 10)

# ── 7. Quantify WNN Improvement over Single Modality ──────────────
# Compare silhouette scores between RNA, ATAC, and WNN embeddings
# Higher silhouette = better cluster separation

library(cluster)

compute_silhouette <- function(object, reduction, clusters) {
  coords  <- Embeddings(object, reduction = reduction)[, 1:2]
  dist_m  <- dist(coords)
  sil     <- silhouette(as.numeric(clusters), dist_m)
  mean(sil[, 3])
}

message("Computing silhouette scores for each embedding...")

sil_rna  <- compute_silhouette(pbmc, "umap.rna",  pbmc$wnn_clusters)
sil_atac <- compute_silhouette(pbmc, "umap.atac", pbmc$wnn_clusters)
sil_wnn  <- compute_silhouette(pbmc, "wnn.umap",  pbmc$wnn_clusters)

silhouette_comparison <- data.frame(
  Embedding  = c("RNA only", "ATAC only", "WNN (joint)"),
  Silhouette = round(c(sil_rna, sil_atac, sil_wnn), 4)
)

message("\nSilhouette Score Comparison:")
print(silhouette_comparison)
write.csv(silhouette_comparison,
          "figures/04_silhouette_comparison.csv",
          row.names = FALSE)

p_silhouette <- ggplot(
  silhouette_comparison,
  aes(x = reorder(Embedding, Silhouette),
      y = Silhouette,
      fill = Embedding)
) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("RNA only" = "#4393C3",
                               "ATAC only" = "#D6604D",
                               "WNN (joint)" = "#1A7837")) +
  coord_flip() +
  labs(title = "WNN Integration Improves Cluster Separation",
       subtitle = "Silhouette score comparison across embeddings",
       x = NULL, y = "Mean Silhouette Score") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave("figures/04d_silhouette_comparison.pdf",
       p_silhouette, width = 7, height = 4)

# ── 8. Save WNN-integrated Object ─────────────────────────────────
message("Saving WNN-integrated object...")
saveRDS(pbmc, "data/pbmc_multiome_wnn_integrated.rds")
message("WNN integration complete. Proceed to 05_cell_annotation.R")
