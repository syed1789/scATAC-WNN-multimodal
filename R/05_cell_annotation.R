####################################################################
# 05_cell_annotation.R
# Multimodal Single-Cell Epigenomics Pipeline
# Author: Syed Nurul Hasan, Ph.D.
# NYU Grossman School of Medicine
#
# Purpose: Immune cell type annotation using marker genes,
#          gene activity scores, and TF motif signatures
#          Focus: high-resolution monocyte subtype identification
####################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(viridis)
})

set.seed(42)

# ── Load WNN-integrated object ─────────────────────────────────────
message("Loading WNN-integrated object...")
pbmc <- readRDS("data/pbmc_multiome_wnn_integrated.rds")

# ── 1. Broad Cell Type Annotation ─────────────────────────────────
# Based on canonical marker gene expression and gene activity scores

# Marker gene sets for major immune populations
cell_type_markers <- list(
  "Classical Monocyte"     = c("CD14", "LYZ", "S100A8", "S100A9", "VCAN"),
  "Non-classical Monocyte" = c("FCGR3A", "CX3CR1", "CDKN1C", "PECAM1", "RHOC"),
  "CD14+CD16+ Monocyte"   = c("CD14", "FCGR3A", "CD64", "CX3CR1"),
  "CD1C+ DC"               = c("CD1C", "FCER1A", "CLEC10A", "HLA-DQA1"),
  "pDC"                    = c("CLEC4C", "IL3RA", "LILRA4", "TCF4"),
  "CD8+ T cell"            = c("CD8A", "CD8B", "GZMB", "PRF1"),
  "CD4+ T cell"            = c("CD4", "IL7R", "CCR7", "LEF1"),
  "NK cell"                = c("GNLY", "NKG7", "KLRD1", "NCAM1"),
  "B cell"                 = c("CD19", "MS4A1", "CD79A", "IGHM"),
  "Naive B cell"           = c("MS4A1", "IL4R", "IGHD", "TCL1A"),
  "HSPC"                   = c("CD34", "HOXA9", "AVP", "SPINK2")
)

# Score each cell for each cell type
DefaultAssay(pbmc) <- "SCT"
pbmc <- AddModuleScore(
  object   = pbmc,
  features = cell_type_markers,
  name     = "CellType_Score"
)

# Manual annotation based on marker expression + cluster inspection
# (In practice: iterate between DotPlot examination and annotation)

cluster_annotations <- c(
  "0"  = "Classical Monocyte",
  "1"  = "CD4+ T cell",
  "2"  = "CD8+ T cell",
  "3"  = "NK cell",
  "4"  = "B cell",
  "5"  = "Classical Monocyte",
  "6"  = "CD4+ T cell (Memory)",
  "7"  = "Non-classical Monocyte",
  "8"  = "CD8+ T cell (Effector)",
  "9"  = "HSPC",
  "10" = "CD1C+ DC",
  "11" = "CD4+ T cell (Naive)",
  "12" = "CD14+CD16+ Monocyte",
  "13" = "pDC",
  "14" = "Naive B cell"
)

pbmc$cell_type <- cluster_annotations[as.character(pbmc$wnn_clusters)]

# ── 2. Annotated WNN UMAP ─────────────────────────────────────────
cell_type_colors <- c(
  "Classical Monocyte"     = "#E41A1C",
  "Non-classical Monocyte" = "#FF7F00",
  "CD14+CD16+ Monocyte"   = "#FDBF6F",
  "CD1C+ DC"               = "#A65628",
  "pDC"                    = "#F781BF",
  "CD8+ T cell"            = "#377EB8",
  "CD8+ T cell (Effector)" = "#1F78B4",
  "CD4+ T cell"            = "#4DAF4A",
  "CD4+ T cell (Memory)"  = "#33A02C",
  "CD4+ T cell (Naive)"   = "#B2DF8A",
  "NK cell"                = "#984EA3",
  "B cell"                 = "#6A3D9A",
  "Naive B cell"           = "#CAB2D6",
  "HSPC"                   = "#FFFF33"
)

p_annotated <- DimPlot(
  object     = pbmc,
  reduction  = "wnn.umap",
  group.by   = "cell_type",
  label      = TRUE,
  repel      = TRUE,
  pt.size    = 0.5,
  cols       = cell_type_colors
) +
  ggtitle("WNN UMAP — Annotated Immune Cell Types") +
  theme_classic() +
  theme(legend.text = element_text(size = 9),
        plot.title  = element_text(face = "bold", hjust = 0.5))

ggsave("figures/05a_wnn_annotated_celltypes.pdf",
       p_annotated, width = 10, height = 8)

# ── 3. Monocyte Subtype Resolution — Key Figure ────────────────────
# Highlight the three monocyte populations central to B-ALL biology
# Classical (CD14+), Non-classical (CD16+), Intermediate (CD14+CD16+)

mono_cells <- subset(pbmc,
  cell_type %in% c("Classical Monocyte",
                   "Non-classical Monocyte",
                   "CD14+CD16+ Monocyte"))

p_mono_zoom <- DimPlot(
  object    = mono_cells,
  reduction = "wnn.umap",
  group.by  = "cell_type",
  pt.size   = 1.5,
  cols      = cell_type_colors
) +
  ggtitle("Monocyte Subpopulations — WNN Resolution") +
  theme_classic()

# Marker validation in monocyte subtypes
DefaultAssay(mono_cells) <- "SCT"
p_mono_markers <- DotPlot(
  object   = mono_cells,
  features = c("CD14", "FCGR3A", "CX3CR1", "PECAM1",
               "S100A8", "CDKN1C", "CSF1R", "ITGA4",
               "HLA-DRA", "CD64"),
  group.by = "cell_type"
) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Monocyte Subtype Marker Expression")

p_mono_combined <- p_mono_zoom / p_mono_markers
ggsave("figures/05b_monocyte_subtypes_wnn.pdf",
       p_mono_combined, width = 10, height = 14)

# ── 4. Coverage Plots at Key Disease-Relevant Loci ────────────────
# Visualize chromatin accessibility at monocyte identity loci
# This is a key figure type for the Aifantis lab application

DefaultAssay(pbmc) <- "ATAC"
Idents(pbmc)       <- "cell_type"

# Coverage at FCGR3A locus (CD16 — NCM identity gene)
p_fcgr3a_cov <- CoveragePlot(
  object    = pbmc,
  region    = "FCGR3A",
  idents    = c("Classical Monocyte",
                "Non-classical Monocyte",
                "CD14+CD16+ Monocyte"),
  extend.upstream   = 3000,
  extend.downstream = 3000
) +
  ggtitle("Chromatin Accessibility at FCGR3A (CD16) Locus")

# Coverage at CX3CR1 locus (NCM patrolling gene)
p_cx3cr1_cov <- CoveragePlot(
  object    = pbmc,
  region    = "CX3CR1",
  idents    = c("Classical Monocyte",
                "Non-classical Monocyte",
                "CD14+CD16+ Monocyte"),
  extend.upstream   = 5000,
  extend.downstream = 3000
) +
  ggtitle("Chromatin Accessibility at CX3CR1 Locus")

ggsave("figures/05c_coverage_fcgr3a.pdf",   p_fcgr3a_cov, width = 10, height = 6)
ggsave("figures/05d_coverage_cx3cr1.pdf",   p_cx3cr1_cov, width = 10, height = 6)

# ── 5. Cell Proportion Summary ────────────────────────────────────
cell_proportions <- pbmc@meta.data %>%
  group_by(cell_type) %>%
  summarise(n = n()) %>%
  mutate(proportion = round(n / sum(n) * 100, 2)) %>%
  arrange(desc(proportion))

write.csv(cell_proportions,
          "figures/05_cell_type_proportions.csv",
          row.names = FALSE)

message("\nCell type proportions:")
print(cell_proportions)

# ── 6. Save Annotated Object ──────────────────────────────────────
message("Saving annotated object...")
saveRDS(pbmc, "data/pbmc_multiome_annotated.rds")
message("Annotation complete. Proceed to 06_tf_motif_analysis.R")
