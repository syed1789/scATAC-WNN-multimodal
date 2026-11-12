####################################################################
# 06_tf_motif_analysis.R
# Multimodal Single-Cell Epigenomics Pipeline
# Author: Syed Nurul Hasan, Ph.D.
# NYU Grossman School of Medicine
#
# Purpose: Deep TF motif analysis comparing monocyte subtypes
#          Identifies regulatory programs driving NCM identity
#          — directly applicable to B-ALL immune microenvironment
####################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(viridis)
  library(pheatmap)
})

set.seed(42)

pbmc <- readRDS("data/pbmc_multiome_annotated.rds")
Idents(pbmc) <- "cell_type"

# ── 1. Differential TF Motif Activity: NCM vs Classical Mono ──────
# Core analysis: what TF programs distinguish NCM from CM?
# Directly mirrors the question for B-ALL relapse specimens

message("Computing differential TF motif activity...")

DefaultAssay(pbmc) <- "chromvar"

# NCM vs Classical Monocyte
motif_diff_ncm_vs_cm <- FindMarkers(
  object    = pbmc,
  ident.1   = "Non-classical Monocyte",
  ident.2   = "Classical Monocyte",
  only.pos  = FALSE,
  test.use  = "wilcox",
  min.pct   = 0.1
)

motif_diff_ncm_vs_cm$motif_name <- ConvertMotifID(
  object = pbmc,
  id     = rownames(motif_diff_ncm_vs_cm),
  assay  = "ATAC"
)

motif_diff_ncm_vs_cm <- motif_diff_ncm_vs_cm %>%
  arrange(p_val_adj) %>%
  filter(p_val_adj < 0.05)

write.csv(motif_diff_ncm_vs_cm,
          "figures/06_ncm_vs_cm_tf_motifs.csv")

message("\nTop TF motifs enriched in Non-classical Monocytes:")
print(head(motif_diff_ncm_vs_cm %>% filter(avg_log2FC > 0) %>%
             arrange(p_val_adj), 15))

message("\nTop TF motifs depleted in Non-classical Monocytes:")
print(head(motif_diff_ncm_vs_cm %>% filter(avg_log2FC < 0) %>%
             arrange(p_val_adj), 10))

# ── 2. TF Heatmap: Monocyte Subtypes ──────────────────────────────
# Show top differential TF motifs across all three monocyte subtypes

top_ncm_motifs <- motif_diff_ncm_vs_cm %>%
  filter(avg_log2FC > 0) %>%
  head(20) %>%
  rownames()

top_cm_motifs <- motif_diff_ncm_vs_cm %>%
  filter(avg_log2FC < 0) %>%
  head(20) %>%
  rownames()

all_top_motifs <- c(top_ncm_motifs, top_cm_motifs)

# Average chromVAR deviation scores per cell type
mono_types <- intersect(c("Classical Monocyte", "Non-classical Monocyte",
                          "CD14+CD16+ Monocyte"), unique(pbmc$cell_type))

avg_motif_scores <- AverageExpression(
  object   = subset(pbmc, cell_type %in% mono_types),
  features = all_top_motifs,
  assays   = "chromvar",
  group.by = "cell_type"
)$chromvar

# Convert motif IDs to TF names for display
motif_names <- ConvertMotifID(
  object = pbmc,
  id     = rownames(avg_motif_scores),
  assay  = "ATAC"
)
rownames(avg_motif_scores) <- motif_names

# Heatmap
pdf("figures/06a_monocyte_tf_motif_heatmap.pdf", width = 8, height = 12)
pheatmap(
  mat              = avg_motif_scores,
  scale            = "row",
  color            = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
  cluster_rows     = TRUE,
  cluster_cols     = FALSE,
  show_rownames    = TRUE,
  fontsize_row     = 8,
  fontsize_col     = 10,
  main             = "TF Motif Activity: Monocyte Subtypes\n(chromVAR deviation scores)"
)
dev.off()

# ── 3. Key TF Motif UMAP Overlay ──────────────────────────────────
# Show spatial distribution of TF activity on WNN UMAP

# TFs of interest for myeloid/NCM biology
key_tfs <- list(
  "SPI1_MA0080"   = "SPI1 (PU.1) — Master myeloid TF",
  "KLF4_MA0039"   = "KLF4 — NCM identity",
  "NR4A1_MA0160"  = "NR4A1 — Patrol monocyte",
  "CEBPB_MA0466"  = "CEBPB — Inflammatory monocyte",
  "IRF8_MA0663"   = "IRF8 — Monocyte vs DC fate",
  "RELA_MA0105"   = "RELA (NF-κB) — Inflammatory signaling"
)

DefaultAssay(pbmc) <- "chromvar"

# Map TF names to the JASPAR motif IDs present in the chromVAR assay
chromvar_ids   <- rownames(pbmc[["chromvar"]])
chromvar_names <- ConvertMotifID(pbmc, id = chromvar_ids, assay = "ATAC")
key_tfs <- key_tfs[sub("_.*", "", names(key_tfs)) %in% chromvar_names]

tf_plots <- lapply(names(key_tfs), function(tf_id) {
  motif_id <- chromvar_ids[chromvar_names == sub("_.*", "", tf_id)][1]
  FeaturePlot(
    object    = pbmc,
    features  = motif_id,
    reduction = "wnn.umap",
    pt.size   = 0.3,
    min.cutoff = "q10",
    max.cutoff = "q90"
  ) +
    scale_color_viridis_c(option = "plasma") +
    ggtitle(key_tfs[[tf_id]]) +
    theme_classic() +
    theme(plot.title  = element_text(size = 9, face = "bold"),
          axis.title  = element_blank(),
          axis.text   = element_blank(),
          axis.ticks  = element_blank(),
          legend.key.size = unit(0.3, "cm"))
})

p_tf_grid <- wrap_plots(tf_plots, ncol = 3)
ggsave("figures/06b_key_tf_motifs_wnn_umap.pdf",
       p_tf_grid, width = 15, height = 10)

# ── 4. Peak-to-Gene Linkage Analysis ──────────────────────────────
# Links distal regulatory elements (peaks) to target genes
# by correlating chromatin accessibility with gene expression
# across cells — identifies enhancer-gene pairs

message("Computing peak-to-gene linkage...")

DefaultAssay(pbmc) <- "ATAC"

pbmc <- LinkPeaks(
  object            = pbmc,
  peak.assay        = "ATAC",
  expression.assay  = "SCT",
  genes.use         = c("FCGR3A", "CX3CR1", "PECAM1",
                        "CD14", "S100A8", "CDKN1C",
                        "CSF1R", "ITGA4")
)

# Visualize peak-to-gene links at NCM identity loci
p_links_fcgr3a <- CoveragePlot(
  object    = pbmc,
  region    = "FCGR3A",
  features  = "FCGR3A",
  idents    = c("Classical Monocyte", "Non-classical Monocyte"),
  expression.assay = "SCT",
  extend.upstream  = 20000,
  extend.downstream = 5000
) +
  ggtitle("Peak-to-Gene Links at FCGR3A\n(CD16 — NCM Identity Gene)")

p_links_cx3cr1 <- CoveragePlot(
  object    = pbmc,
  region    = "CX3CR1",
  features  = "CX3CR1",
  idents    = c("Classical Monocyte", "Non-classical Monocyte"),
  expression.assay = "SCT",
  extend.upstream  = 20000,
  extend.downstream = 5000
) +
  ggtitle("Peak-to-Gene Links at CX3CR1\n(NCM Patrolling Gene)")

ggsave("figures/06c_peak_gene_links_fcgr3a.pdf",
       p_links_fcgr3a, width = 10, height = 8)
ggsave("figures/06d_peak_gene_links_cx3cr1.pdf",
       p_links_cx3cr1, width = 10, height = 8)

# ── 5. Save Final Object ──────────────────────────────────────────
message("Saving final analyzed object...")
saveRDS(pbmc, "data/pbmc_multiome_final.rds")
message("TF motif analysis complete.")
message("Pipeline complete — see figures/ directory for all outputs.")
