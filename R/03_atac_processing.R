####################################################################
# 03_atac_processing.R
# Multimodal Single-Cell Epigenomics Pipeline
# Author: Syed Nurul Hasan, Ph.D.
# NYU Grossman School of Medicine
#
# Purpose: scATAC-seq processing including LSI dimensionality
#          reduction, chromVAR TF motif deviation scoring,
#          and gene activity matrix computation
#
# Key output: TF motif deviation scores per cell — the core
# analytical layer for identifying regulatory programs driving
# immune cell state transitions in disease
####################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(chromVAR)
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(patchwork)
  library(BiocParallel)
})

set.seed(42)
register(MulticoreParam(4))   # Parallel processing

# ── Load RNA-processed object ──────────────────────────────────────
message("Loading RNA-processed object...")
pbmc <- readRDS("data/pbmc_multiome_rna_processed.rds")
DefaultAssay(pbmc) <- "ATAC"

# ── 1. Peak Calling with MACS2 ────────────────────────────────────
# Call peaks per cell type group rather than globally
# This improves sensitivity for cell-type-specific regulatory elements
message("Calling peaks with MACS2...")

# Use RNA-based clusters as grouping variable for peak calling
# This ensures peaks are called in biologically meaningful groups
peaks <- CallPeaks(
  object     = pbmc,
  group.by   = "rna_clusters",
  macs2.path = "/path/to/macs2"   # Update with your MACS2 path
)

# Remove peaks on non-standard chromosomes and blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(
  x     = peaks,
  ranges = blacklist_hg38_unified,
  invert = TRUE
)

# Quantify fragment counts in new peak set
message("Quantifying ATAC counts in called peaks...")
macs2_counts <- FeatureMatrix(
  fragments  = Fragments(pbmc),
  features   = peaks,
  cells      = colnames(pbmc)
)

# Replace default ATAC assay with MACS2-called peaks
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts      = macs2_counts,
  fragments   = Fragments(pbmc),
  annotation  = Annotation(pbmc[["ATAC"]])
)

# ── 2. LSI Dimensionality Reduction ───────────────────────────────
# LSI (Latent Semantic Indexing): TF-IDF followed by SVD
# Standard dimensionality reduction for scATAC-seq
# LSI dim 1 often correlates with sequencing depth — exclude it
message("Running LSI dimensionality reduction...")

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")  # Use all peaks
pbmc <- RunSVD(pbmc)

# Depth correlation plot — confirms dim 1 is depth-correlated
p_depth_cor <- DepthCor(pbmc) +
  geom_hline(yintercept = c(-0.75, 0.75),
             linetype = "dashed", color = "red") +
  ggtitle("LSI Dimension Correlation with Sequencing Depth") +
  theme_classic()

ggsave("figures/03a_lsi_depth_correlation.pdf", p_depth_cor,
       width = 8, height = 5)

# UMAP from ATAC (LSI dims 2:30, excluding dim 1)
pbmc <- RunUMAP(
  object         = pbmc,
  reduction      = "lsi",
  dims           = 2:30,
  reduction.name = "umap.atac",
  reduction.key  = "atacUMAP_",
  verbose        = FALSE
)

p_atac_umap <- DimPlot(
  object    = pbmc,
  reduction = "umap.atac",
  group.by  = "rna_clusters",   # Color by RNA clusters for comparison
  label     = TRUE,
  repel     = TRUE,
  pt.size   = 0.5
) +
  ggtitle("ATAC-only UMAP (colored by RNA clusters)") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("figures/03b_atac_umap.pdf", p_atac_umap, width = 7, height = 6)

# ── 3. Gene Activity Score Matrix ─────────────────────────────────
# Compute gene activity scores: aggregate ATAC signal in gene body
# and promoter regions as a proxy for gene expression
# Used for cross-modality label transfer and visualization
message("Computing gene activity scores...")

gene_activities <- GeneActivity(pbmc)

pbmc[["Activity"]] <- CreateAssayObject(counts = gene_activities)

pbmc <- NormalizeData(
  object         = pbmc,
  assay          = "Activity",
  normalization.method = "LogNormalize",
  scale.factor   = median(pbmc$nCount_ATAC)
)

# Visualize gene activity at monocyte markers
p_gene_activity <- FeaturePlot(
  object    = pbmc,
  features  = c("CD14", "FCGR3A", "CX3CR1", "S100A8"),
  reduction = "umap.atac",
  ncol      = 2,
  pt.size   = 0.3
) &
  theme_classic() &
  ggtitle("Gene Activity Scores — Monocyte Markers")

ggsave("figures/03c_gene_activity_monocyte.pdf",
       p_gene_activity, width = 10, height = 8)

# ── 4. chromVAR: TF Motif Deviation Scores ────────────────────────
# THE CORE ANALYTICAL LAYER for disease immunology applications
#
# chromVAR computes per-cell deviation scores for TF motifs:
# - Positive deviation = TF binding sites more accessible than expected
# - Negative deviation = TF binding sites less accessible than expected
# - This identifies which TF programs are active in each cell
#
# Application: Identifying TF networks driving NCM reprogramming
# in B-ALL (e.g., NF-κB, IRF, PU.1/SPI1, CEBPB activity per cell)

message("Running chromVAR TF motif deviation analysis...")
message("This is the key step for identifying regulatory programs...")

# Get JASPAR 2020 TF motif database (human)
pfm <- getMatrixSet(
  x      = JASPAR2020,
  opts   = list(collection = "CORE",
                tax_group  = "vertebrates",
                all_versions = FALSE)
)

# Add motif information to ATAC assay
pbmc <- AddMotifs(
  object   = pbmc,
  genome   = BSgenome.Hsapiens.UCSC.hg38,
  pfm      = pfm
)

# Run chromVAR to compute per-cell TF deviation scores
pbmc <- RunChromVAR(
  object  = pbmc,
  genome  = BSgenome.Hsapiens.UCSC.hg38
)

# ── 5. Visualize Key TF Motif Activities ──────────────────────────
# Focus on TFs relevant to monocyte biology and disease immunology

# Key TFs for myeloid/monocyte biology:
# SPI1 (PU.1) — master myeloid TF
# CEBPB — emergency myelopoiesis, inflammatory monocytes
# IRF8 — monocyte vs DC fate decision
# NF-κB (RELA) — inflammatory signaling
# KLF4 — non-classical monocyte identity
# NR4A1 — patrol monocyte program

key_motifs <- c(
  "MA0080.4",   # SPI1/PU.1
  "MA0466.2",   # CEBPB
  "MA0663.1",   # IRF8
  "MA0105.4",   # RELA (NF-κB)
  "MA0039.4",   # KLF4
  "MA0160.1"    # NR4A1
)

DefaultAssay(pbmc) <- "chromvar"

p_tf_umap <- FeaturePlot(
  object    = pbmc,
  features  = key_motifs,
  reduction = "umap.atac",
  ncol      = 3,
  pt.size   = 0.3,
  min.cutoff = "q5",
  max.cutoff = "q95"
) &
  scale_color_viridis_c() &
  theme_classic() &
  theme(axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank())

ggsave("figures/03d_chromvar_tf_motifs_umap.pdf",
       p_tf_umap, width = 12, height = 8)

# TF deviation score violin plots per cluster
# Shows which clusters have highest activity for each TF
p_tf_violin <- VlnPlot(
  object   = pbmc,
  features = key_motifs,
  group.by = "rna_clusters",
  pt.size  = 0,
  ncol     = 3
) &
  theme_classic() &
  theme(axis.title.x = element_blank())

ggsave("figures/03e_chromvar_tf_violin.pdf",
       p_tf_violin, width = 15, height = 10)

# ── 6. Differential Motif Analysis between Monocyte Subtypes ──────
# Find TF motifs differentially enriched in NCMs vs classical monocytes
# This is the analytical approach for the B-ALL immune microenvironment project

DefaultAssay(pbmc) <- "ATAC"

# Identify monocyte clusters from RNA annotation
# (will be refined after WNN integration in script 05)
classical_mono_cells    <- WhichCells(pbmc, idents = c("0", "1"))  # Update with actual cluster IDs
nonclassical_mono_cells <- WhichCells(pbmc, idents = c("7"))       # Update with actual cluster IDs

# Find differentially accessible peaks
da_peaks_mono <- FindMarkers(
  object     = pbmc,
  ident.1    = nonclassical_mono_cells,
  ident.2    = classical_mono_cells,
  only.pos   = TRUE,
  test.use   = "LR",
  latent.vars = "nCount_ATAC",
  min.pct    = 0.05
)

# Find enriched TF motifs in NCM-specific open chromatin
enriched_motifs <- FindMotifs(
  object   = pbmc,
  features = rownames(da_peaks_mono)[da_peaks_mono$avg_log2FC > 0.5]
)

# Save motif enrichment results
write.csv(enriched_motifs,
          "figures/03_ncm_enriched_tf_motifs.csv",
          row.names = FALSE)

message("Top TF motifs enriched in non-classical monocytes:")
print(head(enriched_motifs[order(enriched_motifs$p.adjust), ], 15))

# Motif enrichment plot
p_motif_enrichment <- MotifPlot(
  object  = pbmc,
  motifs  = head(rownames(enriched_motifs), 6),
  assay   = "ATAC"
)

ggsave("figures/03f_ncm_enriched_motifs.pdf",
       p_motif_enrichment, width = 12, height = 8)

# ── 7. Save ATAC-processed Object ─────────────────────────────────
DefaultAssay(pbmc) <- "ATAC"
message("Saving ATAC-processed object...")
saveRDS(pbmc, "data/pbmc_multiome_atac_processed.rds")
message("ATAC processing complete. Proceed to 04_wnn_integration.R")
