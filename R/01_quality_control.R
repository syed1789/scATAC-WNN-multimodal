####################################################################
# 01_quality_control.R
# Multimodal Single-Cell Epigenomics Pipeline
# Author: Syed Nurul Hasan, Ph.D.
# NYU Grossman School of Medicine
#
# Purpose: Quality control and filtering for 10x Multiome data
#          (simultaneous scRNA-seq + scATAC-seq from same cells)
#
# Input:  10x Multiome PBMC dataset (filtered_feature_bc_matrix/
#         and atac_fragments.tsv.gz)
# Output: Filtered Seurat object saved as RDS
####################################################################

# ── Libraries ──────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ── 1. Load Data ───────────────────────────────────────────────────
message("Loading 10x Multiome data...")

# Load RNA counts
# Dataset prefix; default is the 10x PBMC unsorted 3k Multiome set (laptop-sized).
# Use MULTIOME_PREFIX=pbmc_granulocyte_sorted_10k for the 10k set.
prefix <- Sys.getenv("MULTIOME_PREFIX", "pbmc_unsorted_3k")
counts <- Read10X_h5(sprintf("data/%s_filtered_feature_bc_matrix.h5", prefix))

# Extract RNA and ATAC matrices
rna_counts  <- counts$`Gene Expression`
atac_counts <- counts$Peaks

# Load fragment file for ATAC
fragpath <- sprintf("data/%s_atac_fragments.tsv.gz", prefix)

# Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))
genome(annotation)    <- "hg38"

# ── 2. Create Seurat Object ────────────────────────────────────────
message("Creating Seurat object...")

# Create object with RNA assay
pbmc <- CreateSeuratObject(
  counts    = rna_counts,
  assay     = "RNA",
  project   = "PBMC_Multiome"
)

# Add ATAC assay using ChromatinAssay
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts      = atac_counts,
  sep         = c(":", "-"),
  fragments   = fragpath,
  annotation  = annotation
)

message(paste("Total cells before QC:", ncol(pbmc)))

# ── 3. RNA Quality Control ─────────────────────────────────────────
message("Computing RNA QC metrics...")

DefaultAssay(pbmc) <- "RNA"

# Mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Ribosomal gene percentage
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# ── 4. ATAC Quality Control ────────────────────────────────────────
message("Computing ATAC QC metrics...")

DefaultAssay(pbmc) <- "ATAC"

# Nucleosome signal: ratio of mono-nucleosomal to sub-nucleosomal fragments
# Low signal = good open chromatin enrichment
pbmc <- NucleosomeSignal(pbmc)

# TSS enrichment: signal at transcription start sites vs. background
# High TSS = good ATAC library quality
pbmc <- TSSEnrichment(pbmc, fast = FALSE)

# Blacklist ratio: fraction of reads in ENCODE blacklist regions
pbmc[["blacklist_ratio"]] <- FractionCountsInRegion(
  object  = pbmc,
  assay   = "ATAC",
  regions = blacklist_hg38_unified
)

# Peak region fragments per cell
pbmc[["peak_region_fragments"]] <- colSums(
  GetAssayData(pbmc, assay = "ATAC", layer = "counts")
)

# ── 5. Visualize QC Metrics ────────────────────────────────────────
message("Plotting QC metrics...")

# RNA QC violin plots
p_rna_qc <- VlnPlot(
  object   = pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol     = 3,
  pt.size  = 0
) &
  theme(axis.title.x = element_blank())

# ATAC QC violin plots
p_atac_qc <- VlnPlot(
  object   = pbmc,
  features = c("nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol     = 3,
  pt.size  = 0
) &
  theme(axis.title.x = element_blank())

# Grouping variables used by the TSS and fragment-length plots
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 3, "High", "Low")
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 2, "NS > 2", "NS < 2")

# TSS enrichment plot — key ATAC QC metric
p_tss <- TSSPlot(pbmc, group.by = "high.tss") +
  ggtitle("TSS Enrichment") +
  theme_classic()

# Nucleosome signal plot — shows nucleosomal banding pattern
p_ns <- FragmentHistogram(pbmc, group.by = "nucleosome_group") +
  ggtitle("Fragment Length Distribution") +
  theme_classic()

# Save QC plots
ggsave("figures/01a_rna_qc_metrics.pdf",  p_rna_qc,  width = 12, height = 5)
ggsave("figures/01b_atac_qc_metrics.pdf", p_atac_qc, width = 12, height = 5)
ggsave("figures/01c_tss_enrichment.pdf",  p_tss,     width = 6,  height = 5)
ggsave("figures/01d_nucleosome_signal.pdf", p_ns,    width = 6,  height = 5)

# ── 6. Apply QC Filters ────────────────────────────────────────────
message("Applying QC filters...")

# RNA filters
rna_min_features   <- 200     # Minimum genes detected
rna_max_features   <- 5000    # Maximum genes (doublet proxy)
rna_max_mt         <- 20      # Maximum mitochondrial %

# ATAC filters
atac_min_peaks     <- 1000    # Minimum peaks detected
atac_max_peaks     <- 30000   # Maximum peaks (doublet proxy)
atac_min_tss       <- 2       # Minimum TSS enrichment score
atac_max_ns        <- 4       # Maximum nucleosome signal
atac_max_blacklist <- 0.05    # Maximum blacklist fraction

# Peaks detected per cell, as plain metadata (the ATAC assay has no normalized
# "data" layer yet, so Seurat cannot derive these on the fly for subset()).
atac_mat <- GetAssayData(pbmc, assay = "ATAC", layer = "counts")
pbmc$nFeature_ATAC <- Matrix::colSums(atac_mat > 0)
pbmc$nCount_ATAC   <- Matrix::colSums(atac_mat)

md <- pbmc@meta.data
keep <- md$nFeature_RNA      > rna_min_features   &
        md$nFeature_RNA      < rna_max_features   &
        md$percent.mt        < rna_max_mt         &
        md$nFeature_ATAC     > atac_min_peaks     &
        md$nFeature_ATAC     < atac_max_peaks     &
        md$TSS.enrichment    > atac_min_tss       &
        md$nucleosome_signal < atac_max_ns        &
        md$blacklist_ratio   < atac_max_blacklist
pbmc_filtered <- pbmc[, rownames(md)[keep & !is.na(keep)]]

message(paste("Cells after QC filtering:", ncol(pbmc_filtered)))
message(paste("Cells removed:", ncol(pbmc) - ncol(pbmc_filtered)))

# ── 7. QC Summary Table ────────────────────────────────────────────
qc_summary <- data.frame(
  Metric    = c("Total cells (pre-QC)",
                "Cells passing QC",
                "Cells removed",
                "Median genes/cell",
                "Median ATAC peaks/cell",
                "Median TSS enrichment",
                "Median nucleosome signal"),
  Value     = c(ncol(pbmc),
                ncol(pbmc_filtered),
                ncol(pbmc) - ncol(pbmc_filtered),
                round(median(pbmc_filtered$nFeature_RNA)),
                round(median(pbmc_filtered$nFeature_ATAC)),
                round(median(pbmc_filtered$TSS.enrichment), 2),
                round(median(pbmc_filtered$nucleosome_signal), 2))
)

write.csv(qc_summary, "figures/01_qc_summary.csv", row.names = FALSE)
print(qc_summary)

# ── 8. Save Filtered Object ────────────────────────────────────────
message("Saving filtered Seurat object...")
saveRDS(pbmc_filtered, "data/pbmc_multiome_qc_filtered.rds")
message("QC complete. Proceed to 02_rna_processing.R")
