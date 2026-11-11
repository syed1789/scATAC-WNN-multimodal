#!/usr/bin/env bash
# Run the whole analysis in containers. Works with Docker (laptop, EC2) or
# Apptainer (HPC, where Docker is not permitted).
#
#   ./run_analysis.sh prefetch   # download input data (needs internet)
#   ./run_analysis.sh run        # run the analysis (default)
#   RUNTIME=apptainer ./run_analysis.sh run
#   DRY_RUN=1 ./run_analysis.sh  # print the commands without running them
set -euo pipefail
cd "$(dirname "$0")"
RUNTIME="${RUNTIME:-docker}"
DRY_RUN="${DRY_RUN:-0}"

run() { if [ "$DRY_RUN" = 1 ]; then echo "+ $*"; else "$@"; fi; }

# build_image <docker-tag> <build-context>
#   docker:    builds the image.
#   apptainer: expects ./<tag>.sif (build once: see HPC_AWS.md); builds it here
#              from the local Docker daemon if Docker happens to be available.
build_image() {
  local tag="$1" ctx="$2" sif="${1%%:*}.sif"
  case "$RUNTIME" in
    docker) run docker build -t "$tag" "$ctx" ;;
    apptainer)
      if [ ! -f "$sif" ]; then
        if command -v docker >/dev/null 2>&1; then
          run docker build -t "$tag" "$ctx"
          run apptainer build "$sif" "docker-daemon://$tag"
        else
          echo "Missing $sif. Build it on a machine with Docker (see HPC_AWS.md)." >&2
          [ "$DRY_RUN" = 1 ] || exit 1
        fi
      fi ;;
    *) echo "RUNTIME must be docker or apptainer" >&2; exit 1 ;;
  esac
}

# ctr <docker-ref-or-tag> <apptainer-source> <command...>
#   apptainer-source: a .sif path, or docker://<ref> for public registry images.
ctr() {
  local ref="$1" src="$2"; shift 2
  case "$RUNTIME" in
    docker) run docker run --rm -e HOME=/tmp -e MPLCONFIGDIR=/tmp -e NUMBA_CACHE_DIR=/tmp \
              -v "$PWD":/work -w /work "$ref" "$@" ;;
    apptainer) run apptainer exec --bind "$PWD":/work --pwd /work "$src" "$@" ;;
  esac
}

step="${1:-run}"

IMG=scatac-wnn:latest
build_image "$IMG" .
case "$step" in
  prefetch)
    echo "Download the 10x PBMC granulocyte-sorted 10k Multiome files into data/ (see README):" >&2
    echo "  pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5" >&2
    echo "  pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz (+ .tbi)" ;;
  run)
    mkdir -p figures data
    for s in R/01_quality_control.R R/02_rna_processing.R R/03_atac_processing.R \
             R/04_wnn_integration.R R/05_cell_annotation.R R/06_tf_motif_analysis.R; do
      ctr "$IMG" scatac-wnn.sif Rscript "$s"
    done ;;
  *) echo "usage: $0 [prefetch|run]" >&2; exit 1 ;;
esac
