FROM bioconductor/bioconductor_docker:RELEASE_3_20
RUN apt-get update && apt-get install -y --no-install-recommends python3-pip \
    && pip3 install --no-cache-dir --break-system-packages macs2 \
    && rm -rf /var/lib/apt/lists/*
RUN R -q -e "install.packages(c('Seurat','Signac','ggplot2','patchwork','viridis','pheatmap','dplyr','cluster'), repos='https://cloud.r-project.org')" \
 && R -q -e "BiocManager::install(c('chromVAR','motifmatchr','JASPAR2020','TFBSTools','BSgenome.Hsapiens.UCSC.hg38','EnsDb.Hsapiens.v86','BiocParallel'), ask=FALSE)"
WORKDIR /work
