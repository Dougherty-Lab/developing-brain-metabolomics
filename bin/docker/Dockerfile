FROM bioconductor/bioconductor_docker:RELEASE_3_18
 
 ##----------------------------------------------------------------------------##
 ## DOWNLOAD SYSTEM PACKAGES
 ##----------------------------------------------------------------------------##
 
 RUN apt-get update --yes && \
   apt-get upgrade --yes
 
 RUN apt-get install --yes \
   python3-pip \
   python3-dev \
   python3-pip \
   python3-tk
 
 RUN apt-get install --yes \
   libcurl4-openssl-dev \
   libssl-dev \
   libcairo2-dev libxt-dev \
   libxml2-dev \
   libudunits2-dev \
   libhdf5-dev \
   libv8-dev \
   libgdal-dev \
   xorg libx11-dev libglu1-mesa-dev \
   samtools libboost-all-dev \
   libgsl-dev
 
 RUN apt-get install -y \
 	wget \
 	unzip \
 	git
 
 ##----------------------------------------------------------------------------##
 ##  DOWNLOAD R PACKAGES
 ##----------------------------------------------------------------------------##
 
 RUN R -e 'install.packages(c("tidyverse", "remotes",  "BiocManager", "devtools", "knitr", "lintr", "styler", "here", "markdown", "rmarkdown"))'
 RUN R -e 'install.packages(c("Seurat", "scCustomize", "qs", "viridis", "ggpubr", "clustree", "randomcoloR", "gt", "anndata", "SoupX", "pheatmap", "hdf5r"))'
 RUN R -e 'BiocManager::install(c("ComplexHeatmap", "DropletUtils", "AnnotationDbi", "clusterProfiler", "ReactomePA", "apeglm", "DESeq2", "AnnotationDbi", "org.Mm.eg.db"))'
 RUN R -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
 RUN R -e 'devtools::install_github("immunogenomics/harmony")'
 RUN R -e 'devtools::install_github("immunogenomics/presto")'
 RUN R -e 'install.packages(c("grr", "ashr"))'
 RUN R -e 'install.packages("Signac")'
 RUN R -e 'install.packages("svglite")'
 RUN R -e 'install.packages("ggfortify")'
 RUN R -e 'BiocManager::install(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs"))'
 RUN R -e 'devtools::install_github("xia-lab/MetaboAnalystR")'
 RUN R -e 'devtools::install_github("aphalo/ggpp")'
 RUN R -e 'install.packages(c("kableExtra", "ggfortify", "ggrepel", "svglite", "pracma", "cowplot", "UpSetR", "patchwork", "ggtext", "glue", "Matrix", "scales", "harmony", "msigdbr", "arrow"))'
 RUN R -e 'BiocManager::install(c("Biobase", "limma", "ropls", "fgsea", "KEGGREST", "edgeR", "biomaRt", "MOFA2"))'
