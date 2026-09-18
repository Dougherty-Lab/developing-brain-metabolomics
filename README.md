# Metabolite-gene expression interactions in the developing human cortex
## Sneha M. Chaturvedi, Kelsey Hennick, Emma Jones, Rencheng Wang, Jaenyeon Kim, Minsoo Son, Leah P. Shriver, Young Ah Goo, Tomasz Nowakowski, Joseph D. Dougherty

Analysis repo for our Chaturvedi et al. 2026 manuscript examining sample-matched metabolite values with single-nucleus gene expression data.

## Abstract
The molecular environment of the early cortex plays a large role on developmental trajectory, as changes in gene expression and cellular signaling both impact neurodevelopmental outcomes. Small metabolites and nutrients are products of specific genes’ functions (e.g., metabolic enzymes) and directly regulate gene expression (e.g., via nuclear hormone receptors). While several such nutrients are essential for brain development, there is no comprehensive catalog of the metabolites present in the developing brain, the degree and nature of relationships between such small molecules and transcriptional regulation, and whether they are global or vary by cell type. In this study, we combine comprehensive metabolomics with single nucleus (sn)RNAseq to test the interrelationships between gene expression and metabolite abundance in the human developing cortex during mid-gestation, and with cell-type specificity. We detected 525 distinct metabolites, with significant changes across gestational week and subtle sex effects, and provide the first measures of testosterone, progesterone, and estradiol directly from the developing brain. Integration of metabolomics with snRNAseq revealed nearly half of metabolites are significantly associated with expression of 542 genes across eleven cell types, via cell-type specific gene regulation. Such widespread relationships provide a resource for new understanding of metabolic impacts during human cortical development.  


## Repo Contents

```
├── bin/                    # Container build files (Dockerfile, Singularity def)
├── data/                   # Input data and cached intermediates (see data/README.md)
│   └── cache/              # Pseudobulk aggregates and lightweight Seurat subsets
├── doc/                    # Annotation files, gene lists, supplementary tables
│   ├── targeted/           # Hormone targeted concentration data
│   ├── untargeted/         # Metabolomics annotation files
│   └── gene_lists/         # Curated gene sets 
├── results/                # Analysis outputs, organized by analysis type
│   ├── targeted/           # Hormone analysis results
│   ├── untargeted/         # Metabolomics QC and analysis results
│   ├── gene-hormone/       # Hormone–gene association results
│   ├── gene-metabolite/    # Metabolite–gene association results
│   └── mofa/               # MOFA2 multi-omics integration results
└── src/                    # All analysis scripts (see src/README.md)
    ├── targeted/           # Gonadal hormone analyses
    ├── untargeted/         # Metabolomics processing and statistical analysis
    ├── gene-hormone/       # Hormone–gene limma-voom associations
    ├── gene-metabolite/    # Metabolite–gene limma-voom associations
    └── mofa/               # MOFA2 multi-omics factor analysis
```

## System Requirements

### Hardware Requirements

Most scripts run on a standard workstation. The following are exceptions that require high-performance computing:

| Script | RAM | Time | Notes |
|--------|-----|------|-------|
| `filter_recluster.R` | ~200 GB | ~2–4 hrs | Loads full Seurat object; Harmony reclustering |
| `prep_data.R` | ~200 GB | ~30 min | Loads filtered Seurat; caches small derived objects |
| `metabolite_gene_limma.R` | ~32 GB | ~4–8 hrs | Loops metabolites through limma-voom |

All other scripts run with ≤16 GB RAM in under 30 minutes.

### Software Requirements

#### Container (recommended)

A Singularity container encapsulating the full computational environment can be found at https://hub.docker.com/r/emmafjones/developing-brain-metabolomics 

The container was built and tested on:
- **Linux:** Ubuntu 22.04 
- **Singularity:** v3.8.0

#### Manual installation

Without the container, install:
- **R** ≥ 4.3
- **Python** ≥ 3.9 with `mofapy2` (MOFA2 backend only)

Key R packages (exact versions recorded in `sessionInfo()` at the end of each script):

| Category | Packages |
|----------|----------|
| Core | `tidyverse`, `readxl`, `arrow` |
| Single-cell | `Seurat`, `edgeR`, `limma` |
| Metabolomics | `Biobase`, `mixOmics` |
| Enrichment | `fgsea`, `msigdbr`, `KEGGREST` |
| Multi-omics | `MOFA2`, `basilisk`, `reticulate` |
| Visualization | `patchwork`, `cowplot`, `ggpubr`, `svglite`, `UpSetR` |

## Installation Guide

### Using the container


```bash
# Pull the container from Docker Hub (~2.58 GB)
docker pull emmafjones/developing-brain-metabolomics:latest

# Or convert to Singularity for HPC use
singularity pull developing-brain-metabolomics.sif \
  docker://emmafjones/developing-brain-metabolomics:latest

# Clone the repository
git clone https://github.com/[org]/developing-brain-metabolomics.git
cd developing-brain-metabolomics
```

### Without the container

```bash
git clone https://github.com/[org]/developing-brain-metabolomics.git
cd developing-brain-metabolomics

# Install R dependencies (from an R session)
install.packages(c("tidyverse", "readxl", "arrow", "patchwork", "cowplot",
                    "ggpubr", "svglite", "UpSetR", "mixOmics", "kableExtra"))
BiocManager::install(c("Seurat", "edgeR", "limma", "fgsea", "msigdbr",
                        "KEGGREST", "MOFA2", "basilisk"))
```

## License

This project is covered under the [MIT License](https://github.com/Dougherty-Lab/developing-brain-metabolomics?tab=MIT-1-ov-file#).


