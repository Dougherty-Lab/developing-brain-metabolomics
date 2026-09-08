## filter_recluster.R
##
## Standalone script: full filtering pipeline + reclustering
## Run via SLURM batch job:
## sbatch run_filter_recluster.sh
##
## Outputs:
##   ../../data/SSD_RNAonly_filtered.rds - filtered + reclustered Seurat object

Packages <- c("tidyverse", "Seurat", "harmony")
lapply(Packages, library, character.only = TRUE)
setwd("/scratch/jdlab/emma/developing-brain-metabolomics/src/gene-hormone/")

set.seed(123)

sample_col <- "Sample"

# ---- 1. Load ----------------------------------------------------------------
cat("Loading SSD_RNAonly.rds...\n")
t0 <- Sys.time()
SSD_data <- readRDS("../../data/SSD_RNAonly.rds")
cat("Load time:", format(Sys.time() - t0), "\n")
cat("Cells loaded:", ncol(SSD_data), "\n")

# ---- 2. Metabolomics sample exclusion ---------------------------------------
cat("\nExcluding metabolomics-absent samples...\n")
excluded_samples <- c(
  "SSD14_tissue14", "SSD16_tissue16",
  "SSD28_tissue28", "SSD29_tissue29"
)
keep <- colnames(SSD_data)[!SSD_data@meta.data[[sample_col]] %in% excluded_samples]
SSD_data <- SSD_data[, keep]
cat("Cells after metabolomics exclusion:", ncol(SSD_data), "\n")

# ---- 3. XIST genotype filtering ---------------------------------------------
cat("\nFiltering ambiguous genotype cells (XIST threshold = 0.2)...\n")
DefaultAssay(SSD_data) <- "RNA"
SSD_data <- JoinLayers(SSD_data)
xist_pos <- WhichCells(SSD_data, expression = XIST > 0.2)
SSD_data$XIST_genotype <- ifelse(colnames(SSD_data) %in% xist_pos, "XX", "XY")

cat("Genotype vs XIST concordance:\n")
print(table(SSD_data@meta.data$genotype, SSD_data$XIST_genotype,
  dnn = c("genotype", "XIST_genotype")
))

SSD_data@graphs    <- list()
SSD_data@neighbors <- list()
keep <- colnames(SSD_data)[
  (SSD_data@meta.data$genotype == "XX" & SSD_data$XIST_genotype == "XX") |
  (SSD_data@meta.data$genotype == "XY" & SSD_data$XIST_genotype == "XY")
]
SSD_data <- SSD_data[, keep]
cat("Cells after XIST filtering:", ncol(SSD_data), "\n")

# ---- 4. Sample exclusion ----------------------------------------------------
cat("\nRemoving low-cell-count samples...\n")
remove_samples <- c("SSD13_tissue13", "SSD04_tissue4")
flag_samples   <- c(
  "SSD42_tissue42", "SSD48_tissue48",
  "SSD14_tissue14", "SSD15_tissue15"
)

SSD_data <- JoinLayers(SSD_data)
keep <- colnames(SSD_data)[!SSD_data@meta.data[[sample_col]] %in% remove_samples]
SSD_data <- SSD_data[, keep]
SSD_data@graphs    <- list()
SSD_data@neighbors <- list()
keep <- colnames(SSD_data)[!SSD_data@meta.data[[sample_col]] %in% remove_samples]
SSD_data <- SSD_data[, keep]

SSD_data@meta.data$low_cell_flag <- SSD_data@meta.data[[sample_col]] %in% flag_samples

cat("Removed:", paste(remove_samples, collapse = ", "), "\n")
cat("Flagged:", paste(flag_samples, collapse = ", "), "\n")
cat("Cells after sample exclusion:", ncol(SSD_data), "\n")

# ---- 5. Reclustering --------------------------------------------------------
cat("\nReclustering...\n")

DefaultAssay(SSD_data) <- "RNA"
SSD_data <- NormalizeData(SSD_data,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
SSD_data <- FindVariableFeatures(SSD_data,
  selection.method = "vst",
  nfeatures = 2000
)
SSD_data <- ScaleData(SSD_data, features = rownames(SSD_data))
SSD_data <- SCTransform(SSD_data, vars.to.regress = "percent.mt", verbose = FALSE)
SSD_data <- RunPCA(SSD_data, features = VariableFeatures(object = SSD_data))
SSD_data <- RunHarmony(SSD_data,
  group.by.vars  = "orig.ident",
  reduction      = "pca",
  reduction.save = "harmony_rna"
)
SSD_data <- RunUMAP(SSD_data,
  reduction      = "harmony_rna",
  dims           = 1:50,
  reduction.name = "umap.rna",
  reduction.key  = "rnaUMAP_"
)
DefaultAssay(SSD_data) <- "SCT"
SSD_data <- FindNeighbors(SSD_data,
  assay     = "SCT",
  reduction = "harmony_rna",
  dims      = 1:50
)
SSD_data <- FindClusters(SSD_data,
  verbose    = FALSE,
  resolution = 0.7,
  graph.name = "SCT_snn"
)
cat("Reclustering complete.\n")

# ---- 6. Save ----------------------------------------------------------------
cat("\nSaving SSD_RNAonly_filtered.rds...\n")
saveRDS(SSD_data, "../../data/SSD_RNAonly_filtered.rds")
cat("Done. Final cell count:", ncol(SSD_data), "\n")
cat("Total time:", format(Sys.time() - t0), "\n")
