# prep_data.R
# ---------------------------------------------------------------------------
# Loads the filtered Seurat object once and caches small derived objects for
# all downstream scripts: pseudobulk counts at three annotation levels,
# batch proportions, lightweight Seurat subsets for UMAP and receptor plots,
# and a pseudobulk comparison subset.
#
# Inputs:  data/SSD_RNAonly_filtered.rds
# Outputs: data/cache/pb_base.rds, pb_base_subclass.rds, pb_base_class.rds,
#          batch_props_sample.rds, SSD_meta_umap.rds, receptor_expr_object.rds,
#          pb_comparison_subset.rds
#
# Upstream:  filter_recluster.R
# Downstream: All QMD scripts in src/gene-hormone/
#
# Run: Rscript prep_data.R (or via sbatch; ~200GB RAM recommended)
# ---------------------------------------------------------------------------
Packages <- c("tidyverse", "Seurat", "harmony")
lapply(Packages, library, character.only = TRUE)
<<<<<<< HEAD
setwd("/scratch/jdlab/emma/developing-brain-metabolomics/src/gene-hormone/")
=======
>>>>>>> f939efd08d50f1978b54af5677f8e75f12309a7b

set.seed(123)

sample_col <- "Sample"

# ---- 1. Load ----------------------------------------------------------------
cat("Loading SSD_RNAonly_filtered.rds (this is the only full load)...\n")
t0 <- Sys.time()
SSD_data <- readRDS(file.path(root, "data/SSD_RNAonly_filtered.rds"))
cat("Load time:", format(Sys.time() - t0), "\n")

celltype_col <- "celltype"
sample_col   <- "Sample"

# Bootstrap project root
find_project_root <- function(marker = ".git") {
  d <- normalizePath(getwd())
  repeat {
    if (dir.exists(file.path(d, marker))) return(d)
    parent <- dirname(d)
    if (parent == d) stop("Project root not found (no ", marker, " above ", getwd(), ")")
    d <- parent
  }
}
root <- find_project_root()

source(file.path(root, "src/gene-hormone/pseudobulk_functions.R"))

cache_dir <- file.path(root, "data/cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 1. pb_base: shared pseudobulk (Gene_Covariate_Testing, TT, E2/P4) -----
cat("\nBuilding pseudobulk...\n")
pb_base <- build_pseudobulk(
  SSD_data,
  sample_col   = sample_col,
  celltype_col = celltype_col
)
saveRDS(pb_base, file.path(cache_dir, "pb_base.rds"))
cat("Saved pb_base.rds\n")

pb_base_subclass <- build_pseudobulk(
  SSD_data,
  sample_col   = sample_col,
  celltype_col = "subclass"
)
saveRDS(pb_base_subclass, file.path(cache_dir, "pb_base_subclass.rds"))
cat("Saved pb_base_subclass.rds\n")

pb_base_class <- build_pseudobulk(
  SSD_data,
  sample_col   = sample_col,
  celltype_col = "class"
)
saveRDS(pb_base_class, file.path(cache_dir, "pb_base_class.rds"))
cat("Saved pb_base_class.rds\n")

# ---- 2. batch_props_sample (Gene_Covariate_Testing batch section) ---------
batch_props_sample <- SSD_data@meta.data %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(
    batch1_frac = mean(batch == "batch1"),
    .groups     = "drop"
  )
saveRDS(batch_props_sample, file.path(cache_dir, "batch_props_sample.rds"))
cat("Saved batch_props_sample.rds\n")

# ---- 3. SSD_meta_umap: meta.data + umap.rna only (TT / E2 / P4 scripts) ---
# Needed for add_batch_fraction() (reads @meta.data$Sample, $batch) and
# plot_analysis_umap() (reads @meta.data, Embeddings(), DimPlot()).

cat("\nBuilding minimal meta+UMAP object...\n")
DefaultAssay(SSD_data) <- "RNA"
SSD_meta_umap <- DietSeurat(
  SSD_data,
  layers    = "counts",
  assays    = "RNA",
  dimreducs = "umap.rna"
)
SSD_meta_umap <- subset(SSD_meta_umap, features = rownames(SSD_meta_umap)[1])
saveRDS(SSD_meta_umap, file.path(cache_dir, "SSD_meta_umap.rds"))
cat("Saved SSD_meta_umap.rds  (", format(object.size(SSD_meta_umap), units = "MB"), ")\n")

# ---- 4. receptor_expr_object: 18 receptor genes only -----------------------
receptor_genes <- c(
  "ESR1", "ESR2", "PGR", "AR",
  "FSHR", "LHCGR", "GNRHR",
  "PTGER1", "PTGER2", "PTGER3", "PTGER4",
  "THRA", "THRB",
  "PAQR7", "PAQR8", "PAQR5", "PGRMC1", "PGRMC2"
)
receptor_genes_present <- intersect(receptor_genes, rownames(SSD_data))
cat(
  "\nReceptor genes present:", length(receptor_genes_present),
  "/", length(receptor_genes), "\n"
)
cat("Missing:", setdiff(receptor_genes, receptor_genes_present), "\n")

receptor_expr_object <- DietSeurat(
  SSD_data,
  layers    = c("counts", "data"),
  assays    = "RNA",
  dimreducs = "umap.rna"
)
receptor_expr_object <- subset(receptor_expr_object, features = receptor_genes_present)
saveRDS(receptor_expr_object, file.path(cache_dir, "receptor_expr_object.rds"))
cat("Saved receptor_expr_object.rds  (", format(object.size(receptor_expr_object), units = "MB"), ")\n")

# ---- 5. pb_comparison_subset: raw per-cell counts for the comparison ------
TARGET_CELLTYPE <- "EN-IT-Immature"
TARGET_SAMPLES  <- c("SSD10_tissue10", "SSD43_tissue43")

meta <- SSD_data@meta.data
cell_idx <- which(
  meta$celltype == TARGET_CELLTYPE &
    meta$Sample %in% TARGET_SAMPLES
)

counts_sub    <- GetAssayData(SSD_data, assay = "RNA", layer = "counts")[, cell_idx, drop = FALSE]
sample_labels <- meta$Sample[cell_idx]

saveRDS(
  list(counts_sub = counts_sub, sample_labels = sample_labels),
  file.path(cache_dir, "pb_comparison_subset.rds")
)
cat("Saved pb_comparison_subset.rds\n")

cat("\nDone. All cached objects written to", cache_dir, "\n")
# ---- AI assistance disclosure ------------------------------------------------
# Code in this script was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.

# ---- session info ------------------------------------------------------------
cat("\n\n---- Session Info ----\n")
print(sessionInfo())
