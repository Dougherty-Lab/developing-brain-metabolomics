# prep_data.R
#
# Loads SSD_RNAonly_filtered.rds (~44G) ONCE and caches the small derived
# objects that subsequent scripts need

# Outputs (all in ../../data/cache/):
#   pb_base.rds              - build_pseudobulk() output, celltype level
#   pb_base_subclass.rds     - build_pseudobulk() output, subclass level
#   pb_base_class.rds        - build_pseudobulk() output, class level
#   batch_props_sample.rds   - per-sample batch1_frac (Gene_Covariate_Testing)
#   SSD_meta_umap.rds        - DietSeurat: meta.data + umap.rna only, 1 gene
#                              (testosterone / estradiol / progesterone —
#                               add_batch_fraction() + plot_analysis_umap())
#   receptor_expr_object.rds - DietSeurat: meta.data + umap.rna + counts/data
#                              for the 18 receptor genes only
#   pb_comparison_subset.rds - counts_sub + sample_labels for
#                              IN-CGE-Immature x SSD07/SSD47
#                              (pseudobulk_cell_comparison.qmd)
Packages <- c("tidyverse", "Seurat", "harmony")
lapply(Packages, library, character.only = TRUE)
setwd("/scratch/jdlab/sneha/developing-brain-metabolomics/src/gene-hormone/")

set.seed(123)

sample_col <- "Sample"

# ---- 1. Load ----------------------------------------------------------------
cat("Loading SSD_RNAonly_filtered.rds (this is the only full load)...\n")
t0 <- Sys.time()
SSD_data <- readRDS("../../data/SSD_RNAonly_filtered.rds")
cat("Load time:", format(Sys.time() - t0), "\n")

celltype_col <- "celltype"
sample_col   <- "Sample"

source("pseudobulk_functions.R")

cache_dir <- "../../data/cache"
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