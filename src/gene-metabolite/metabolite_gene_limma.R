#!/usr/bin/env Rscript
# metabolite_gene_limma.R
# ---------------------------------------------------------------------------
# Loop ~569 metabolites through limma-voom gene-ASSOCIATION models and stream
# one parquet file per metabolite to results/gene-metabolite/parquet/.
#
#   response  : pseudobulk gene counts (voom)
#   predictor : standardised metabolite abundance  (log2 -> min/5 impute -> z)
#   design    : ~ metab + GW + Sex + batch1_frac   (split by cell type)
#
# Reuses cached pb_base.rds + batch_props_sample.rds, so the full ~45 GB Seurat
# object is never loaded. Writes RAW, UNFILTERED results only -- thresholding,
# global-FDR, annotation, and figures happen in metabolite_gene_associations.qmd
# so screen changes never re-run this loop.
#
# Run:  Rscript metabolite_gene_limma.R
#       (light enough to run interactively; wrap in an sbatch script if you'd
#        rather detach it -- no large-memory allocation is required.)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(arrow)
})

source("pseudobulk_functions.R")     # build_pseudobulk, merge_hormone_metadata, .safe_vars
source("metabolite_functions.R")     # metabolite transform + analysis helpers

set.seed(123)

# ---- paths / parameters ------------------------------------------------------
cache_dir   <- "../../data/cache"
metab_path  <- "../../results/untargeted/batch2_peak_area_clean.csv"
hormone_xlsx <- "../../doc/targeted/targeted_hormones.xlsx"

# Predictor transform for this run. "int" (leverage-robust, current run) writes
# to .../parquet; any other method writes to .../parquet-<method> so a second
# run cannot clobber the first.
# log2_na: log2(x+1) -> 0.75*IQR outlier removal -> z-score, no imputation.
# Missing and outlier samples stay NA and are dropped per-metabolite by
# analyze_metabolite()'s !is.na(md$metab) filter. MIN_SAMPLES=10 post-NA-drop.
METHOD  <- "log2_na"                           # "int" | "zscore" | "zscore_trim" | "log2_na"
out_dir <- if (METHOD == "int") {
  "../../results/gene-metabolite/parquet"
} else {
  sprintf("../../results/gene-metabolite/parquet-%s", METHOD)
}
MIN_CELLS   <- 30    # min cells per (sample x cell type), matches hormone pipeline
MIN_SAMPLES <- 10    # min samples per cell-type group to fit a model

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- cached inputs -----------------------------------------------------------
pb_base     <- readRDS(file.path(cache_dir, "pb_base.rds"))
batch_props <- readRDS(file.path(cache_dir, "batch_props_sample.rds"))

# ---- targeted hormone sheet -> Sex / GW source -------------------------------
targeted_comparison_raw <- readxl::read_xlsx(hormone_xlsx)
targeted_comparison_raw$E2 <- targeted_comparison_raw$`17B-E2 (ng/mL)`
targeted_comparison_raw$TT <- targeted_comparison_raw$`TT (ng/mL)`
targeted_comparison_raw$P4 <- targeted_comparison_raw$`P4 (ng/mL)`

# ---- covariate metadata (Sex / GW / batch1_frac on pb sample set) ------------
pb <- build_metab_covariates(
  pb_base, targeted_comparison_raw, batch_props,
  min_cells = MIN_CELLS
)
cat(sprintf("Pseudobulk samples: %d | cell types: %d\n",
            nrow(pb$metadata), dplyr::n_distinct(pb$metadata$cell_type)))

# ---- metabolite matrix: transform per METHOD (see top) -----------------------
metab_mat   <- read_metabolite_matrix(metab_path)
metab_z     <- transform_metabolite_matrix(metab_mat, method = METHOD)
name_lookup <- attr(metab_z, "name_lookup")
cat(sprintf("Transform: %s | metabolites: %d | metabolite samples: %d\n",
            METHOD, nrow(metab_z), ncol(metab_z)))

# Map metabolite SampleNN columns -> pb sample order (NA where a pb sample has
# no metabolite measurement). Built once; reused for every metabolite.
metab_col_for_pb <- metab_sample_to_pb(colnames(metab_z), pb$metadata$sample_id)
if (all(is.na(metab_col_for_pb)))
  stop("No metabolite columns mapped to pb sample_ids -- check Sample<->tissue keys.")
cat(sprintf("pb samples matched to a metabolite column: %d / %d\n",
            sum(!is.na(metab_col_for_pb)), length(metab_col_for_pb)))

# ---- loop over metabolites ---------------------------------------------------
compound_ids <- rownames(metab_z)
n_written    <- 0L
t0           <- Sys.time()

for (i in seq_along(compound_ids)) {
  cid <- compound_ids[i]

  # Standardised value for THIS metabolite, aligned to pb sample order.
  # match() -> integer positions with numeric NA for pb samples that have no
  # metabolite column; matrix indexing tolerates numeric NA (character NA errors),
  # yielding NA metab values that analyze_metabolite() drops via !is.na(md$metab).
  col_idx <- match(metab_col_for_pb, colnames(metab_z))
  pb$metadata$metab <- as.numeric(metab_z[cid, col_idx])

  res <- analyze_metabolite(
    pb$counts, pb$metadata,
    compound_id = cid,
    name        = unname(name_lookup[cid]),
    min_samples = MIN_SAMPLES
  )

  if (!is.null(res) && nrow(res) > 0) {
    fn <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", cid), ".parquet"))
    arrow::write_parquet(res, fn)
    n_written <- n_written + 1L
  }

  if (i %% 50 == 0 || i == length(compound_ids))
    cat(sprintf("[%s] %d/%d metabolites processed (%d with results)\n",
                format(Sys.time() - t0), i, length(compound_ids), n_written))
}

cat(sprintf("\nDone. %d metabolite parquet files written to %s\n",
            n_written, out_dir))
