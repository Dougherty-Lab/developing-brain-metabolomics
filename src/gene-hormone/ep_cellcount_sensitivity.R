# ep_cellcount_sensitivity.R
#
# Sensitivity analysis: E2 and P4 hormone-gene associations with cell_count
# included as an additive covariate in all four models, at the celltype level,
# sex-combined only.
#
# Rationale: pseudobulks aggregated from fewer cells are noisier estimates of
# the true cell-type expression profile, independent of voom's depth
# normalisation. Including cell_count tests whether any hormone signal is
# attributable to aggregation reliability rather than hormone concentration.
# Raw cell_count is used (not log-transformed) because the justification is
# aggregation noise, not library depth (which voom already accounts for via
# precision weights; Law et al. 2014, Genome Biol).
#
# Models (per hormone):
#   cc_base           ~ hormone + cell_count
#   cc_sex            ~ hormone + Sex + cell_count
#   cc_batch          ~ hormone + batch1_frac + cell_count
#   cc_sex_batch      ~ hormone + Sex + batch1_frac + cell_count
#
# Outputs (../../results/gene-hormone/cell_count_sensitivity/):
#   estradiol/   — CSVs + summary plots for E2 x 4 models
#   progesterone/ — CSVs + summary plots for P4 x 4 models

# ---- Packages ---------------------------------------------------------------
Packages <- c(
  "tidyverse", "Seurat", "Signac", "glue", "readxl", "cowplot",
  "edgeR", "limma", "patchwork", "svglite", "Matrix", "ggpubr"
)
lapply(Packages, library, character.only = TRUE)

setwd("/scratch/jdlab/sneha/developing-brain-metabolomics/src/gene-hormone/")
source("pseudobulk_functions.R")

set.seed(123)

# ---- Output directories -----------------------------------------------------
out_base <- "../../results/gene-hormone/cell_count_sensitivity"
out_e2   <- file.path(out_base, "estradiol")
out_p4   <- file.path(out_base, "progesterone")

for (d in c(out_e2, out_p4,
            file.path(out_e2, "svg"),
            file.path(out_p4, "svg"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- Local model wrappers (cell_count added) --------------------------------
# cell_count must be present in metadata; it is populated by build_pseudobulk()
# and survives merge_hormone_metadata() (minimum value = 30 after the min_cells
# filter). No transformation applied — see header rationale.

run_voom_lm_cc <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  formula_vars <- .safe_vars(c(predictor, "cell_count"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v      <- limma::voom(counts, design)
  fit    <- limma::eBayes(limma::lmFit(v, design))
  limma::topTable(fit, coef = colnames(design)[2],
                  number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_sex_cc <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  formula_vars <- .safe_vars(c(predictor, "Sex", "cell_count"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v      <- limma::voom(counts, design)
  fit    <- limma::eBayes(limma::lmFit(v, design))
  limma::topTable(fit, coef = colnames(design)[2],
                  number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_batch_cc <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"batch1_frac" %in% colnames(metadata))
    stop("batch1_frac not found. Run add_batch_fraction() first.")
  formula_vars <- .safe_vars(c(predictor, "batch1_frac", "cell_count"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v      <- limma::voom(counts, design)
  fit    <- limma::eBayes(limma::lmFit(v, design))
  limma::topTable(fit, coef = colnames(design)[2],
                  number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_sex_batch_cc <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"batch1_frac" %in% colnames(metadata))
    stop("batch1_frac not found. Run add_batch_fraction() first.")
  formula_vars <- .safe_vars(c(predictor, "Sex", "batch1_frac", "cell_count"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v      <- limma::voom(counts, design)
  fit    <- limma::eBayes(limma::lmFit(v, design))
  limma::topTable(fit, coef = colnames(design)[2],
                  number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

# ---- Load data --------------------------------------------------------------
cat("Loading cached objects...\n")
SSD_data <- readRDS("../../data/cache/SSD_meta_umap.rds")
pb_base  <- readRDS("../../data/cache/pb_base.rds")

targeted_comparison_raw <- read_xlsx("../../doc/targeted/targeted_hormones.xlsx")
targeted_comparison_raw$E2 <- targeted_comparison_raw$`17B-E2 (ng/mL)`
targeted_comparison_raw$TT <- targeted_comparison_raw$`TT (ng/mL)`
targeted_comparison_raw$P4 <- targeted_comparison_raw$`P4 (ng/mL)`

pb_sample_pattern <- "tissue\\d+"

# ---- Outlier removal (identical to primary analysis) ------------------------
detect_outliers <- function(data, variable) {
  outlier_limit <- quantile(data[[variable]], 0.75, na.rm = TRUE) +
    1.5 * IQR(data[[variable]], na.rm = TRUE)
  data %>% dplyr::filter(.data[[variable]] <= outlier_limit)
}

targeted_E2 <- detect_outliers(targeted_comparison_raw, "E2")
targeted_P4 <- detect_outliers(targeted_comparison_raw, "P4")

cat("E2 — before:", nrow(targeted_comparison_raw),
    "after:", nrow(targeted_E2), "\n")
cat("P4 — before:", nrow(targeted_comparison_raw),
    "after:", nrow(targeted_P4), "\n")

# ---- Pseudobulk + batch fraction (identical to primary analysis) ------------
pb_E2 <- merge_hormone_metadata(
  pb_base, targeted_E2,
  pb_sample_pattern = pb_sample_pattern, min_cells = 30
) %>% add_batch_fraction(SSD_data)

pb_P4 <- merge_hormone_metadata(
  pb_base, targeted_P4,
  pb_sample_pattern = pb_sample_pattern, min_cells = 30
) %>% add_batch_fraction(SSD_data)

cat("E2 pseudobulk samples:", nrow(pb_E2$metadata),
    "| NAs in cell_count:", sum(is.na(pb_E2$metadata$cell_count)), "\n")
cat("P4 pseudobulk samples:", nrow(pb_P4$metadata),
    "| NAs in cell_count:", sum(is.na(pb_P4$metadata$cell_count)), "\n")

# ---- Helper: run all four cc models for one hormone -------------------------
run_cc_models <- function(pb, hormone, out_dir) {

  models <- list(
    list(label = "cc_base",      fn = run_voom_lm_cc,          tag = paste0("~ ", hormone, " + cell_count")),
    list(label = "cc_sex",       fn = run_voom_lm_sex_cc,      tag = paste0("~ ", hormone, " + Sex + cell_count")),
    list(label = "cc_batch",     fn = run_voom_lm_batch_cc,    tag = paste0("~ ", hormone, " + batch1_frac + cell_count")),
    list(label = "cc_sex_batch", fn = run_voom_lm_sex_batch_cc,tag = paste0("~ ", hormone, " + Sex + batch1_frac + cell_count"))
  )

  results_all <- list()

  for (m in models) {
    cat("\n--- Running", hormone, m$label, "---\n")

    res <- analyze_hormone_combined(
      pb$counts, pb$metadata,
      hormone     = hormone,
      model       = m$fn,
      min_samples = 21
    )

    summary_df <- summarize_sig_genes(res$combined_results)

    # Save CSV
    write_csv(
      res$combined_results,
      file.path(out_dir, paste0(hormone, "_", m$label, "_results.csv"))
    )
    write_csv(
      summary_df,
      file.path(out_dir, paste0(hormone, "_", m$label, "_summary.csv"))
    )

    # Save summary plot
    save_dual_format(
      plot_sig_gene_summary(summary_df, paste0(hormone, " — ", m$tag)),
      output_dir    = out_dir,
      filename_base = paste0(hormone, "_", m$label, "_sig_gene_summary")
    )

    results_all[[m$label]] <- list(results = res, summary = summary_df,
                                   tag = m$tag)
    cat("Sig genes (FDR < 0.10, |logFC| >= 0.25):",
        sum(summary_df$num_sig_genes), "\n")
  }

  results_all
}

# ---- Run E2 -----------------------------------------------------------------
cat("\n===== ESTRADIOL =====\n")
res_E2 <- run_cc_models(pb_E2, "E2", out_e2)

# Cross-model summary for E2
cross_E2 <- purrr::map_dfr(names(res_E2), function(nm) {
  res_E2[[nm]]$summary %>%
    dplyr::mutate(model = res_E2[[nm]]$tag)
})
write_csv(cross_E2, file.path(out_e2, "E2_cross_model_summary.csv"))
cat("\nE2 cross-model totals:\n")
print(
  cross_E2 %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(total_sig = sum(num_sig_genes), .groups = "drop"),
  n = Inf
)

# ---- Run P4 -----------------------------------------------------------------
cat("\n===== PROGESTERONE =====\n")
res_P4 <- run_cc_models(pb_P4, "P4", out_p4)

# Cross-model summary for P4
cross_P4 <- purrr::map_dfr(names(res_P4), function(nm) {
  res_P4[[nm]]$summary %>%
    dplyr::mutate(model = res_P4[[nm]]$tag)
})
write_csv(cross_P4, file.path(out_p4, "P4_cross_model_summary.csv"))
cat("\nP4 cross-model totals:\n")
print(
  cross_P4 %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(total_sig = sum(num_sig_genes), .groups = "drop"),
  n = Inf
)

cat("\nDone. Results written to", out_base, "\n")
