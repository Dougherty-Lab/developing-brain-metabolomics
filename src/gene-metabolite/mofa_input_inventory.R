#!/usr/bin/env Rscript
# mofa_input_inventory.R
# ---------------------------------------------------------------------------
# Inventory the sample x cell-type grid that feeds the gene-metabolite analyses,
# ahead of any MOFA/MEFISTO integrative model. Answers:
#
#   1. Exact n of samples with BOTH usable pseudobulk and metabolomics.
#   2. Which cell types are dense enough to be separate MOFA "views" vs.
#      needing collapse to class level.
#   3. Whether the per-cell-type hit ranking is a power artifact.
#   4. What complete-case (imputation removed) will cost in usable n.
#
# Reuses the same cached inputs and the same MIN_CELLS / MIN_SAMPLES gates as
# metabolite_gene_limma.R, so the numbers describe exactly the data that
# pipeline saw. Read-only -- writes CSVs, touches nothing else.
#
# Run from the same src/ subdir as metabolite_gene_limma.R:
#   Rscript mofa_input_inventory.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

source("pseudobulk_functions.R")     # merge_hormone_metadata, .safe_vars
source("metabolite_functions.R")     # read_metabolite_matrix, metab_sample_to_pb, ...

set.seed(123)

# ---- config (mirrors metabolite_gene_limma.R) -------------------------------
cache_dir    <- "../../data/cache"
metab_path   <- "../../results/untargeted/batch2_peak_area_clean.csv"
hormone_xlsx <- "../../doc/targeted/targeted_hormones.xlsx"
out_dir      <- "../../results/gene-metabolite/inventory"

MIN_CELLS   <- 30           # min cells per (sample x cell type)
MIN_SAMPLES <- 10           # min samples per cell type to fit a model
N_MODEL_TERMS <- 5L         # intercept + metab + GW + Sex + batch1_frac

# Which completed run to inventory hits from. Must match the METHOD used in
# metabolite_gene_imputation_sensitivity.qmd for the run being examined.
METHOD_HITS <- "zscore_trim"    # "int" | "zscore" | "zscore_trim"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- cached inputs -----------------------------------------------------------
pb_base     <- readRDS(file.path(cache_dir, "pb_base.rds"))
batch_props <- readRDS(file.path(cache_dir, "batch_props_sample.rds"))

targeted_comparison_raw <- readxl::read_xlsx(hormone_xlsx)
targeted_comparison_raw$E2 <- targeted_comparison_raw$`17B-E2 (ng/mL)`
targeted_comparison_raw$TT <- targeted_comparison_raw$`TT (ng/mL)`
targeted_comparison_raw$P4 <- targeted_comparison_raw$`P4 (ng/mL)`

pb <- build_metab_covariates(
  pb_base, targeted_comparison_raw, batch_props,
  min_cells = MIN_CELLS
)

# ---- metabolite side ---------------------------------------------------------
# Raw matrix: NA = not detected / not picked. No transform -- we only care about
# which columns exist and where the NAs are.
metab_mat <- read_metabolite_matrix(metab_path)

metab_col_for_pb <- metab_sample_to_pb(colnames(metab_mat), pb$metadata$sample_id)

meta <- pb$metadata %>%
  tibble::rownames_to_column("pseudobulk_id") %>%
  dplyr::mutate(
    metab_col = metab_col_for_pb,
    has_metab = !is.na(metab_col)
  )

n_both <- dplyr::n_distinct(meta$sample_id[meta$has_metab])

cat("\n================ HEADLINE ================\n")
cat(sprintf("Pseudobulk rows (sample x cell type, >= %d cells): %d\n",
            MIN_CELLS, nrow(meta)))
cat(sprintf("Distinct samples in pseudobulk:                   %d\n",
            dplyr::n_distinct(meta$sample_id)))
cat(sprintf("Metabolite matrix: %d metabolites x %d samples\n",
            nrow(metab_mat), ncol(metab_mat)))
cat(sprintf("Distinct samples with BOTH assays:                %d   <-- the n\n", n_both))
cat(sprintf("Cell types surviving MIN_CELLS:                   %d\n",
            dplyr::n_distinct(meta$cell_type)))


# ================================================================================
# 1. Per-sample table
# ================================================================================
metab_observed <- tibble::tibble(
  metab_col     = colnames(metab_mat),
  n_observed    = colSums(!is.na(metab_mat)),
  n_missing     = colSums(is.na(metab_mat)),
  frac_observed = n_observed / nrow(metab_mat)
)

sample_tbl <- meta %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    GW              = dplyr::first(GW),
    Sex             = dplyr::first(Sex),
    batch1_frac     = dplyr::first(batch1_frac),
    has_metab       = dplyr::first(has_metab),
    metab_col       = dplyr::first(metab_col),
    n_celltypes_ok  = dplyr::n(),
    total_cells     = sum(cell_count),
    median_cells_ct = median(cell_count),
    .groups = "drop"
  ) %>%
  dplyr::left_join(metab_observed, by = "metab_col") %>%
  dplyr::arrange(dplyr::desc(has_metab), GW, Sex)

readr::write_csv(sample_tbl, file.path(out_dir, "sample_qc_table.csv"))

cat("\n================ PER-SAMPLE ================\n")
print(as.data.frame(sample_tbl), row.names = FALSE)

orphan_metab <- setdiff(colnames(metab_mat), na.omit(meta$metab_col))
if (length(orphan_metab) > 0)
  cat(sprintf("\nMetabolite samples with no usable pseudobulk (%d): %s\n",
              length(orphan_metab), paste(orphan_metab, collapse = ", ")))

# NOTE: sample-level coverage-vs-covariate correlations were removed here. They
# were underpowered (n ~ 22) and blunt. The per-metabolite retention-window table
# in section 5b answers the same question at the resolution that matters: which
# specific metabolites end up tested over a truncated GW range.


# ================================================================================
# 2. Per-cell-type table
# ================================================================================
ct_tbl <- meta %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    n_samples_pb   = dplyr::n(),
    n_samples_both = sum(has_metab),
    n_male_both    = sum(has_metab & Sex == "M"),
    n_female_both  = sum(has_metab & Sex == "F"),
    gw_min         = suppressWarnings(min(GW[has_metab], na.rm = TRUE)),
    gw_max         = suppressWarnings(max(GW[has_metab], na.rm = TRUE)),
    n_gw_levels    = dplyr::n_distinct(GW[has_metab]),
    median_cells   = median(cell_count[has_metab]),
    min_cells      = suppressWarnings(min(cell_count[has_metab])),
    total_cells    = sum(cell_count[has_metab]),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    passes_min_samples = n_samples_both >= MIN_SAMPLES,
    resid_df           = n_samples_both - N_MODEL_TERMS
  ) %>%
  dplyr::arrange(dplyr::desc(n_samples_both))

readr::write_csv(ct_tbl, file.path(out_dir, "celltype_qc_table.csv"))

cat("\n================ PER-CELL-TYPE ================\n")
print(as.data.frame(ct_tbl), row.names = FALSE)
cat(sprintf("\nCell types passing MIN_SAMPLES >= %d: %d of %d\n",
            MIN_SAMPLES, sum(ct_tbl$passes_min_samples), nrow(ct_tbl)))


# ================================================================================
# 3. Presence matrix + view completeness   <-- the view-granularity decision
# ================================================================================
presence <- meta %>%
  dplyr::filter(has_metab) %>%
  dplyr::select(cell_type, sample_id, cell_count) %>%
  tidyr::pivot_wider(names_from = sample_id, values_from = cell_count) %>%
  dplyr::arrange(cell_type)

readr::write_csv(presence, file.path(out_dir, "celltype_sample_presence.csv"))

completeness <- presence %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    n_present    = sum(!is.na(dplyr::c_across(-cell_type))),
    frac_present = n_present / n_both
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(cell_type, n_present, frac_present) %>%
  dplyr::arrange(dplyr::desc(frac_present))

readr::write_csv(completeness, file.path(out_dir, "view_completeness.csv"))

cat("\n================ VIEW COMPLETENESS ================\n")
cat(sprintf("(fraction of the %d dual-assay samples in which each cell type clears %d cells)\n\n",
            n_both, MIN_CELLS))
print(as.data.frame(completeness), row.names = FALSE)
cat("\n  frac_present >= ~0.8  -> viable as its own MOFA view\n")
cat("  frac_present <  ~0.8  -> collapse to class level\n")


# ================================================================================
# 4. Power-confound check (both hit tiers, from the annotated file)
# ================================================================================
suffix    <- if (METHOD_HITS == "int") "" else paste0("-", METHOD_HITS)
hits_path <- file.path("../../results/gene-metabolite",
                       paste0("stability-imputation", suffix),
                       "metabolite_gene_hits_imputation_annotated.csv")

power_tbl <- NULL

if (!file.exists(hits_path)) {
  cat(sprintf("\n[SKIP] Power check: %s not found.\n", hits_path))
} else {
  hits_annot <- readr::read_csv(hits_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      imp_insensitive = !imputation_sensitive %in% TRUE & !too_few_observed %in% TRUE
    )

  cat(sprintf("\nRead %s hits from %s\n",
              format(nrow(hits_annot), big.mark = ","), basename(hits_path)))

  hits_reported <- hits_annot %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      # Tier 1: full global-FDR screen (~7,427)
      n_associations = dplyr::n(),
      n_metabolites  = dplyr::n_distinct(Compound.ID),
      n_genes        = dplyr::n_distinct(gene),
      # Tier 2: imputation-insensitive -- what pathway analysis sees (~2,264)
      n_assoc_clean  = sum(imp_insensitive),
      n_metab_clean  = dplyr::n_distinct(Compound.ID[imp_insensitive]),
      n_genes_clean  = dplyr::n_distinct(gene[imp_insensitive]),
      n_too_few      = sum(too_few_observed %in% TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(retention = n_assoc_clean / n_associations)

  power_tbl <- ct_tbl %>%
    dplyr::inner_join(hits_reported, by = "cell_type") %>%
    dplyr::mutate(
      assoc_per_metab       = n_associations / n_metabolites,
      assoc_per_metab_clean = n_assoc_clean / pmax(n_metab_clean, 1L)
    )

  readr::write_csv(power_tbl, file.path(out_dir, "hits_vs_power.csv"))

  if (nrow(power_tbl) == 0) {
    cat("\n[WARN] cell_type names did not join. ct_tbl has:\n  ",
        paste(ct_tbl$cell_type, collapse = ", "), "\n")
  } else {
    cat("\n================ SCREEN vs CLEAN BY CELL TYPE ================\n")
    print(as.data.frame(
      power_tbl %>%
        dplyr::select(cell_type, n_samples_both, median_cells,
                      n_associations, n_assoc_clean, retention, n_too_few) %>%
        dplyr::arrange(dplyr::desc(n_associations))
    ), row.names = FALSE)

    cat(sprintf("\nRetention range: %.1f%% - %.1f%%  (overall %.1f%%)\n",
                100 * min(power_tbl$retention), 100 * max(power_tbl$retention),
                100 * sum(power_tbl$n_assoc_clean) / sum(power_tbl$n_associations)))
    cat("(Wide spread => the filter reshuffles the ranking, and Figure 5A's order\n")
    cat(" is not the order the pathway analysis rests on.)\n")

    rank_ct <- suppressWarnings(
      cor.test(power_tbl$n_associations, power_tbl$n_assoc_clean, method = "spearman")
    )
    cat(sprintf("\n  rank(screen) ~ rank(clean)      rho = %+.3f   p = %.3g\n",
                rank_ct$estimate, rank_ct$p.value))

    cat("\n================ POWER CONFOUND ================\n")
    for (y in c("n_associations", "n_assoc_clean", "n_metabolites", "retention")) {
      for (v in c("n_samples_both", "median_cells", "total_cells")) {
        ct <- suppressWarnings(
          cor.test(power_tbl[[v]], power_tbl[[y]], method = "spearman")
        )
        cat(sprintf("  %-14s ~ %-15s  rho = %+.3f   p = %.3g\n",
                    y, v, ct$estimate, ct$p.value))
      }
    }
    cat(sprintf("\n  NOTE: n = %d cell types. |rho| must reach ~0.62 for p<0.05.\n",
                nrow(power_tbl)))
    cat("  Treat as effect-size eyeballs, not tests -- they exist to catch a\n")
    cat("  STRONG confound, the only kind that changes the figure decision.\n")
  }
}


# ================================================================================
# 5. Complete-case projection (imputation removed)
# ================================================================================
# With min/5 imputation dropped, each metabolite's usable n becomes the number of
# dual-assay samples in which it was actually observed -- and that n now varies
# per metabolite. This forecasts how many metabolites survive MIN_SAMPLES and
# retain residual df, per cell type.

dual_cols <- na.omit(unique(meta$metab_col))
metab_dual <- metab_mat[, colnames(metab_mat) %in% dual_cols, drop = FALSE]

metab_missing <- tibble::tibble(
  Compound.ID   = rownames(metab_dual),
  n_obs_dual    = rowSums(!is.na(metab_dual)),
  n_miss_dual   = rowSums(is.na(metab_dual)),
  frac_obs_dual = n_obs_dual / ncol(metab_dual),
  median_obs    = apply(metab_dual, 1, median, na.rm = TRUE)
) %>%
  dplyr::arrange(n_obs_dual)

readr::write_csv(metab_missing, file.path(out_dir, "metabolite_missingness.csv"))

cat("\n================ COMPLETE-CASE PROJECTION ================\n")
cat(sprintf("Dual-assay metabolite columns: %d\n", ncol(metab_dual)))
cat(sprintf("Metabolites fully observed:    %d / %d (%.1f%%)\n",
            sum(metab_missing$n_miss_dual == 0), nrow(metab_missing),
            100 * mean(metab_missing$n_miss_dual == 0)))
cat("\nDistribution of per-metabolite observed n (complete-case):\n")
print(summary(metab_missing$n_obs_dual))

cat(sprintf("\nMetabolites with n_obs < MIN_SAMPLES (%d): %d  -> dropped entirely\n",
            MIN_SAMPLES, sum(metab_missing$n_obs_dual < MIN_SAMPLES)))
cat(sprintf("Metabolites with n_obs < %d terms + 1: %d  -> zero residual df\n",
            N_MODEL_TERMS, sum(metab_missing$n_obs_dual < N_MODEL_TERMS + 1L)))

# Missingness mechanism: left-censoring vs. peak-picking failure.
# If low-abundance metabolites are the ones with missing values, missingness is
# left-censoring (MNAR below LOD) and complete-case truncates the low end of the
# range -> attenuates slopes. If missingness is unrelated to abundance, it is a
# detection/alignment artifact and complete-case is the right call.
partial <- metab_missing %>% dplyr::filter(n_obs_dual >= 3)
if (nrow(partial) >= 10) {
  mech <- suppressWarnings(
    cor.test(partial$median_obs, partial$frac_obs_dual, method = "spearman")
  )
  cat(sprintf("\n---- Missingness mechanism ----\n"))
  cat(sprintf("  frac_observed ~ median observed intensity   rho = %+.3f   p = %.3g\n",
              mech$estimate, mech$p.value))
  cat("  Strong POSITIVE rho => low-abundance metabolites are the missing ones\n")
  cat("    => left-censoring dominates; complete-case truncates the low range\n")
  cat("       and will attenuate slopes for exactly those metabolites.\n")
  cat("  Near ZERO rho      => missingness is unrelated to abundance\n")
  cat("    => detection/alignment artifact; dropping is the right call.\n")
}

# ---- 5b. Per-metabolite retention window -------------------------------------
# The number that actually matters. For each metabolite, complete-case tests it
# only in the samples where it was detected. This reports the GW window and sex
# balance of THAT retained set, per metabolite.
#
# Why: if metabolite A is tested over GW16-24 and metabolite B only over GW20-24,
# their slopes are not the same estimand -- B's is a within-window slope. That
# breaks logFC comparability across the ~569 predictors, which is the property
# zscore_trim was chosen to preserve, and therefore breaks MSEA ranking validity.
# This is a per-metabolite column you can filter/report, not a coefficient to
# argue about.

covar_by_col <- sample_tbl %>%
  dplyr::filter(has_metab) %>%
  dplyr::select(metab_col, GW, Sex, batch1_frac)

gw_vec    <- covar_by_col$GW[match(colnames(metab_dual), covar_by_col$metab_col)]
sex_vec   <- covar_by_col$Sex[match(colnames(metab_dual), covar_by_col$metab_col)]
batch_vec <- covar_by_col$batch1_frac[match(colnames(metab_dual), covar_by_col$metab_col)]

obs_mask <- !is.na(metab_dual)

# GW of retained (observed) samples; GW of dropped samples, for contrast.
gw_mat            <- matrix(gw_vec, nrow(metab_dual), ncol(metab_dual), byrow = TRUE)
gw_kept           <- gw_mat; gw_kept[!obs_mask] <- NA
gw_dropped        <- gw_mat; gw_dropped[obs_mask] <- NA

batch_mat         <- matrix(batch_vec, nrow(metab_dual), ncol(metab_dual), byrow = TRUE)
batch_kept        <- batch_mat; batch_kept[!obs_mask] <- NA

is_f              <- matrix(sex_vec == "F", nrow(metab_dual), ncol(metab_dual), byrow = TRUE)

gw_full_range <- diff(range(gw_vec, na.rm = TRUE))

row_min <- function(m) suppressWarnings(apply(m, 1, min, na.rm = TRUE))
row_max <- function(m) suppressWarnings(apply(m, 1, max, na.rm = TRUE))

metab_window <- tibble::tibble(
  Compound.ID   = rownames(metab_dual),
  n_obs_dual    = rowSums(obs_mask),
  frac_obs_dual = rowSums(obs_mask) / ncol(metab_dual),
  gw_min_kept   = row_min(gw_kept),
  gw_max_kept   = row_max(gw_kept),
  gw_mean_kept  = rowMeans(gw_kept, na.rm = TRUE),
  n_gw_kept     = apply(gw_kept, 1, function(z) dplyr::n_distinct(z[!is.na(z)])),
  gw_mean_dropped = rowMeans(gw_dropped, na.rm = TRUE),
  n_F_kept      = rowSums(obs_mask & is_f),
  n_M_kept      = rowSums(obs_mask & !is_f),
  batch1_mean_kept = rowMeans(batch_kept, na.rm = TRUE)
) %>%
  dplyr::mutate(
    gw_min_kept   = ifelse(is.finite(gw_min_kept), gw_min_kept, NA_real_),
    gw_max_kept   = ifelse(is.finite(gw_max_kept), gw_max_kept, NA_real_),
    gw_span_kept  = gw_max_kept - gw_min_kept,
    gw_span_frac  = gw_span_kept / gw_full_range,
    # GW shift: positive => retained samples skew LATE (early samples dropped)
    gw_shift      = gw_mean_kept - mean(gw_vec, na.rm = TRUE),
    frac_F_kept   = n_F_kept / pmax(n_F_kept + n_M_kept, 1L),
    # Flags for filtering / supplement reporting
    truncated_gw  = !is.na(gw_span_frac) & gw_span_frac < 0.75,
    sex_skewed    = !is.na(frac_F_kept) & (frac_F_kept < 0.25 | frac_F_kept > 0.75),
    below_min_n   = n_obs_dual < MIN_SAMPLES
  ) %>%
  dplyr::left_join(
    metab_missing %>% dplyr::select(Compound.ID, median_obs),
    by = "Compound.ID"
  ) %>%
  dplyr::arrange(gw_span_frac, n_obs_dual)

# Attach metabolite names if available.
nl <- attr(metab_mat, "name_lookup")
if (!is.null(nl))
  metab_window <- metab_window %>%
    dplyr::mutate(Name = unname(nl[Compound.ID]), .after = Compound.ID)

readr::write_csv(metab_window, file.path(out_dir, "metabolite_retention_window.csv"))

cat("\n================ PER-METABOLITE RETENTION WINDOW ================\n")
cat(sprintf("Full GW range across dual-assay samples: %.1f - %.1f (span %.1f)\n",
            min(gw_vec, na.rm = TRUE), max(gw_vec, na.rm = TRUE), gw_full_range))

cat("\nDistribution of gw_span_frac (fraction of full GW range each metabolite retains):\n")
print(summary(metab_window$gw_span_frac))

cat(sprintf("\n  Full GW span retained (>=0.95):  %d / %d (%.1f%%)\n",
            sum(metab_window$gw_span_frac >= 0.95, na.rm = TRUE), nrow(metab_window),
            100 * mean(metab_window$gw_span_frac >= 0.95, na.rm = TRUE)))
cat(sprintf("  TRUNCATED window  (<0.75):       %d (%.1f%%)   <-- slopes not comparable\n",
            sum(metab_window$truncated_gw), 100 * mean(metab_window$truncated_gw)))
cat(sprintf("  Sex-skewed (<25%% or >75%% F):     %d (%.1f%%)\n",
            sum(metab_window$sex_skewed), 100 * mean(metab_window$sex_skewed)))
cat(sprintf("  Below MIN_SAMPLES (%d):           %d (%.1f%%)\n",
            MIN_SAMPLES, sum(metab_window$below_min_n),
            100 * mean(metab_window$below_min_n)))

cat("\nDirection of GW shift among metabolites with ANY missingness:\n")
shifted <- metab_window %>% dplyr::filter(n_obs_dual < ncol(metab_dual), n_obs_dual >= 3)
if (nrow(shifted) > 0) {
  cat(sprintf("  n = %d metabolites; mean gw_shift = %+.2f weeks (range %+.2f to %+.2f)\n",
              nrow(shifted), mean(shifted$gw_shift, na.rm = TRUE),
              min(shifted$gw_shift, na.rm = TRUE), max(shifted$gw_shift, na.rm = TRUE)))
  cat(sprintf("  Retained set skews LATE in %d / %d (%.1f%%)\n",
              sum(shifted$gw_shift > 0, na.rm = TRUE), nrow(shifted),
              100 * mean(shifted$gw_shift > 0, na.rm = TRUE)))
  cat("  (Consistently positive => early-GW samples are the ones being dropped,\n")
  cat("   i.e. lower tissue yield -> fewer peaks clear detection.)\n")
}

cat("\n25 most truncated metabolites:\n")
print(as.data.frame(
  metab_window %>%
    dplyr::select(dplyr::any_of(c("Name", "Compound.ID")), n_obs_dual,
                  gw_min_kept, gw_max_kept, gw_span_frac, gw_shift, frac_F_kept) %>%
    head(25)
), row.names = FALSE)

# Per cell type x metabolite: usable n after complete-case.
cc_grid <- meta %>%
  dplyr::filter(has_metab) %>%
  dplyr::select(cell_type, metab_col)

obs_long <- which(!is.na(metab_dual), arr.ind = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    Compound.ID = rownames(metab_dual)[row],
    metab_col   = colnames(metab_dual)[col]
  ) %>%
  dplyr::select(Compound.ID, metab_col)

cc_ct <- cc_grid %>%
  dplyr::inner_join(obs_long, by = "metab_col", relationship = "many-to-many") %>%
  dplyr::count(cell_type, Compound.ID, name = "n_usable") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    n_metab_testable   = sum(n_usable >= MIN_SAMPLES),
    n_metab_zero_df    = sum(n_usable < N_MODEL_TERMS + 1L),
    median_usable_n    = median(n_usable),
    min_usable_n       = min(n_usable),
    .groups = "drop"
  ) %>%
  dplyr::left_join(ct_tbl %>% dplyr::select(cell_type, n_samples_both),
                   by = "cell_type") %>%
  dplyr::arrange(dplyr::desc(n_metab_testable))

readr::write_csv(cc_ct, file.path(out_dir, "complete_case_by_celltype.csv"))

cat("\n---- Testable metabolites per cell type under complete-case ----\n")
print(as.data.frame(cc_ct), row.names = FALSE)
cat(sprintf("\n(Compare n_metab_testable against the %d metabolites the imputed run tested.)\n",
            nrow(metab_missing)))


cat(sprintf("\n\nWrote CSVs to %s\n", out_dir))
cat("  sample_qc_table.csv            -- the n, per-sample coverage\n")
cat("  celltype_qc_table.csv          -- per-cell-type n, depth, resid df\n")
cat("  celltype_sample_presence.csv   -- cell_type x sample grid (cell_count)\n")
cat("  view_completeness.csv          -- MOFA view-granularity decision\n")
cat("  hits_vs_power.csv              -- both hit tiers vs. power\n")
cat("  metabolite_missingness.csv     -- per-metabolite observed n\n")
cat("  metabolite_retention_window.csv -- per-metabolite GW window + sex balance\n")
cat("  complete_case_by_celltype.csv  -- forecast for the imputation-free run\n")
