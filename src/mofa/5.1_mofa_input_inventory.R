#!/usr/bin/env Rscript
# mofa_input_inventory.R
# ---------------------------------------------------------------------------
# Inventory the sample x cell-type grid for MOFA integration. Reports exact
# sample counts with paired pseudobulk + metabolomics, cell type density,
# view completeness, and hit-tier power-confound checks.
#
# Inputs:  data/cache/pb_base.rds, data/cache/batch_props_sample.rds,
#          results/untargeted/peak_area_clean.csv, doc/targeted/targeted_hormones.xlsx
# Outputs: CSVs in results/gene-metabolite/inventory/
#
# Upstream:  prep_data.R, Filtering.Rmd
# Downstream: mofa_integration.Rmd
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

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
source(file.path(root, "src/gene-hormone/pseudobulk_functions.R"))     # merge_hormone_metadata, .safe_vars
source(file.path(root, "src/gene-metabolite/metabolite_functions.R"))     # read_metabolite_matrix, metab_sample_to_pb, ...

set.seed(123)

# ---- config (mirrors metabolite_gene_limma.R) -------------------------------
cache_dir    <- file.path(root, "data/cache")
metab_path   <- file.path(root, "results/untargeted/peak_area_clean.csv")
hormone_xlsx <- file.path(root, "doc/targeted/targeted_hormones.xlsx")
out_dir      <- file.path(root, "results/gene-metabolite/inventory")

MIN_CELLS   <- 30           # min cells per (sample x cell type)
MIN_SAMPLES <- 10           # min samples per cell type to fit a model
N_MODEL_TERMS <- 5L         # intercept + metab + GW + Sex + batch1_frac

# Which completed run to inventory hits from. Must match the METHOD used in
# metabolite_gene_imputation_sensitivity.qmd for the run being examined.
METHOD_HITS <- "log2_na"    # "int" | "zscore" | "zscore_trim"

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
  dplyr::arrange(dplyr::desc(has_metab), GW, Sex)

readr::write_csv(sample_tbl, file.path(out_dir, "sample_qc_table.csv"))

cat("\n================ PER-SAMPLE ================\n")
print(as.data.frame(sample_tbl), row.names = FALSE)

orphan_metab <- setdiff(colnames(metab_mat), na.omit(meta$metab_col))
if (length(orphan_metab) > 0)
  cat(sprintf("\nMetabolite samples with no usable pseudobulk (%d): %s\n",
              length(orphan_metab), paste(orphan_metab, collapse = ", ")))


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
hits_path <- file.path(file.path(root, "results/gene-metabolite"),
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


cat(sprintf("\n\nWrote CSVs to %s\n", out_dir))
cat("  sample_qc_table.csv            -- the n, per-sample coverage\n")
cat("  celltype_qc_table.csv          -- per-cell-type n, depth, resid df\n")
cat("  celltype_sample_presence.csv   -- cell_type x sample grid (cell_count)\n")
cat("  view_completeness.csv          -- MOFA view-granularity decision\n")
cat("  hits_vs_power.csv              -- both hit tiers vs. power\n")

# ---- AI assistance disclosure ------------------------------------------------
# Code in this script was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.

# ---- session info ------------------------------------------------------------
cat("\n\n---- Session Info ----\n")
print(sessionInfo())

