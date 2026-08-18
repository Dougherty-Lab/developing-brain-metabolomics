# pseudobulk_functions.R
#
# Shared functions for pseudobulk + limma-voom DE analysis of snRNA-seq data.
# Sourced by manual hormone_gene_limma_sexseparate.Rmd and sex-combined scripts.
#
# Design:
#   build_pseudobulk()         -> pure aggregation; no hormone/phenotype knowledge.
#                                 Returns counts (samples x genes), logcounts, metadata.
#   merge_hormone_metadata()   -> joins hormone phenotype data onto pb$metadata and
#                                 re-aligns counts/logcounts to the matched samples.
#   run_voom_lm()              -> ~ hormone. Sex-separate.
#   run_voom_lm_gw()           -> ~ hormone + GW. Sex-separate.
#   run_voom_lm_sex()          -> ~ hormone + Sex. Sex-combined.
#   run_voom_lm_gw_sex()       -> ~ hormone + GW + Sex. Sex-combined.
#   analyze_hormone()          -> per (cell_type x Sex) split; manual gene filter;
#                                 per-group min-sample threshold. Sex-separate.
#   analyze_hormone_combined() -> per cell_type split only; Sex enters as covariate.
#                                 Sex-combined.
#   plot_top_genes_from_subsets(), summarize_sig_genes(), plot_sig_gene_summary()
#                              -> plotting / summary helpers shared by both modes.
#
# Gene filtering:
#   Manual row-sum filter applied per group inside analyze_hormone* functions.
#   Keeps genes with >= count_threshold counts in >= gene_sample_threshold samples.
#   Defaults: count_threshold = 2, gene_sample_threshold = 3.
#
# Sample thresholding:
#   min_samples is a required, caller-supplied absolute floor (no default).
#   Groups with fewer than min_samples non-zero-hormone samples are dropped
#   before model fitting. The floor must exceed the model's coefficient count
#   to leave residual df (e.g. 10 for TT, 21 for E2/P4).


# ---- 1. Pseudobulk aggregation -----------------------------------------------

#' Aggregate a Seurat object to pseudobulk counts by (sample, cell_type).
#'
#' @param seurat_obj  Seurat object.
#' @param sample_col  Name of the sample-ID column in @meta.data (e.g. "Sample").
#' @param celltype_col Name of the cell-type column in @meta.data
#'                    (e.g. "v1_celltype", "celltype", "subclass", "class").
#' @param assay       Assay to pull counts from. Default "RNA".
#' @param layer       Layer/slot to pull. Default "counts".
#'
#' @return A list with:
#'   counts    : matrix (pseudobulk_samples x genes) of summed raw counts
#'   logcounts : matrix (pseudobulk_samples x genes), log1p(counts)
#'   metadata  : data.frame, rownames = pseudobulk_id ("sample.celltype"),
#'               columns: sample_id, cell_type, cell_count
build_pseudobulk <- function(seurat_obj,
                             sample_col,
                             celltype_col,
                             assay = "RNA",
                             layer = "counts") {

  stopifnot(sample_col   %in% colnames(seurat_obj@meta.data))
  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  counts <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
  meta   <- seurat_obj@meta.data

  # Build grouping IDs from the columns the caller chose; do not overwrite anything.
  group_ids    <- interaction(meta[[sample_col]], meta[[celltype_col]], drop = TRUE)
  group_levels <- levels(group_ids)

  aggregated_list <- list()
  metadata_list   <- list()

  for (gid in group_levels) {
    group_cells <- which(group_ids == gid)
    if (length(group_cells) < 1) next
    group_counts <- Matrix::rowSums(counts[, group_cells, drop = FALSE])
    aggregated_list[[gid]] <- group_counts
    metadata_list[[gid]]   <- data.frame(
      sample_id  = meta[[sample_col]][group_cells[1]],
      cell_type  = meta[[celltype_col]][group_cells[1]],
      cell_count = length(group_cells),
      row.names  = gid
    )
  }

  pseudobulk_mat     <- do.call(rbind, aggregated_list)
  metadata_df        <- do.call(rbind, metadata_list)
  pseudobulk_mat_log <- log1p(as.matrix(pseudobulk_mat))

  list(
    counts    = pseudobulk_mat,
    logcounts = pseudobulk_mat_log,
    metadata  = metadata_df
  )
}


# ---- 2. Merge hormone phenotype data and align ------------------------------

#' Join hormone (or other phenotype) metadata onto a pseudobulk list and
#' re-align counts / logcounts to the matched samples.
#'
#' Assumes the hormone data frame has a `Sample` column whose values, when
#' prefixed with `sample_id_prefix`, match the `sample_id` field in pb$metadata.
#'
#' @param pb                Output of build_pseudobulk().
#' @param hormone_df        Data frame with at minimum: Sample, Sex, TT, E2, P4, GW.
#' @param sample_id_prefix  String prepended to hormone_df$Sample to construct the
#'                          join key on the hormone side. Default "tissue".
#' @param pb_sample_pattern Optional regex. If non-NULL, the join key on the
#'                          pseudobulk side is str_extract(sample_id, pb_sample_pattern)
#'                          rather than the full sample_id. Use this when the
#'                          Seurat object stores sample IDs with extra prefixes
#'                          (e.g. "SSD01_tissue1") and you need to extract the
#'                          matching portion (e.g. "tissue1"). Default NULL
#'                          (use sample_id as-is).
#' @param min_cells         Minimum cells per pseudobulk sample to retain. Default 10.
#'
#' @return Updated pb list (counts, logcounts, metadata) restricted to the
#'         pseudobulk samples that survived the join and the min_cells filter.
merge_hormone_metadata <- function(pb,
                                   hormone_df,
                                   sample_id_prefix = "tissue",
                                   pb_sample_pattern = NULL,
                                   min_cells = 10) {

  hormone_df$Sample <- paste0(sample_id_prefix, hormone_df$Sample)

  pb_metadata <- pb$metadata %>%
    as.data.frame() %>%
    tibble::rownames_to_column("pseudobulk_id") %>%
    dplyr::mutate(
      .join_key = if (is.null(pb_sample_pattern)) sample_id
                  else stringr::str_extract(sample_id, pb_sample_pattern)
    ) %>%
    dplyr::left_join(
      hormone_df %>% dplyr::select(Sample, Sex, TT, E2, P4, GW),
      by = c(".join_key" = "Sample")
    ) %>%
    dplyr::select(-.join_key) %>%
    dplyr::filter(!is.na(Sex), cell_count >= min_cells) %>%
    dplyr::mutate(group = paste(cell_type, Sex, sep = "_")) %>%
    tibble::column_to_rownames("pseudobulk_id")

  common_ids <- intersect(rownames(pb$logcounts), rownames(pb_metadata))

  pb$counts    <- pb$counts[common_ids, , drop = FALSE]
  pb$logcounts <- pb$logcounts[common_ids, , drop = FALSE]
  pb$metadata  <- pb_metadata[common_ids, , drop = FALSE]

  pb
}


# ---- 3. limma-voom model wrappers (sex-separate) ----------------------------

# Helper: backtick-quote variable names that are not syntactically valid.
.safe_vars <- function(x) ifelse(make.names(x) != x, paste0("`", x, "`"), x)

run_voom_lm <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  formula_vars <- .safe_vars(c(predictor))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_gw <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"GW" %in% colnames(metadata)) metadata$GW <- NA
  formula_vars <- .safe_vars(c(predictor, "GW"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}


# ---- 3b. limma-voom model wrappers (sex-combined) ---------------------------
#
# Identical to the sex-separate wrappers above, except that Sex is included as
# an additive covariate. Including Sex controls for sex-related expression
# differences when estimating the hormone coefficient across the pooled
# (male + female) sample. See Law et al. 2018 (F1000Research) for the general
# limma-edgeR workflow and the rationale for blocking or covarying on known
# batch/biological factors before testing.

run_voom_lm_sex <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  formula_vars <- .safe_vars(c(predictor, "Sex"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_gw_sex <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"GW" %in% colnames(metadata)) metadata$GW <- NA
  formula_vars <- .safe_vars(c(predictor, "GW", "Sex"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}


# ---- 3c. limma-voom model wrappers (batch-adjusted) -------------------------
#
# Rationale: both batches draw cells from the same biological samples but may
# capture different cell subsets at different depths, producing residual
# variation in pseudobulk library composition that is orthogonal to hormone
# concentration (confirmed by Gene_Covariate_Testing.qmd: |r| <= 0.13).
# Including batch as an additive covariate absorbs this technical variation
# without conflating it with the hormone signal.
#
# Encoding: add_batch_fraction() (see below) computes batch1_frac — the
# fraction of each pseudobulk sample's cells drawn from batch1 — and adds it to
# $metadata as a continuous covariate. Across the current samples this fraction
# is distributed continuously between 0 and 1 (not bimodal), so a continuous
# term captures the batch gradient better than dichotomizing at 0.5, which
# would collapse near-identical intermediate samples into opposite categories
# and discard information (Royston, Altman & Sauerbrei 2006, Stat Med 25:127).
# The coefficient is the expression change per unit increase in batch1
# fraction, holding hormone constant. batch1_frac is a proportion bounded
# [0, 1] and well-defined for every sample (all-batch2 = 0, all-batch1 = 1), so
# it has none of the zero-denominator problem of a batch1/batch2 ratio.
#
# batch1_frac must be present in pb$metadata before calling these wrappers.
# add_batch_fraction() issues a warning for any pseudobulk samples with NA
# batch information; those samples will be silently dropped by model.matrix().
#
# run_voom_lm_gw_batch()   -> ~ hormone + GW + batch1_frac. Sex-separate.
# run_voom_lm_batch()      -> ~ hormone + batch1_frac. Sex-separate.
# run_voom_lm_sex_batch()  -> ~ hormone + Sex + batch1_frac. Sex-combined.

run_voom_lm_gw_batch <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"GW" %in% colnames(metadata)) metadata$GW <- NA
  if (!"batch1_frac" %in% colnames(metadata))
    stop("batch1_frac not found in metadata. Run add_batch_fraction() first.")
  formula_vars <- .safe_vars(c(predictor, "GW", "batch1_frac"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_batch <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"batch1_frac" %in% colnames(metadata))
    stop("batch1_frac not found in metadata. Run add_batch_fraction() first.")
  formula_vars <- .safe_vars(c(predictor, "batch1_frac"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_sex_batch <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"batch1_frac" %in% colnames(metadata))
    stop("batch1_frac not found in metadata. Run add_batch_fraction() first.")
  formula_vars <- .safe_vars(c(predictor, "Sex", "batch1_frac"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}


#' Compute per-sample batch1 fraction, then join onto pb$metadata.
#'
#' Call this AFTER merge_hormone_metadata() so the join operates on the
#' already-filtered sample set. One column is added to $metadata:
#'   batch1_frac : numeric, fraction of batch1 cells. Continuous batch
#'                 covariate used by run_voom_lm_gw_batch(),
#'                 run_voom_lm_batch(), and run_voom_lm_sex_batch().
#'
#' @param pb           Output of merge_hormone_metadata(). Must contain
#'                     $metadata with a sample_id column.
#' @param seurat_obj   Seurat object whose @meta.data contains sample_col
#'                     and batch_col.
#' @param sample_col   Name of the sample ID column in @meta.data. Default "Sample".
#' @param batch_col    Name of the batch column in @meta.data. Default "batch".
#' @param batch1_label Value in batch_col representing batch 1. Default "batch1".
#'
#' @return Updated pb list with batch1_frac added to $metadata.
add_batch_fraction <- function(pb,
                               seurat_obj,
                               sample_col   = "Sample",
                               batch_col    = "batch",
                               batch1_label = "batch1") {

  batch_props <- seurat_obj@meta.data %>%
    dplyr::group_by(.data[[sample_col]]) %>%
    dplyr::summarise(
      batch1_frac = mean(.data[[batch_col]] == batch1_label, na.rm = TRUE),
      .groups     = "drop"
    )

  pb$metadata <- pb$metadata %>%
    tibble::rownames_to_column("pseudobulk_id") %>%
    dplyr::left_join(batch_props, by = c("sample_id" = sample_col)) %>%
    tibble::column_to_rownames("pseudobulk_id")

  n_missing <- sum(is.na(pb$metadata$batch1_frac))
  if (n_missing > 0)
    warning(sprintf(
      "add_batch_fraction: %d pseudobulk sample(s) have no batch info (NA batch1_frac).",
      n_missing
    ))

  pb
}


# ---- 4. Per-group analysis (cell_type x Sex) -- sex-separate -----------------
#
# Gene filter: keeps genes with >= count_threshold counts in >=
#   gene_sample_threshold samples within each (cell_type x Sex) group.
#   Defaults (count_threshold = 2, gene_sample_threshold = 3) match the manual
#   filter used across the project Rmd scripts.
#
# Sample threshold: min_samples is a fixed absolute floor supplied by the
#   caller (required, no default). A group is dropped unless at least
#   min_samples samples have a non-zero hormone value. The floor must exceed
#   the number of model coefficients or the fit has no residual df; callers set
#   it per model (e.g. 10 for the male-only TT models). A fixed floor keeps the
#   inclusion rule consistent across cell types rather than scaling it to each
#   type's sample count.

analyze_hormone <- function(pb_counts, metadata_merged, hormone,
                            min_samples,
                            count_threshold       = 2,
                            gene_sample_threshold = 3,
                            model                 = run_voom_lm) {

  if (missing(min_samples))
    stop("min_samples must be specified (absolute minimum samples per group).")

  split_ids <- split(rownames(metadata_merged),
                     interaction(metadata_merged$cell_type, metadata_merged$Sex,
                                 drop = TRUE))

  results_by_group <- purrr::map(names(split_ids), function(group_name) {
    ids          <- split_ids[[group_name]]
    counts_sub   <- t(pb_counts)[, ids, drop = FALSE]
    metadata_sub <- metadata_merged[ids, , drop = FALSE]

    overall_min_samples <- min_samples

    cat("\nGroup:", group_name, "\n")
    cat("Initial number of samples:", nrow(metadata_sub), "\n")
    cat("Min samples required (fixed floor):", overall_min_samples, "\n")

    keep_genes <- rowSums(counts_sub >= count_threshold) >= gene_sample_threshold
    counts_sub <- counts_sub[keep_genes, , drop = FALSE]
    cat("Genes passing manual filter:", sum(keep_genes), "\n")

    if (nrow(counts_sub) == 0) {
      cat("Group removed: no genes passed the filter\n")
      return(NULL)
    }

    nonzero_hormone_mask <- metadata_sub[[hormone]] != 0 &
                            !is.na(metadata_sub[[hormone]])
    cat("Samples with non-zero hormone:", sum(nonzero_hormone_mask), "\n")

    if (sum(nonzero_hormone_mask) < overall_min_samples) {
      cat("Group removed: insufficient non-zero hormone values",
          sum(nonzero_hormone_mask), "\n")
      return(NULL)
    }

    final_counts_sub   <- counts_sub[, nonzero_hormone_mask, drop = FALSE]
    final_metadata_sub <- metadata_sub[nonzero_hormone_mask, , drop = FALSE]

    cat("Final sample size:", nrow(final_metadata_sub), "\n")

    result <- model(final_counts_sub, final_metadata_sub, hormone)
    if (is.null(result)) {
      cat("Model returned NULL result\n")
      return(NULL)
    }

    result$cell_type <- final_metadata_sub$cell_type[1]
    result$sex       <- final_metadata_sub$Sex[1]

    list(
      limma_results = result,
      counts_sub    = final_counts_sub,
      metadata_sub  = final_metadata_sub
    )
  })

  results_by_group <- results_by_group[!sapply(results_by_group, is.null)]
  combined_results <- purrr::map_dfr(results_by_group, ~ .x$limma_results)
  if (nrow(combined_results) == 0) cat("No significant results found\n")

  list(combined_results = combined_results, by_group_data = results_by_group)
}


# ---- 4b. Per-group analysis (cell_type only) -- sex-combined -----------------
#
# Mirrors analyze_hormone() with two key differences:
#   1. Splits only by cell_type (not cell_type x Sex); Sex enters the model
#      as a covariate via the run_voom_lm_*_sex() wrappers above. This is the
#      standard approach for controlling a known biological factor when pooling
#      groups -- see Ritchie et al. 2015 (limma paper) and Law et al. 2018.
#   2. result$sex is set to "combined" rather than a single-sex value.
#
# Sex guard: if fewer than 2 sexes are present after hormone filtering, the
# group is skipped -- the Sex covariate would be rank-deficient in the design.
#
# Sample threshold: min_samples is a fixed absolute floor supplied by the
#   caller (required, no default). A group is dropped unless at least
#   min_samples samples have a non-zero hormone value. It must exceed the
#   number of model coefficients or the fit has no residual df; callers set it
#   per model (e.g. 21 for the sex-combined E2/P4 models).

analyze_hormone_combined <- function(pb_counts, metadata_merged, hormone,
                                     min_samples,
                                     count_threshold       = 2,
                                     gene_sample_threshold = 3,
                                     model                 = run_voom_lm_sex) {

  if (missing(min_samples))
    stop("min_samples must be specified (absolute minimum samples per group).")

  split_ids <- split(rownames(metadata_merged),
                     metadata_merged$cell_type)   # cell_type only; no Sex split

  results_by_group <- purrr::map(names(split_ids), function(group_name) {
    ids          <- split_ids[[group_name]]
    counts_sub   <- t(pb_counts)[, ids, drop = FALSE]
    metadata_sub <- metadata_merged[ids, , drop = FALSE]

    overall_min_samples <- min_samples

    cat("\nGroup:", group_name, "\n")
    cat("Initial number of samples:", nrow(metadata_sub), "\n")
    cat("Min samples required (fixed floor):", overall_min_samples, "\n")

    keep_genes <- rowSums(counts_sub >= count_threshold) >= gene_sample_threshold
    counts_sub <- counts_sub[keep_genes, , drop = FALSE]
    cat("Genes passing manual filter:", sum(keep_genes), "\n")

    if (nrow(counts_sub) == 0) {
      cat("Group removed: no genes passed the filter\n")
      return(NULL)
    }

    nonzero_hormone_mask <- metadata_sub[[hormone]] != 0 &
                            !is.na(metadata_sub[[hormone]])
    cat("Samples with non-zero hormone:", sum(nonzero_hormone_mask), "\n")

    if (sum(nonzero_hormone_mask) < overall_min_samples) {
      cat("Group removed: insufficient non-zero hormone values",
          sum(nonzero_hormone_mask), "\n")
      return(NULL)
    }

    final_counts_sub   <- counts_sub[, nonzero_hormone_mask, drop = FALSE]
    final_metadata_sub <- metadata_sub[nonzero_hormone_mask, , drop = FALSE]

    if (length(unique(final_metadata_sub$Sex)) < 2) {
      cat("Group removed: only one sex present after hormone filtering",
          "(Sex covariate would be rank-deficient)\n")
      return(NULL)
    }

    cat("Final sample size:", nrow(final_metadata_sub), "\n")

    result <- model(final_counts_sub, final_metadata_sub, hormone)
    if (is.null(result)) {
      cat("Model returned NULL result\n")
      return(NULL)
    }

    result$cell_type <- final_metadata_sub$cell_type[1]
    result$sex       <- "combined"

    list(
      limma_results = result,
      counts_sub    = final_counts_sub,
      metadata_sub  = final_metadata_sub
    )
  })

  results_by_group <- results_by_group[!sapply(results_by_group, is.null)]
  combined_results <- purrr::map_dfr(results_by_group, ~ .x$limma_results)
  if (nrow(combined_results) == 0) cat("No significant results found\n")

  list(combined_results = combined_results, by_group_data = results_by_group)
}


# ---- 5. Plotting & summary helpers ------------------------------------------

#' Plot expression-vs-hormone trajectories for significant genes.
#'
#' @param color_by  Point colour aesthetic. "cell_count" (default, viridis)
#'   preserves the original behaviour; "GW" uses the gestational-week gradient
#'   shared with the sample PCA plots (#C8E6C9 -> #1B5E20).
#' @param single_fit  If TRUE, draw ONE `lm` fit across the whole subset instead
#'   of one per sex. Set TRUE for sex-combined models (e.g. ~ hormone + GW + Sex),
#'   where the model estimates a single hormone slope and sex enters only as an
#'   intercept offset -- a per-sex fit would misrepresent the tested effect.
#'   Sex is still encoded by point shape.
plot_top_genes_from_subsets <- function(by_group_data, output_dir, hormone,
                                        fdr_cutoff = 0.1,
                                        logfc_cutoff = 0.25,
                                        top_n = 15,
                                        color_by = c("cell_count", "GW"),
                                        single_fit = FALSE) {
  color_by <- match.arg(color_by)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("Number of groups to plot: ", length(by_group_data))

  purrr::walk(by_group_data, function(group_list) {
    if (is.null(group_list)) {
      message("- Skipping NULL group")
      return()
    }

    limma_res    <- group_list$limma_results
    counts_sub   <- group_list$counts_sub        # genes x samples
    metadata_sub <- group_list$metadata_sub

    if (nrow(metadata_sub) < 1 || ncol(counts_sub) < 1) {
      message("- Skipping empty subset: ",
              metadata_sub$cell_type[1], "_", metadata_sub$Sex[1])
      return()
    }

    # Fall back to cell_count if GW was requested but is absent/all-NA for this
    # subset (run_voom_* fills GW with NA when it is missing from metadata).
    color_var <- color_by
    if (color_by == "GW" &&
        (!"GW" %in% colnames(metadata_sub) || all(is.na(metadata_sub$GW)))) {
      message("- GW unavailable for this subset; colouring by cell_count")
      color_var <- "cell_count"
    }

    counts_cpm <- edgeR::cpm(counts_sub, log = FALSE)
    counts_log <- t(log1p(counts_cpm))           # samples x genes

    # Apply the same FDR + logFC filter as summarize_sig_genes(), so the number
    # of trajectories plotted matches the bar-chart count for this group.
    # top_n is a safety cap (defaults to Inf = no cap).
    sig_genes <- limma_res %>% dplyr::filter(adj.P.Val < fdr_cutoff)
    if (!is.null(logfc_cutoff)) {
      sig_genes <- sig_genes %>% dplyr::filter(abs(logFC) >= logfc_cutoff)
    }
    sig_genes <- sig_genes %>%
      dplyr::arrange(-abs(t)) %>%
      dplyr::slice_head(n = top_n)

    top_genes       <- sig_genes %>% dplyr::pull(gene)
    top_genes_pvals <- sig_genes %>%
      dplyr::select(gene, adj.P.Val) %>%
      as.data.frame()

    # Determine sex label for file naming: auto-detect "combined" when both
    # sexes are present (sex-combined analysis), otherwise use the single value.
    label_sex  <- if (length(unique(metadata_sub$Sex)) == 1) metadata_sub$Sex[1] else "combined"

    message(" Plotting group: ", metadata_sub$cell_type[1], "_", label_sex)
    message("- Significant genes (FDR < ", fdr_cutoff,
            if (!is.null(logfc_cutoff)) paste0(", |logFC| >= ", logfc_cutoff) else "",
            "): ", length(top_genes))

    plot_list <- purrr::map(seq_len(nrow(top_genes_pvals)), function(i) {
      gene     <- top_genes_pvals$gene[i]
      adj_pval <- top_genes_pvals$adj.P.Val[i]

      if (!gene %in% colnames(counts_log)) {
        message("    Gene ", gene, " not in counts -- skipping")
        return(NULL)
      }

      expr_vals <- counts_log[, gene]
      if (all(is.na(expr_vals)) || all(expr_vals == 0)) {
        message("    Gene ", gene, " has all NA or zero expression -- skipping")
        return(NULL)
      }

      message("    Plotting gene: ", gene,
              " (adj. p-value: ", format(adj_pval, digits = 3), ")")
      flush.console()

      df <- metadata_sub %>% dplyr::mutate(expr = expr_vals)

      ggplot(df, aes(x = !!sym(hormone), y = expr,
                     color = !!sym(color_var), shape = Sex)) +
        geom_point(size = 2, alpha = 0.8) +
        # aes(group = 1) overrides the grouping that `shape = Sex` would
        # otherwise pass to geom_smooth; colour/fill are fixed so the fit is not
        # split by the continuous colour aesthetic either.
        (if (single_fit) {
          geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                      color = "grey25", fill = "grey70")
        } else {
          geom_smooth(method = "lm", se = TRUE)
        }) +
        (if (color_var == "GW") {
          scale_color_gradient(name = "Gestational Week",
                               low = "#C8E6C9", high = "#1B5E20")
        } else {
          scale_color_viridis_c()
        }) +
        scale_shape_manual(values = c("M" = 16, "F" = 17)) +
        labs(title = paste0(gene,
                            " (adj. p-value: ", format(adj_pval, digits = 3),
                            ") in ", metadata_sub$cell_type[1]),
             x = hormone, y = "log1p(CPM)") +
        theme_minimal(base_size = 14)
    }) %>% purrr::compact()

    if (length(plot_list) == 0) {
      message("No valid plots for group")
      return()
    }

    combined_plot <- patchwork::wrap_plots(plotlist = plot_list, ncol = 2)
    clean_name    <- gsub("/", "-",
                          paste0(metadata_sub$cell_type[1], "_", label_sex))

    save_dual_format(combined_plot,
                     output_dir,
                     paste0(clean_name, "_top_genes"),
                     width  = 10,
                     height = 6 + ceiling(length(plot_list) / 2))
    message("  Saved plot for ", clean_name)
  })
}

summarize_sig_genes <- function(results, fdr_cutoff = 0.1, logfc_cutoff = 0.25) {
  sig <- results %>% dplyr::filter(adj.P.Val < fdr_cutoff)
  if (!is.null(logfc_cutoff)) {
    sig <- sig %>% dplyr::filter(abs(logFC) >= logfc_cutoff)
  }
  sig %>%
    dplyr::count(cell_type, sex, name = "num_sig_genes") %>%
    dplyr::arrange(desc(num_sig_genes))
}

plot_sig_gene_summary <- function(summary_df, title) {
  ggplot(summary_df, aes(x = cell_type, y = num_sig_genes, fill = sex)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = title, x = "Cell Type", y = "# Significant Genes") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # "combined" uses medium grey (#888888): neutral alongside steelblue (M) and
    # violetred3 (F); grey is distinguishable across deuteranopia and protanopia.
    scale_fill_manual(values = c("M" = "steelblue", "F" = "violetred3",
                                 "combined" = "#888888"))
}


# ---- 6. Covariate regression utilities --------------------------------------
#
# Functions for testing whether sample-level covariates predict pseudobulk
# gene expression at the whole-transcriptome level. Used by
# Gene_Covariate_Testing.qmd to justify covariate inclusion in limma-voom
# hormone ~ gene models.
#
# Design:
#   perform_regression()         -> gene-wise lm for any covariate (continuous
#                                   or categorical). Uses lapply + do.call(rbind)
#                                   to avoid the O(n^2) cost of growing a data
#                                   frame row-by-row in a loop. Returns tidy
#                                   results with BH-adjusted p-values.
#   summarize_regression_results() -> prints counts and top-10 table.
#   plot_covariate_summary()     -> 4-panel ggplot2 diagnostic: volcano,
#                                   effect-size histogram, p-value histogram,
#                                   R² histogram. Handles empty results
#                                   gracefully (no base-R xlim crash).
#   plot_top_genes_continuous()  -> scatter + lm ribbon for top n genes.
#   plot_top_genes_categorical() -> violin + jitter for top n genes.
#   save_dual_format()           -> saves a ggplot as PNG (output_dir/) and
#                                   SVG (output_dir/svg/).
#
# Input:
#   All perform_regression / plot_top_genes functions expect the pb_object
#   returned by merge_hormone_metadata() -- i.e. a list with:
#     $logcounts  matrix (pseudobulk_samples x genes), log1p(counts)
#     $metadata   data.frame, one row per pseudobulk sample


#' Run gene-wise linear regression against a single covariate.
#'
#' @param pb_object       Output of merge_hormone_metadata(). Must contain
#'                        $logcounts (samples x genes) and $metadata.
#' @param covariate       Column name in pb_object$metadata to test.
#' @param covariate_label Human-readable label for messages and plots.
#'                        Defaults to covariate.
#' @param analysis_name   String used in progress messages.
#'
#' @details
#'   Continuous covariates (is.numeric): fits lm(expr ~ cov).
#'   Categorical covariates (character/factor): fits lm(expr ~ factor(cov)).
#'   The reported estimate is the second factor level vs. the reference
#'   (first level alphabetically). Genes with < 3 valid observations are
#'   skipped. If the model matrix collapses to a single column (e.g. only
#'   one group present after NA removal), the gene is skipped via the
#'   nrow(ct) < 2 guard.
#'
#' @return Data frame with columns: gene, estimate, std_error, t_value,
#'   p_value, r_squared, n_samples, adj_p_value. Ordered by adj_p_value.
perform_regression <- function(pb_object,
                               covariate,
                               covariate_label = NULL,
                               analysis_name   = "Pseudobulk") {

  if (is.null(covariate_label)) covariate_label <- covariate

  logcounts <- pb_object$logcounts
  metadata  <- pb_object$metadata

  if (!covariate %in% colnames(metadata))
    stop(sprintf(
      "Column '%s' not found in metadata. Available: %s",
      covariate, paste(colnames(metadata), collapse = ", ")
    ))

  cov_vals      <- metadata[[covariate]]
  is_continuous <- is.numeric(cov_vals)

  if (!is_continuous) {
    cov_vals  <- factor(cov_vals)
    ref_level <- levels(cov_vals)[1]
    cat(sprintf(
      "Categorical covariate '%s' -- reference level: '%s'\n",
      covariate_label, ref_level
    ))
  }

  genes   <- colnames(logcounts)
  n_genes <- length(genes)
  cat(sprintf(
    "Running regression: %d genes ~ '%s' (%s) in %s\n",
    n_genes, covariate_label,
    ifelse(is_continuous, "continuous", "categorical"),
    analysis_name
  ))

  results_list <- lapply(seq_along(genes), function(i) {
    if (i %% 1000 == 0)
      cat(sprintf("  %d / %d (%.0f%%)\n", i, n_genes, 100 * i / n_genes))

    expr  <- logcounts[, genes[i]]
    valid <- !is.na(expr) & !is.na(cov_vals)
    if (sum(valid) < 3) return(NULL)

    tryCatch({
      m  <- lm(expr[valid] ~ cov_vals[valid])
      ms <- summary(m)
      ct <- coef(ms)
      if (nrow(ct) < 2) return(NULL)
      data.frame(
        gene          = genes[i],
        estimate      = ct[2L, "Estimate"],
        std_error     = ct[2L, "Std. Error"],
        t_value       = ct[2L, "t value"],
        p_value       = ct[2L, "Pr(>|t|)"],
        r_squared     = ms$r.squared,
        n_samples     = sum(valid),
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)
  })

  results <- do.call(rbind, Filter(Negate(is.null), results_list))

  if (!is.null(results) && nrow(results) > 0) {
    results$adj_p_value <- p.adjust(results$p_value, method = "fdr")
    results             <- results[order(results$adj_p_value), ]
    rownames(results)   <- NULL
  } else {
    warning(sprintf("No valid regression results for '%s'.", covariate_label))
    results <- data.frame(
      gene = character(), estimate = numeric(), std_error = numeric(),
      t_value = numeric(), p_value = numeric(), r_squared = numeric(),
      n_samples = integer(), adj_p_value = numeric()
    )
  }

  cat(sprintf(
    "Done. %d / %d genes returned valid results.\n",
    nrow(results), n_genes
  ))
  results
}


#' Print a text summary of perform_regression() output.
#'
#' @param results         Output of perform_regression().
#' @param covariate_label Label printed in the header.
#' @param alpha           FDR threshold for significance counts.
summarize_regression_results <- function(results, covariate_label, alpha = 0.05) {
  n_sig <- sum(results$adj_p_value < alpha, na.rm = TRUE)
  n_pos <- sum(results$adj_p_value < alpha & results$estimate > 0, na.rm = TRUE)
  n_neg <- sum(results$adj_p_value < alpha & results$estimate < 0, na.rm = TRUE)

  cat(sprintf("\n=== %s ===\n",         toupper(covariate_label)))
  cat(sprintf("  Genes tested      : %d\n",     nrow(results)))
  cat(sprintf("  Sig. FDR < %.2f  : %d (%.1f%%)\n",
              alpha, n_sig, 100 * n_sig / max(nrow(results), 1)))
  cat(sprintf("  Positive effect   : %d\n",     n_pos))
  cat(sprintf("  Negative effect   : %d\n\n",   n_neg))
  cat("Top 10:\n")
  print(
    head(results[, c("gene", "estimate", "p_value", "adj_p_value", "r_squared")], 10),
    row.names = FALSE, digits = 4
  )
}


#' ggplot2 4-panel diagnostic plot for perform_regression() output.
#'
#' Panels: volcano, effect-size histogram, p-value histogram, R² histogram.
#' Returns the combined patchwork plot invisibly if results are empty.
#'
#' @param results         Output of perform_regression().
#' @param covariate_label Character label for axis and plot titles.
#' @param alpha           FDR threshold for significance colouring.
#' @param output_dir      If non-NULL (with filename_base), saves via
#'                        save_dual_format().
#' @param filename_base   File stem without extension.
#' @param width,height    Dimensions in inches passed to save_dual_format().
#'
#' @return A patchwork ggplot object.
plot_covariate_summary <- function(results,
                                   covariate_label,
                                   alpha         = 0.05,
                                   output_dir    = NULL,
                                   filename_base = NULL,
                                   width         = 10,
                                   height        = 7) {

  if (nrow(results) == 0) {
    warning("No results to plot for: ", covariate_label)
    return(invisible(NULL))
  }

  res    <- results %>% dplyr::mutate(significant = adj_p_value < alpha)
  n_sig  <- sum(res$significant,   na.rm = TRUE)
  n_tot  <- nrow(res)
  med_r2 <- median(res$r_squared,  na.rm = TRUE)
  mu_eff <- mean(res$estimate,     na.rm = TRUE)

  p_volcano <- ggplot(res, aes(x = estimate,
                                y = -log10(p_value),
                                color = significant)) +
    geom_point(alpha = 0.35, size = 0.7) +
    scale_color_manual(
      values = c("FALSE" = "grey60", "TRUE" = "#E74C3C"),
      labels = c("FALSE" = "Non-sig.",
                 "TRUE"  = sprintf("FDR < %.2f", alpha))
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "steelblue", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.3) +
    labs(title = "Volcano",
         x     = sprintf("%s coefficient", covariate_label),
         y     = expression(-log[10](p)),
         color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  p_effect <- ggplot(res, aes(x = estimate)) +
    geom_histogram(bins = 60, fill = "steelblue", color = NA, alpha = 0.85) +
    geom_vline(xintercept = 0,      linetype = "dashed",
               color = "#E74C3C",  linewidth = 0.6) +
    geom_vline(xintercept = mu_eff, linetype = "dashed",
               color = "navy",     linewidth = 0.6) +
    annotate("text", x = Inf, y = Inf,
             label  = sprintf("Mean = %.4f", mu_eff),
             hjust  = 1.1, vjust = 1.5, size = 3.2) +
    labs(title = "Effect Size Distribution",
         x     = sprintf("%s coefficient", covariate_label),
         y     = "Genes") +
    theme_minimal(base_size = 11)

  p_pval <- ggplot(res, aes(x = p_value)) +
    geom_histogram(bins = 50, fill = "#2ECC71", color = NA, alpha = 0.85) +
    geom_vline(xintercept = 0.05, linetype = "dashed",
               color = "#E74C3C", linewidth = 0.6) +
    labs(title = "P-value Distribution", x = "P-value", y = "Genes") +
    theme_minimal(base_size = 11)

  p_r2 <- ggplot(res, aes(x = r_squared)) +
    geom_histogram(bins = 50, fill = "#F39C12", color = NA, alpha = 0.85) +
    geom_vline(xintercept = med_r2, linetype = "dashed",
               color = "darkgreen", linewidth = 0.6) +
    annotate("text", x = Inf, y = Inf,
             label  = sprintf("Median R\u00b2 = %.3f", med_r2),
             hjust  = 1.1, vjust = 1.5, size = 3.2) +
    labs(title = expression(R^2 ~ "Distribution"),
         x     = expression(R^2), y = "Genes") +
    theme_minimal(base_size = 11)

  combined <- (p_volcano + p_effect) / (p_pval + p_r2) +
    patchwork::plot_annotation(
      title    = sprintf(
        "%s Effect on Pseudobulk Gene Expression", covariate_label
      ),
      subtitle = sprintf(
        "n\u2009=\u2009%d genes; %d (%.1f%%) significant at FDR\u2009<\u2009%.2f",
        n_tot, n_sig, 100 * n_sig / n_tot, alpha
      )
    )

  if (!is.null(output_dir) && !is.null(filename_base))
    save_dual_format(combined, output_dir, filename_base, width, height)

  combined
}


#' Scatter + lm ribbon plots for the top n genes (continuous covariate).
#'
#' @param pb_object       Output of merge_hormone_metadata().
#' @param results         Output of perform_regression().
#' @param covariate       Column name in pb_object$metadata (must be numeric).
#' @param covariate_label Axis / title label.
#' @param n_genes         Number of top genes to plot.
#' @param output_dir,filename_base,width,height  Passed to save_dual_format().
#'
#' @return A ggplot object.
plot_top_genes_continuous <- function(pb_object,
                                      results,
                                      covariate,
                                      covariate_label,
                                      n_genes       = 6,
                                      output_dir    = NULL,
                                      filename_base = NULL,
                                      width         = 12,
                                      height        = 8) {

  logcounts <- pb_object$logcounts
  metadata  <- pb_object$metadata

  top_genes <- head(
    results$gene[results$gene %in% colnames(logcounts)], n_genes
  )

  plot_data <- purrr::map_dfr(top_genes, function(g) {
    expr  <- logcounts[, g]
    cov   <- as.numeric(metadata[[covariate]])
    valid <- !is.na(expr) & !is.na(cov)
    gr    <- results[results$gene == g, ][1L, ]
    data.frame(
      gene          = g,
      expr          = expr[valid],
      covariate_val = cov[valid],
      facet_label   = sprintf(
        "%s\nFDR\u2009=\u2009%.2e  R\u00b2\u2009=\u2009%.3f",
        g, gr$adj_p_value, gr$r_squared
      ),
      stringsAsFactors = FALSE
    )
  })

  p <- ggplot(plot_data, aes(x = covariate_val, y = expr)) +
    geom_point(alpha = 0.7, color = "steelblue", size = 1.8) +
    geom_smooth(method = "lm", se = TRUE,
                color  = "#E74C3C", linewidth = 0.8,
                fill   = "#E74C3C", alpha = 0.15) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 3) +
    labs(title = sprintf("Top %d Genes \u2014 %s", n_genes, covariate_label),
         x     = covariate_label,
         y     = "log1p expression") +
    theme_minimal(base_size = 11)

  if (!is.null(output_dir) && !is.null(filename_base))
    save_dual_format(p, output_dir, filename_base, width, height)

  p
}


#' Violin + jitter plots for the top n genes (categorical covariate).
#'
#' @param pb_object       Output of merge_hormone_metadata().
#' @param results         Output of perform_regression().
#' @param covariate       Column name in pb_object$metadata (character/factor).
#' @param covariate_label Axis / title label.
#' @param n_genes         Number of top genes to plot.
#' @param output_dir,filename_base,width,height  Passed to save_dual_format().
#'
#' @return A ggplot object.
plot_top_genes_categorical <- function(pb_object,
                                       results,
                                       covariate,
                                       covariate_label,
                                       n_genes       = 6,
                                       output_dir    = NULL,
                                       filename_base = NULL,
                                       width         = 12,
                                       height        = 8) {

  logcounts <- pb_object$logcounts
  metadata  <- pb_object$metadata

  top_genes <- head(
    results$gene[results$gene %in% colnames(logcounts)], n_genes
  )

  plot_data <- purrr::map_dfr(top_genes, function(g) {
    expr  <- logcounts[, g]
    cov   <- factor(metadata[[covariate]])
    valid <- !is.na(expr) & !is.na(cov)
    gr    <- results[results$gene == g, ][1L, ]
    data.frame(
      gene        = g,
      expr        = expr[valid],
      group       = cov[valid],
      facet_label = sprintf(
        "%s\nFDR\u2009=\u2009%.2e  R\u00b2\u2009=\u2009%.3f",
        g, gr$adj_p_value, gr$r_squared
      ),
      stringsAsFactors = FALSE
    )
  })

  p <- ggplot(plot_data, aes(x = group, y = expr, fill = group)) +
    geom_violin(alpha = 0.4, draw_quantiles = 0.5, trim = TRUE) +
    geom_jitter(width = 0.12, alpha = 0.75, size = 1.5, color = "grey30") +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 3) +
    labs(title = sprintf("Top %d Genes \u2014 %s", n_genes, covariate_label),
         x     = covariate_label,
         y     = "log1p expression") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  if (!is.null(output_dir) && !is.null(filename_base))
    save_dual_format(p, output_dir, filename_base, width, height)

  p
}


#' Save a ggplot as PNG and SVG (SVG written to output_dir/svg/).
#'
#' Mirrors the save pattern used in other pipeline scripts. Both files share
#' the same dimensions. The svg/ subdirectory is created if absent.
#'
#' @param plot         A ggplot object.
#' @param output_dir   Directory for the PNG file.
#' @param filename_base File stem without extension.
#' @param width,height Plot dimensions in inches.
#' @param dpi          PNG resolution (default 300).
#'
#' @return The plot, invisibly.
save_dual_format <- function(plot,
                             output_dir,
                             filename_base,
                             width  = 10,
                             height = 7,
                             dpi    = 300) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  svg_dir <- file.path(output_dir, "svg")
  dir.create(svg_dir, showWarnings = FALSE, recursive = TRUE)

  png_path <- file.path(output_dir, paste0(filename_base, ".png"))
  svg_path <- file.path(svg_dir,    paste0(filename_base, ".svg"))

  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  dpi = dpi)
  ggplot2::ggsave(svg_path, plot = plot, width = width, height = height,
                  device = svglite::svglite)

  message("Saved: ", basename(png_path),
          "  +  svg/", basename(svg_path))
  invisible(plot)
}


# ---- 7. Analysis UMAP -------------------------------------------------------
#
# Generates a UMAP showing which cells were included in a hormone-gene analysis.
# Cells from groups that passed all thresholds (outlier removal, min-cells,
# gene filter, non-zero hormone mask, and — for sex-combined — the sex guard)
# are colored by cell type, matching the color scheme of the QC UMAP.
# All other cells are drawn in grey underneath the colored layer.
#
# Cell type labels are placed at cluster centroids computed over all cells
# (both included and excluded), so label positions exactly match the total UMAP.
# Included cell type labels are drawn in their cell type color; excluded cell
# type labels are drawn in grey.
#
# Color matching:
#   A silent reference DimPlot is built with the same group.by and reduction as
#   the QC UMAP; ggplot_build() extracts the auto-assigned colors so the
#   analysis UMAP is pixel-consistent with the QC celltype UMAP without
#   requiring any hard-coded palette.
#
# Sex handling:
#   sex_separate = TRUE  (testosterone, analyze_hormone):
#     Only XY cells from (sample x cell_type) groups where Sex == male_hormone
#     are colored. XX cells — even for included cell types — are grey, because
#     they do not contribute to the testosterone model.
#   sex_separate = FALSE (E2/P4, analyze_hormone_combined):
#     All cells from (sample x cell_type) groups that ran the model are colored,
#     regardless of sex.


#' Plot a UMAP highlighting cells included in the hormone-gene analysis.
#'
#' @param seurat_obj     Seurat object. Must contain the reduction specified
#'                       by `reduction` and the columns `celltype_col`,
#'                       `sample_col`, and (for sex_separate) `seurat_sex_col`.
#' @param by_group_data  The `by_group_data` element from `analyze_hormone()`
#'                       or `analyze_hormone_combined()`. Each non-NULL element
#'                       must have a `metadata_sub` data frame with columns
#'                       sample_id, cell_type, and Sex.
#' @param title          Plot title string.
#' @param sex_separate   Logical. TRUE for sex-separate analysis (testosterone):
#'                       only XY cells from included (sample x cell_type) groups
#'                       are colored; XX cells are grey even for included cell
#'                       types. FALSE for sex-combined analysis (E2, P4): all
#'                       cells from included groups are colored.
#' @param celltype_col   Column in Seurat @meta.data holding cell type labels.
#'                       Default "celltype".
#' @param sample_col     Column in Seurat @meta.data holding sample IDs.
#'                       Default "Sample". Values must match metadata_sub$sample_id.
#' @param seurat_sex_col Column in Seurat @meta.data holding sex genotype.
#'                       Default "genotype" (values "XX" / "XY").
#' @param male_genotype  Value of seurat_sex_col corresponding to male.
#'                       Default "XY".
#' @param male_hormone   Value of Sex in metadata_sub corresponding to male.
#'                       Default "M".
#' @param excluded_color Color for excluded cell points. Default "grey80".
#' @param label_excl_color Color for excluded cell type text labels.
#'                       Default "grey50".
#' @param reduction      Seurat reduction to use. Default "umap.rna".
#' @param pt_size        Point size passed to geom_point. Default 0.3.
#' @param label_size     Text size for cell type labels. Default 3.
#' @param output_dir,filename_base,width,height
#'                       Passed to save_dual_format(). If either is NULL,
#'                       no file is saved.
#'
#' @return A ggplot object, invisibly.
plot_analysis_umap <- function(seurat_obj,
                               by_group_data,
                               title,
                               sex_separate     = FALSE,
                               celltype_col     = "celltype",
                               sample_col       = "Sample",
                               seurat_sex_col   = "genotype",
                               male_genotype    = "XY",
                               male_hormone     = "M",
                               excluded_color   = "grey80",
                               label_excl_color = "grey50",
                               reduction        = "umap.rna",
                               pt_size          = 0.3,
                               label_size       = 3,
                               output_dir       = NULL,
                               filename_base    = NULL,
                               width            = 10,
                               height           = 8) {

  # ---- 1. Extract reference cell type colors (match QC UMAP) ---------------
  # Build a silent DimPlot with the same grouping as the QC UMAP, then use
  # ggplot_build() to retrieve the auto-assigned color per factor level.
  # This guarantees color consistency without hard-coding any palette.
  p_ref    <- suppressMessages(
    DimPlot(seurat_obj, reduction = reduction, group.by = celltype_col)
  )
  ref_data <- ggplot_build(p_ref)$data[[1]]

  # Factor levels in the order Seurat uses internally (alphabetical for
  # character columns; existing levels for factor columns).
  ct_levels <- levels(factor(seurat_obj@meta.data[[celltype_col]]))

  # ggplot_build group integers are 1-based and match the factor level order.
  ref_map <- ref_data %>%
    dplyr::distinct(group, colour) %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::arrange(group) %>%
    dplyr::mutate(cell_type = ct_levels[group])

  ct_color_map <- setNames(ref_map$colour, ref_map$cell_type)

  # ---- 2. Identify included (sample x cell_type) pairs ----------------------
  meta <- seurat_obj@meta.data

  # Unpack every non-NULL group from by_group_data into a flat table of
  # (sample_id, cell_type, sex_meta) rows — one row per sample per group.
  included_df <- purrr::map_dfr(by_group_data, function(g) {
    if (is.null(g)) return(NULL)
    ms <- g$metadata_sub
    data.frame(
      sample_id = ms$sample_id,
      cell_type = ms$cell_type[1],          # constant within a group
      sex_meta  = if ("Sex" %in% colnames(ms)) ms$Sex[1] else NA_character_,
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::distinct()

  # ---- 3. Compute per-cell inclusion flag -----------------------------------
  if (sex_separate) {
    # Testosterone: retain only the XY groups (sex_meta == male_hormone) and
    # further require the Seurat cell to carry male_genotype.
    incl_pairs <- included_df %>%
      dplyr::filter(sex_meta == male_hormone) %>%
      dplyr::select(sample_id, cell_type) %>%
      dplyr::distinct() %>%
      dplyr::mutate(pair_key = paste(sample_id, cell_type))

    cell_key    <- paste(meta[[sample_col]], meta[[celltype_col]])
    is_included <- cell_key %in% incl_pairs$pair_key &
                   meta[[seurat_sex_col]] == male_genotype

  } else {
    # E2 / P4: include all cells whose (sample x cell_type) group ran the model.
    incl_pairs <- included_df %>%
      dplyr::select(sample_id, cell_type) %>%
      dplyr::distinct() %>%
      dplyr::mutate(pair_key = paste(sample_id, cell_type))

    cell_key    <- paste(meta[[sample_col]], meta[[celltype_col]])
    is_included <- cell_key %in% incl_pairs$pair_key
  }

  message(sprintf(
    "plot_analysis_umap [%s]: %d included cells (%.1f%%), %d excluded",
    title, sum(is_included),
    100 * mean(is_included),
    sum(!is_included)
  ))

  # ---- 4. Assemble per-cell data frame -------------------------------------
  umap_embed           <- as.data.frame(Embeddings(seurat_obj, reduction = reduction))
  colnames(umap_embed) <- c("UMAP_1", "UMAP_2")

  plot_df <- umap_embed %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::mutate(
      cell_type   = meta[cell_id, celltype_col],
      is_included = is_included[match(cell_id, rownames(meta))]
    )

  plot_df_excl <- dplyr::filter(plot_df, !is_included)
  plot_df_incl <- dplyr::filter(plot_df,  is_included)

  # ---- 5. Cell type label positions ----------------------------------------
  # Centroids computed over ALL cells so positions match the total UMAP exactly.
  centroid_df <- plot_df %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      UMAP_1      = median(UMAP_1),
      UMAP_2      = median(UMAP_2),
      any_included = any(is_included),
      .groups     = "drop"
    )

  centroid_incl <- dplyr::filter(centroid_df,  any_included)
  centroid_excl <- dplyr::filter(centroid_df, !any_included)

  # ---- 6. Build ggplot2 UMAP -----------------------------------------------
  # Layer order: excluded (grey, bottom) -> included (colored, top) -> labels.
  # scale_color_manual is shared by geom_point and the included geom_text layer;
  # the excluded geom_text layer uses a fixed color outside aes() so it does
  # not conflict with the scale.
  # The `breaks` argument limits legend entries to included cell types only.
  p <- ggplot() +
    # Excluded cells: grey, drawn first so colored cells render on top
    geom_point(data    = plot_df_excl,
               mapping = aes(x = UMAP_1, y = UMAP_2),
               color   = excluded_color,
               size    = pt_size,
               alpha   = 0.5) +
    # Included cells: colored by cell type
    geom_point(data    = plot_df_incl,
               mapping = aes(x = UMAP_1, y = UMAP_2, color = cell_type),
               size    = pt_size,
               alpha   = 0.8) +
    # Color scale: values cover all cell types; legend shows included types only
    scale_color_manual(
      values = ct_color_map,
      breaks = centroid_incl$cell_type,
      name   = "Cell Type"
    ) +
    # Labels for excluded cell types at their cluster centroids (grey text)
    geom_text(data     = centroid_excl,
              mapping  = aes(x = UMAP_1, y = UMAP_2, label = cell_type),
              color    = label_excl_color,
              size     = label_size,
              fontface = "bold") +
    # Labels for included cell types; color inherited from scale_color_manual
    geom_text(data        = centroid_incl,
              mapping     = aes(x = UMAP_1, y = UMAP_2,
                                label = cell_type, color = cell_type),
              size        = label_size,
              fontface    = "bold",
              show.legend = FALSE) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_cowplot(12) +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold", size = 13)
    )

  if (!is.null(output_dir) && !is.null(filename_base))
    save_dual_format(p, output_dir, filename_base, width, height)

  invisible(p)
}
