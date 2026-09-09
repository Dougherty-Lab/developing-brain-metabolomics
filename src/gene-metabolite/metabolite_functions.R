# metabolite_functions.R
# ---------------------------------------------------------------------------
# Shared helpers for metabolite-gene ASSOCIATION analysis: snRNA-seq pseudobulk
# gene expression vs. untargeted metabolite abundances.
# Sourced alongside pseudobulk_functions.R (reuses .safe_vars, merge_hormone_metadata).
#
# Model orientation
#   Gene counts are the voom RESPONSE; the metabolite is a continuous PREDICTOR.
#   The reported logFC is therefore d(gene log2-CPM) per unit of the metabolite
#   predictor. Because voom precision weights are defined for count responses
#   (Law et al. 2014, Genome Biology 15:R29), the metabolite cannot be the
#   response here, and regression association is symmetric -- results are
#   ASSOCIATIONS, not directional predictions.
#
# Metabolite transform (matches the Batch2 metabolomics pipeline, extended with
# z-scoring), applied per metabolite:
#   log2(x + 1)  ->  impute NaN with min(non-NA)/5 on the log2 scale
#                ->  z-score (mean 0, sd 1) across all samples
#   Imputation order/maths mirror functions.R::min_value_impute exactly.
#   Z-scoring standardises the predictor so |logFC| is comparable across the
#   ~600 metabolites (Gelman 2008, Stat Med 27:2865). It rescales only the
#   predictor ("per 1 SD" instead of "per 1 log2-unit"); the gene-side
#   numerator is unchanged.
#
# Design:  ~ metabolite + GW + Sex + batch1_frac  (split by cell type only)
#   batch1_frac is the continuous per-sample sequencing-batch fraction (the new
#   annotation's settled batch covariate). Any covariate that is single-level /
#   zero-variance within a cell-type subset is dropped for that subset; the
#   metabolite term is always kept.
# ---------------------------------------------------------------------------


# ---- 1. Imputation (identical maths to functions.R::min_value_impute) --------
#' Replace NaN/NA in a numeric vector with one fifth of its smallest observed
#' value. Applied per metabolite (row-wise below), on the log2 scale.
min_value_impute <- function(column) {
  min_value <- min(column, na.rm = TRUE)
  column[is.na(column)] <- min_value / 5
  column
}


# ---- 2. Load + transform the metabolite matrix -------------------------------
#' Read a cleaned peak-area CSV into a (features x samples) numeric matrix.
#'
#' @param path CSV with columns Compound.ID, Name, Formula, and one column per
#'             sample named "Sample01", "Sample02", ...  NaN = metabolite absent.
#' @return numeric matrix, rownames = Compound.ID, columns = SampleNN.
#'         attr(, "name_lookup") maps Compound.ID -> Name for labelling.
read_metabolite_matrix <- function(path) {
  raw <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot("Compound.ID" %in% colnames(raw), "Name" %in% colnames(raw))

  sample_cols <- grep("^Sample\\d+$", colnames(raw), value = TRUE)
  if (length(sample_cols) == 0)
    stop("No ^Sample\\d+$ columns found in ", path)

  mat <- as.matrix(raw[, sample_cols])
  storage.mode(mat) <- "numeric"        # coerces character/NA cleanly
  rownames(mat)     <- raw$Compound.ID

  attr(mat, "name_lookup") <- stats::setNames(raw$Name, raw$Compound.ID)
  mat
}

#' Rank-based inverse-normal transform (INT), Blom offset, average-rank ties.
#' Maps values to normal quantiles by rank, so no single sample can sit at
#' extreme leverage (bounded to ~+/-2 at n=20). Monotonic -> association
#' directions are unchanged vs. the raw predictor (McCaw et al. 2020,
#' Biometrics 76:1262). Applied to the metabolite because it is a PREDICTOR
#' here; unrelated to the metabolite's role as a response in the sex/GW models.
rank_int <- function(x) {
  n <- sum(!is.na(x))
  r <- rank(x, ties.method = "average", na.last = "keep")
  stats::qnorm((r - 0.375) / (n + 0.25))
}

#' Per-metabolite transform, four selectable methods (row-wise, features x samples):
#'   "int"          min/5 impute -> rank-inverse-normal      (leverage-robust; default)
#'   "zscore"       log2(x+1) -> min/5 impute -> z-score      (original)
#'   "zscore_trim"  log2(x+1) -> IQR-trim observed outliers (excluded, set NA)
#'                  -> min/5 impute the originally-missing -> z-score
#'   "log2_na"      log2(x+1) -> IQR-trim observed outliers (excluded, set NA)
#'                  -> z-score; NO imputation -- originally-missing AND outliers
#'                  stay NA and are dropped by analyze_metabolite()'s
#'                  !is.na(md$metab) filter. This is the primary imputation-free
#'                  run; MIN_SAMPLES=10 is applied post-NA-drop as usual.
#'
#' For "zscore_trim" and "log2_na", outliers are EXCLUDED (kept NA), not
#' winsorized. Trimming runs on observed values BEFORE any imputation so fences
#' reflect true data. IQR fence is Q1 - 0.75*IQR / Q3 + 0.75*IQR, matching
#' the hormone-measurement convention (stricter than Tukey's 1.5x;
#' Tukey 1977, Exploratory Data Analysis).
transform_metabolite_matrix <- function(mat, method = c("int", "zscore", "zscore_trim", "log2_na")) {
  method <- match.arg(method)
  fn <- switch(method,
    int         = .tf_int,
    zscore      = .tf_zscore,
    zscore_trim = .tf_zscore_trim,
    log2_na     = .tf_log2_na
  )
  out <- t(apply(mat, 1, fn))
  dimnames(out)            <- dimnames(mat)
  attr(out, "name_lookup") <- attr(mat, "name_lookup")
  out
}

.zscore_vec <- function(v) {
  s <- stats::sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(v)))
  (v - mean(v, na.rm = TRUE)) / s          # NA positions stay NA -> excluded
}

.tf_int <- function(x) rank_int(min_value_impute(x))

.tf_zscore <- function(x) .zscore_vec(min_value_impute(log2(x + 1)))

.tf_zscore_trim <- function(x) {
  logged  <- log2(x + 1)
  orig_na <- is.na(logged)                 # originally-missing samples
  obs     <- logged[!orig_na]
  if (length(obs) < 2) return(rep(NA_real_, length(x)))

  qs  <- stats::quantile(obs, c(0.25, 0.75), names = FALSE)
  iqr <- qs[2] - qs[1]
  lo  <- qs[1] - 0.75 * iqr
  hi  <- qs[2] + 0.75 * iqr

  v <- logged
  v[!orig_na & (logged < lo | logged > hi)] <- NA   # trim outliers -> excluded
  v[orig_na] <- min(v, na.rm = TRUE) / 5            # impute ONLY original missing
  .zscore_vec(v)
}

# Imputation-free variant: same IQR outlier fence as zscore_trim, but
# originally-missing values are NOT imputed. Both outliers and missing stay NA
# and are dropped per-metabolite by analyze_metabolite()'s !is.na(md$metab)
# filter, so the model only sees genuinely observed, non-extreme values.
.tf_log2_na <- function(x) {
  logged  <- log2(x + 1)
  orig_na <- is.na(logged)                 # originally-missing samples
  obs     <- logged[!orig_na]
  if (length(obs) < 2) return(rep(NA_real_, length(x)))

  qs  <- stats::quantile(obs, c(0.25, 0.75), names = FALSE)
  iqr <- qs[2] - qs[1]
  lo  <- qs[1] - 0.75 * iqr
  hi  <- qs[2] + 0.75 * iqr

  v <- logged
  v[!orig_na & (logged < lo | logged > hi)] <- NA   # trim outliers -> excluded
  # NO imputation: orig_na positions remain NA
  .zscore_vec(v)
}


# ---- 3. Covariate metadata (Sex / GW / batch1_frac on the pb sample set) -----
#' Attach Sex + GW (from the targeted hormone sheet) and batch1_frac (cached) to
#' the pseudobulk sample set, reusing the existing hormone merge path so
#' covariates are identical to the hormone-gene models.
#'
#' @param pb_base      build_pseudobulk() output (cached pb_base.rds).
#' @param hormone_df   targeted_comparison_raw; needs Sample, Sex, GW
#'                     (and TT/E2/P4, which merge_hormone_metadata selects).
#' @param batch_props  batch_props_sample.rds; columns Sample, batch1_frac,
#'                     where Sample matches pb sample_id (e.g. "SSD07_tissue7").
#' @param min_cells    minimum cells per (sample x cell type). Default 30.
#' @return pb list (counts / logcounts / metadata) with Sex, GW, batch1_frac on
#'         $metadata, aligned to the surviving pseudobulk samples.
build_metab_covariates <- function(pb_base, hormone_df, batch_props,
                                   min_cells = 30) {

  pb <- merge_hormone_metadata(
    pb_base, hormone_df,
    pb_sample_pattern = "tissue\\d+",
    min_cells         = min_cells
  )

  pb$metadata <- pb$metadata %>%
    tibble::rownames_to_column("pseudobulk_id") %>%
    dplyr::left_join(batch_props, by = c("sample_id" = "Sample")) %>%
    tibble::column_to_rownames("pseudobulk_id")

  n_missing <- sum(is.na(pb$metadata$batch1_frac))
  if (n_missing > 0)
    warning(sprintf(
      "build_metab_covariates: %d pseudobulk sample(s) missing batch1_frac.",
      n_missing
    ))

  pb
}


# ---- 4. Map SampleNN metabolite columns to pb sample_ids ---------------------
#' Build the per-pb-sample metabolite-column lookup. Metabolite columns are
#' "Sample07"; pb sample_ids are "SSD07_tissue7". Both key on the integer:
#' Sample07 <-> tissue7 (zero-padding stripped).
#'
#' @param metab_cols     colnames of the metabolite matrix ("Sample01", ...).
#' @param pb_sample_ids  pb$metadata$sample_id, in row order.
#' @return character vector (length == length(pb_sample_ids)) giving the
#'         metabolite column name for each pb sample, NA where absent.
metab_sample_to_pb <- function(metab_cols, pb_sample_ids) {
  tissue_of_metab        <- paste0("tissue", as.integer(sub("^Sample", "", metab_cols)))
  names(metab_cols)      <- tissue_of_metab
  pb_tissue              <- stringr::str_extract(pb_sample_ids, "tissue\\d+")
  unname(metab_cols[pb_tissue])
}


# ---- 5. Model wrapper: ~ metab + GW + Sex + batch1_frac ----------------------
#' Fit one voom model for a single metabolite within one already-subset group.
#' Drops any covariate that is single-level / zero-variance in this subset; the
#' metabolite predictor is always retained. Returns NULL if residual df < 1
#' (guards limma's .ebayes() zero-residual-df crash).
run_voom_metab <- function(counts, metadata, predictor = "metab") {
  stopifnot(ncol(counts) == nrow(metadata))

  candidate_cov <- c("GW", "Sex", "batch1_frac")
  keep_cov <- candidate_cov[vapply(candidate_cov, function(cv) {
    if (!cv %in% colnames(metadata)) return(FALSE)
    x <- metadata[[cv]]
    if (all(is.na(x))) return(FALSE)
    length(unique(x[!is.na(x)])) >= 2          # non-constant / >=2 levels
  }, logical(1))]

  formula_vars <- .safe_vars(c(predictor, keep_cov))
  design       <- model.matrix(reformulate(formula_vars), data = metadata)

  # model.matrix silently drops rows with NA in any term; realign counts.
  kept   <- rownames(design)
  counts <- counts[, kept, drop = FALSE]

  if (ncol(counts) < 3 || (ncol(counts) - ncol(design)) < 1) return(NULL)

  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))

  limma::topTable(fit, coef = colnames(design)[2],
                  number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      dropped_covariates = paste(setdiff(candidate_cov, keep_cov), collapse = ";")
    )
}


# ---- 6. One metabolite -> all cell types -------------------------------------
#' Loop a single metabolite across cell types. No hormone-style non-zero
#' predictor filter: all pb samples that carry a (non-NA) metabolite value are
#' used, subject to the gene filter and a required minimum sample count.
#'
#' @param pb_counts  pb$counts (samples x genes), covariate-aligned.
#' @param metadata   pb$metadata (Sex/GW/batch1_frac) with a `metab` column
#'                   holding the standardised value for THIS metabolite.
#' @param compound_id,name  identifiers carried into the output.
#' @param count_threshold,gene_sample_threshold  manual gene filter
#'                   (rowSums(counts >= 2) >= 3). Settled; do not revisit.
#' @param min_samples  required samples per cell-type group (no default; guards
#'                   zero-residual-df fits). Default 10 here for the loop.
#' @return long data.frame: gene, logFC, AveExpr, t, P.Value, adj.P.Val, B,
#'         dropped_covariates, Compound.ID, Name, cell_type, n_samples.
#'         adj.P.Val is the WITHIN (metabolite x cell type) FDR across genes.
analyze_metabolite <- function(pb_counts, metadata,
                               compound_id, name,
                               count_threshold       = 2,
                               gene_sample_threshold = 3,
                               min_samples           = 10,
                               model                 = run_voom_metab) {

  split_ids <- split(rownames(metadata), metadata$cell_type)

  purrr::map_dfr(names(split_ids), function(ct) {
    ids  <- split_ids[[ct]]
    md   <- metadata[ids, , drop = FALSE]
    cnts <- t(pb_counts)[, ids, drop = FALSE]

    # Keep only samples that carry this metabolite (non-NA standardised value).
    ok   <- !is.na(md$metab)
    md   <- md[ok, , drop = FALSE]
    cnts <- cnts[, ok, drop = FALSE]

    if (nrow(md) < min_samples)          return(NULL)   # sample-count guard
    if (length(unique(md$metab)) < 2)    return(NULL)   # constant predictor

    keep <- rowSums(cnts >= count_threshold) >= gene_sample_threshold
    cnts <- cnts[keep, , drop = FALSE]
    if (nrow(cnts) == 0)                 return(NULL)

    res <- model(cnts, md, predictor = "metab")
    if (is.null(res) || nrow(res) == 0)  return(NULL)

    res$Compound.ID <- compound_id
    res$Name        <- name
    res$cell_type   <- ct
    res$n_samples   <- nrow(md)
    res
  })
}

# ---- AI assistance disclosure ------------------------------------------------
# Code in this file was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.
