# qq_stable_associations.R
# ---------------------------------------------------------------------------
# QQ plots of p-value distributions from the DFBETAS-stable association
# subset, checking for residual inflation.
#
# Inputs:  Stable-hit CSVs from metabolite_gene_imputation_sensitivity.Rmd
# Outputs: QQ plot figures in results/gene-metabolite/stability-dfbetas-log2_na/
#
# Upstream:  metabolite_gene_imputation_sensitivity.Rmd
# Downstream: None
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(cowplot)
  library(svglite)
  library(patchwork)
})

# Anchor paths to project root (matches gene_pathway_analysis.qmd convention)
find_project_root <- function(marker = ".git") {
  d <- normalizePath(getwd())
  repeat {
    if (dir.exists(file.path(d, marker))) return(d)
    parent <- dirname(d)
    if (parent == d) stop("Could not locate project root (no ", marker, " found above ", getwd(), ")")
    d <- parent
  }
}
results_dir <- file.path(find_project_root(), "results")

source(file.path(find_project_root(), "src", "gene-metabolite", "pseudobulk_functions.R"))

set.seed(123)

# ── Paths ──────────────────────────────────────────────────────────────────────
METHOD      <- "log2_na"
suffix      <- paste0("-", METHOD)
parquet_dir <- file.path(results_dir, sprintf("gene-metabolite/parquet-%s", METHOD))
csv_dir     <- file.path(results_dir, "gene-metabolite", paste0("csv", suffix))
hits_csv    <- file.path(csv_dir, "metabolite_gene_hits.csv")
qq_dir      <- file.path(csv_dir, "qq-plots")
dir.create(qq_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ─────────────────────────────────────────────────────────────────
MIN_TESTS   <- 100    # skip cell types with fewer tests (QQ shape unreliable)
TAIL_THRESH <- 2.0    # -log10(p) threshold: points above are "tail" and kept in full
THIN_FRAC   <- 0.01   # fraction of null body kept for plotting (reduces overplotting)

# ── 1. Load global-FDR hits ─────────────────────────────────────────────────
fdr_hits <- read_csv(hits_csv, show_col_types = FALSE)
cat(sprintf("Global-FDR hits: %d rows | %d cell types\n",
            nrow(fdr_hits), n_distinct(fdr_hits$cell_type)))

# ── 2. Load full p-value universe from parquet via Arrow ───────────────────────
# One parquet per metabolite; Arrow opens as a lazy dataset.
# Collecting only P.Value + cell_type avoids loading the full result table.
cat("Opening parquet dataset …\n")
ds <- open_dataset(parquet_dir)

cat("Collecting P.Value + cell_type across all metabolites …\n")
all_pvals_df <- ds |>
  dplyr::select(P.Value, cell_type) |>
  dplyr::filter(!is.na(P.Value)) |>
  collect()

cat(sprintf("Loaded %s tests across %d cell types\n",
            format(nrow(all_pvals_df), big.mark = ","),
            n_distinct(all_pvals_df$cell_type)))

# ── 3. Helper functions ────────────────────────────────────────────────────────

#' Genomic inflation factor λ.
#' λ = median(χ²_obs) / 0.4549, where 0.4549 = median of χ²(1) under H0.
#' Values >> 1 indicate model inflation; < 1 indicate overdispersion correction.
compute_lambda <- function(pvals) {
  chisq_obs <- qchisq(1 - pvals[pvals > 0 & pvals < 1], df = 1)
  median(chisq_obs, na.rm = TRUE) / qchisq(0.5, df = 1)
}

#' Build one per-cell-type QQ plot.
#'
#' @param ct       Cell type label.
#' @param pvals    All p-values for this cell type (no NAs).
#' @param hit_p    P-values of global-FDR hits in this cell type (length ≥ 0).
#' @param lambda   Pre-computed λ.
#' @return ggplot object.
make_qq_plot <- function(ct, pvals, hit_p, lambda) {
  n        <- length(pvals)
  p_sorted <- sort(pvals)               # ascending: rank 1 = most significant
  obs_all  <- -log10(p_sorted)
  exp_all  <- -log10(seq_len(n) / (n + 1L))   # expected U(0,1) quantiles

  # ── Thinning ──────────────────────────────────────────────────────────────
  # Keep full tail (obs > TAIL_THRESH) + a random sample of the null body.
  # Global-FDR hits are added back explicitly so none are dropped.
  idx_tail    <- which(obs_all > TAIL_THRESH)
  idx_body    <- which(obs_all <= TAIL_THRESH)
  n_body_keep <- max(200L, ceiling(length(idx_body) * THIN_FRAC))
  idx_body_smp <- sort(sample(idx_body, size = min(n_body_keep, length(idx_body))))

  # Locate global-FDR hits in sorted p-value vector and force-include them
  hit_ranks <- integer(0L)
  if (length(hit_p) > 0L) {
    hit_ranks <- findInterval(sort(hit_p), p_sorted)
    hit_ranks <- pmax(1L, pmin(n, hit_ranks))
  }

  keep_idx <- sort(unique(c(idx_tail, idx_body_smp, hit_ranks)))

  # ── 95% CI ribbon ─────────────────────────────────────────────────────────
  # j-th order statistic of U(0,1) ~ Beta(j, n-j+1).
  # Ribbon sits around the diagonal (y = x), showing expected range under H0.
  # Computed only at kept indices (avoids n calls to qbeta).
  ci_df <- tibble(
    x    = exp_all[keep_idx],
    ymin = -log10(qbeta(0.975, keep_idx, n - keep_idx + 1L)),
    ymax = -log10(qbeta(0.025, keep_idx, n - keep_idx + 1L))
  )

  # ── Plot data frame ────────────────────────────────────────────────────────
  plot_df <- tibble(
    expected = exp_all[keep_idx],
    observed = obs_all[keep_idx],
    is_hit   = keep_idx %in% hit_ranks
  )

  max_val <- max(c(plot_df$observed, plot_df$expected), na.rm = TRUE) * 1.05

  # ── Build ggplot ───────────────────────────────────────────────────────────
  ggplot(plot_df, aes(expected, observed)) +
    # CI ribbon around the null diagonal
    geom_ribbon(
      data        = ci_df,
      aes(x = x, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill        = "steelblue",
      alpha       = 0.20
    ) +
    # Null diagonal y = x
    geom_abline(slope = 1, intercept = 0,
                colour = "grey50", linewidth = 0.45, linetype = "dashed") +
    # Null-body points (thinned)
    geom_point(
      data    = dplyr::filter(plot_df, !is_hit),
      size    = 0.7, alpha = 0.35, colour = "grey55"
    ) +
    # Global-FDR hits highlighted in the tail
    geom_point(
      data    = dplyr::filter(plot_df, is_hit),
      size    = 1.8, alpha = 0.90, colour = "firebrick"
    ) +
    # λ annotation (top-left)
    annotate(
      "text",
      x      = 0.03 * max_val,
      y      = 0.97 * max_val,
      label  = sprintf("\u03bb = %.3f\nn = %s\nFDR hits = %d",
                       lambda,
                       format(n, big.mark = ","),
                       length(hit_p)),
      hjust  = 0, vjust = 1,
      size   = 3.2, colour = "grey20"
    ) +
    coord_fixed(xlim = c(0, max_val), ylim = c(0, max_val)) +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = 0.01)) +
    labs(
      title = ct,
      x     = expression(Expected ~ -log[10](italic(p))),
      y     = expression(Observed ~ -log[10](italic(p)))
    ) +
    theme_cowplot(11) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

# ── 4. Loop per cell type ──────────────────────────────────────────────────────
cell_types     <- sort(unique(all_pvals_df$cell_type))
lambda_records <- list()
qq_plots       <- list()

cat("\nBuilding QQ plots:\n")
for (ct in cell_types) {

  pvals_ct <- all_pvals_df$P.Value[all_pvals_df$cell_type == ct]
  pvals_ct <- pvals_ct[is.finite(pvals_ct) & !is.na(pvals_ct)]

  if (length(pvals_ct) < MIN_TESTS) {
    cat(sprintf("  [skip] '%s'  n = %d < MIN_TESTS (%d)\n",
                ct, length(pvals_ct), MIN_TESTS))
    next
  }

  hit_p_ct  <- fdr_hits$P.Value[fdr_hits$cell_type == ct]
  lambda_ct <- compute_lambda(pvals_ct)

  cat(sprintf("  %-30s  n = %s  |  \u03bb = %.3f  |  FDR hits = %d\n",
              ct, format(length(pvals_ct), big.mark = ","),
              lambda_ct, length(hit_p_ct)))

  lambda_records[[ct]] <- tibble(
    cell_type = ct,
    n_tests   = length(pvals_ct),
    lambda    = lambda_ct,
    n_hits    = length(hit_p_ct)
  )

  p_ct <- make_qq_plot(ct, pvals_ct, hit_p_ct, lambda_ct)
  qq_plots[[ct]] <- p_ct

  safe_ct <- gsub("[^A-Za-z0-9]+", "-", ct)
  save_dual_format(p_ct, qq_dir, paste0("qq_", safe_ct), width = 5, height = 5)
}

# ── 5. λ summary CSV ───────────────────────────────────────────────────────────
lambda_df <- bind_rows(lambda_records) |>
  dplyr::arrange(dplyr::desc(lambda))

write_csv(lambda_df, file.path(qq_dir, "qq_lambda_summary.csv"))
cat("\n\u03bb summary (descending):\n")
print(lambda_df, n = Inf)

# Flag any cell types where λ warrants scrutiny (> 1.5 is a loose alarm)
flagged <- dplyr::filter(lambda_df, lambda > 1.5)
if (nrow(flagged) > 0L) {
  cat(sprintf("\n! %d cell type(s) with \u03bb > 1.5 \u2014 inspect individually:\n", nrow(flagged)))
  print(flagged$cell_type)
}

# ── 6. Combined patchwork grid ──────────────────────────────────────────────────
n_plots <- length(qq_plots)
if (n_plots > 0L) {
  n_col  <- min(3L, n_plots)
  n_row  <- ceiling(n_plots / n_col)

  combo <- wrap_plots(qq_plots, ncol = n_col) +
    plot_annotation(
      title    = "Gene\u2013Metabolite Associations: QQ plots by cell type",
      subtitle = sprintf(
        "Red = global-FDR hits (FDR < 0.10, |logFC| \u2265 0.25) | %s transform | grey ribbon = 95%% CI under H\u2080",
        METHOD
      ),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "grey40")
      )
    )

  save_dual_format(combo, qq_dir, "qq_combined",
                   width  = n_col * 5,
                   height = n_row * 5 + 0.8)
  cat(sprintf("\nCombined grid: %d col \u00d7 %d row saved to %s/qq_combined.png\n",
              n_col, n_row, qq_dir))
}

cat(sprintf("\nDone. All outputs in:\n  %s\n", normalizePath(qq_dir, mustWork = FALSE)))

# ---- AI assistance disclosure ------------------------------------------------
# Code in this script was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.

# ---- session info ------------------------------------------------------------
cat("\n\n---- Session Info ----\n")
print(sessionInfo())
