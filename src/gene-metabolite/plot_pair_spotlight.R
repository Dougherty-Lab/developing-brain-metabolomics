# plot_pair_spotlight.R
# ---------------------------------------------------------------------------
# Generate metabolite-gene scatter plots for spotlight pairs from the
# global-FDR hit list.
#
# Inputs:  Association CSVs, cached pseudobulk and metabolite data
# Outputs: Scatter PNGs/SVGs in results/gene-metabolite/figures-log2_na/spotlight/
#
# Upstream:  metabolite_gene_associations.Rmd
# Downstream: None
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(patchwork)
  library(cowplot)
  library(svglite)
  library(readxl)
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
source(file.path(root, "src/gene-metabolite/pseudobulk_functions.R"))   # save_dual_format(), merge_hormone_metadata()
source(file.path(root, "src/gene-metabolite/metabolite_functions.R"))   # read_metabolite_matrix(), transform_metabolite_matrix(), etc.

# ── CONFIG ────────────────────────────────────────────────────────────────────
METHOD       <- "log2_na"
SAVE_OUTPUT  <- TRUE

cache_dir    <- file.path(root, "data/cache")
metab_path   <- file.path(root, "results/untargeted/peak_area_clean.csv")
hormone_xlsx <- file.path(root, "doc/targeted/targeted_hormones.xlsx")
hits_csv     <- file.path(root, "results/gene-metabolite",
                          paste0("csv-", METHOD),
                          "metabolite_gene_hits.csv")
scatter_dir  <- file.path(root, "results/gene-metabolite",
                          paste0("figures-", METHOD), "spotlight")
# ──────────────────────────────────────────────────────────────────────────────

# ── LOAD HITS CSV (model results) ─────────────────────────────────────────────
hits <- read_csv(hits_csv, show_col_types = FALSE)
cat(sprintf("Loaded %s hits from %s\n\n", format(nrow(hits), big.mark = ","),
            basename(hits_csv)))

# ── LOAD DATA OBJECTS (needed for scatter panels) ─────────────────────────────
pb_base     <- readRDS(file.path(cache_dir, "pb_base.rds"))
batch_props <- readRDS(file.path(cache_dir, "batch_props_sample.rds"))

targeted_comparison_raw <- read_xlsx(hormone_xlsx)
targeted_comparison_raw$E2 <- targeted_comparison_raw$`17B-E2 (ng/mL)`
targeted_comparison_raw$TT <- targeted_comparison_raw$`TT (ng/mL)`
targeted_comparison_raw$P4 <- targeted_comparison_raw$`P4 (ng/mL)`

pb        <- build_metab_covariates(pb_base, targeted_comparison_raw, batch_props, min_cells = 30)
metab_mat <- read_metabolite_matrix(metab_path)
metab_tx  <- transform_metabolite_matrix(metab_mat, method = METHOD)  # log2_na predictor

# Raw log2 matrix (all observed samples, no imputation — panel 1 x-axis)
metab_log <- log2(metab_mat + 1)

# Log2 with IQR outliers flagged NA — panels 2-3 x-axis.
# Same fence as log2_na but skips z-scoring to keep interpretable log2 units.
.log2_trim_only <- function(x) {
  logged  <- log2(x + 1)
  orig_na <- is.na(logged)
  obs     <- logged[!orig_na]
  if (length(obs) < 2) return(rep(NA_real_, length(x)))
  qs  <- quantile(obs, c(0.25, 0.75), names = FALSE)
  iqr <- qs[2] - qs[1]
  v   <- logged
  v[!orig_na & (logged < qs[1] - 0.75 * iqr | logged > qs[2] + 0.75 * iqr)] <- NA
  v
}
metab_log_trim <- t(apply(metab_mat, 1, .log2_trim_only))
dimnames(metab_log_trim) <- dimnames(metab_mat)

# Global z-score x-axis: symmetric range covering 99.8% of z-scores across
# all metabolites. Fixed so z-score panels are directly comparable across plots.
z_vals   <- as.vector(metab_tx)
z_abs    <- quantile(abs(z_vals), 0.999, na.rm = TRUE)
z_x_lim  <- c(-z_abs, z_abs)

# Per-pseudobulk helpers
rn        <- rownames(pb$counts)
row_ct    <- sub("^[^.]+\\.(.*)$", "\\1", rn)          # cell type from rowname
metab_col <- metab_sample_to_pb(colnames(metab_tx),     # SampleNN → pb row mapping
                                 sub("^([^.]+)\\..*$", "\\1", rn))
cc_col    <- intersect(c("cell_count", "n_cells", "ncells"), colnames(pb$metadata))[1]

if (SAVE_OUTPUT) dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)

# ── print_model_info() ────────────────────────────────────────────────────────
# Formats the limma-voom result row from the hits CSV to the console.
print_model_info <- function(row) {
  sep <- paste(rep("─", 62), collapse = "")
  cat(sep, "\n")
  cat(sprintf("Gene:       %s\n", row$gene))
  cat(sprintf("Metabolite: %s  (%s)\n", row$Name, row$Compound.ID))
  cat(sprintf("Cell type:  %s\n", row$cell_type))
  cat(sprintf("Class:      %s › %s › %s\n",
              row$Super.Class %||% "NA", row$Class %||% "NA", row$Sub.Class %||% "NA"))

  cat(sprintf("\nModel  (limma-voom, ~ metab_%s + GW + Sex + batch1_frac):\n", METHOD))
  cat(sprintf("  logFC       = %+.4f   (Δlog2-CPM per 1 SD of metabolite)\n", row$logFC))
  cat(sprintf("  AveExpr     = %+.4f   (mean log2-CPM across samples)\n",     row$AveExpr))
  cat(sprintf("  t           = %+.4f\n",    row$t))
  cat(sprintf("  B           = %+.4f   (log-odds gene is DE)\n",              row$B))
  cat(sprintf("  P.Value     = %.3g\n",     row$P.Value))
  cat(sprintf("  adj.P.Val   = %.3g   (BH within metabolite × cell type)\n", row$adj.P.Val))
  cat(sprintf("  FDR_global  = %.3g   (BH across all tests)\n",              row$FDR_global))

  if ("max_dfb" %in% colnames(row) && !is.na(row$max_dfb)) {
    cat(sprintf("\nInfluence  (DFBETAS gate: flag if |DFBETAS| > 2/√n):\n"))
    cat(sprintf("  max_dfb     = %.4f   cutoff = %.4f   stable = %s\n",
                row$max_dfb, row$dfb_cut,
                ifelse(isTRUE(row$stable), "YES ✓", "NO ✗")))
  }

  if (!is.na(row$Annotation_Confidence))
    cat(sprintf("\nAnnotation confidence: %s\n", row$Annotation_Confidence))
  if (!is.na(row$KEGG) && nchar(as.character(row$KEGG)) > 0)
    cat(sprintf("KEGG:       %s\n", row$KEGG))
  if (!is.na(row$Pathways) && nchar(as.character(row$Pathways)) > 0)
    cat(sprintf("Pathways:   %s\n", row$Pathways))

  cat(sep, "\n\n")
}

# Null-coalescing helper (base R lacks %||%)
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

# ── make_pair_plot() ──────────────────────────────────────────────────────────
# Five-panel scatter mirroring metabolite_gene_stability_dfbetas.qmd exactly.
# Row 1: raw log2 (all samples) | outlier-trimmed log2 | covariate-adjusted log2
# Row 2:           [spacer]     | outlier-trimmed z     | covariate-adjusted z
# One pooled regression line per panel (Sex is a model covariate, not a group).
make_pair_plot <- function(cid, gene, ct, meta_name) {
  idx      <- which(row_ct == ct)
  pos      <- match(metab_col[idx], colnames(metab_tx))
  present  <- !is.na(pos) & gene %in% colnames(pb$counts)
  rows_use <- idx[present]
  pos_p    <- pos[present]

  if (length(rows_use) < 3) {
    warning(sprintf("Fewer than 3 samples with metabolite data for %s ~ %s (%s) — skipping.",
                    gene, meta_name, ct))
    return(NULL)
  }

  expr <- as.numeric(
    log1p(edgeR::cpm(t(pb$counts[rows_use, , drop = FALSE]), log = FALSE)[gene, ])
  )
  md <- pb$metadata[rn[rows_use], ]

  df <- tibble(
    metab_raw     = as.numeric(metab_log[cid, pos_p]),       # raw log2 (NAs for missing)
    metab_logtrim = as.numeric(metab_log_trim[cid, pos_p]),  # log2, outliers = NA
    metab_trim    = as.numeric(metab_tx[cid, pos_p]),        # z-score (model predictor)
    expr          = expr,
    Sex           = md$Sex,
    GW            = md$GW,
    batch1_frac   = md$batch1_frac
  )

  covs <- c("GW", "Sex", "batch1_frac")

  # ---- log2 panels (free x-axis per metabolite) ----
  x_range <- range(df$metab_raw, na.rm = TRUE)
  x_pad   <- diff(x_range) * 0.05
  x_lim   <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  x_lab   <- paste0(meta_name, " [log2 peak area]")

  bl_log <- list(
    geom_point(aes(shape = Sex), size = 2.5, alpha = 0.85),
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8),
    scale_color_viridis_c(name = "GW"),
    scale_shape_manual(values = c("M" = 16, "F" = 17)),
    coord_cartesian(xlim = x_lim),
    theme_minimal(base_size = 13)
  )

  # Panel 1: raw log2 — all observed samples, before trimming
  df_raw <- df |> dplyr::filter(!is.na(metab_raw))
  p_raw <- ggplot(df_raw, aes(metab_raw, expr, color = GW)) + bl_log +
    labs(title = "Raw log2 (all samples)", x = x_lab, y = "log1p(CPM)")

  # Panel 2: outlier-trimmed log2
  dft <- df |> dplyr::filter(!is.na(metab_logtrim))
  keepc <- covs[vapply(covs, function(cv) length(unique(dft[[cv]])) >= 2, logical(1))]
  p_mid <- ggplot(dft, aes(metab_logtrim, expr, color = GW)) + bl_log +
    labs(title = "Outlier-trimmed (log2)", x = x_lab, y = "log1p(CPM)")

  # Panel 3: covariate-adjusted log2 (partial regression)
  fit_log  <- lm(reformulate(c("metab_logtrim", keepc), response = "expr"), data = dft)
  dft$part <- resid(fit_log) + coef(fit_log)[["metab_logtrim"]] * dft$metab_logtrim
  p_adj <- ggplot(dft, aes(metab_logtrim, part, color = GW)) + bl_log +
    labs(title = "Covariate-adjusted (log2)", x = x_lab, y = "log1p(CPM), adjusted")

  # ---- z-score panels (fixed x-axis across all plots) ----
  z_lab <- paste0(meta_name, " [z-score]")

  bl_z <- list(
    geom_point(aes(shape = Sex), size = 2.5, alpha = 0.85),
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8),
    scale_color_viridis_c(name = "GW"),
    scale_shape_manual(values = c("M" = 16, "F" = 17)),
    coord_cartesian(xlim = z_x_lim),
    theme_minimal(base_size = 13)
  )

  # Panel 4: outlier-trimmed z-score
  dft_z   <- df |> dplyr::filter(!is.na(metab_trim))
  keepc_z <- covs[vapply(covs, function(cv) length(unique(dft_z[[cv]])) >= 2, logical(1))]
  p_zmid  <- ggplot(dft_z, aes(metab_trim, expr, color = GW)) + bl_z +
    labs(title = "Outlier-trimmed (z-score)", x = z_lab, y = "log1p(CPM)")

  # Panel 5: covariate-adjusted z-score
  fit_z        <- lm(reformulate(c("metab_trim", keepc_z), response = "expr"), data = dft_z)
  dft_z$part_z <- resid(fit_z) + coef(fit_z)[["metab_trim"]] * dft_z$metab_trim
  p_zadj <- ggplot(dft_z, aes(metab_trim, part_z, color = GW)) + bl_z +
    labs(title = "Covariate-adjusted (z-score)", x = z_lab, y = "log1p(CPM), adjusted")

  (p_raw | p_mid | p_adj) / (plot_spacer() | p_zmid | p_zadj) +
    plot_annotation(
      title = paste0(gene, "  ~  ", meta_name, "   (", ct, ")"),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    ) +
    plot_layout(guides = "collect")
}

# ── MAIN LOOP ─────────────────────────────────────────────────────────────────
for (p in PAIRS) {
  cid <- p[["Compound.ID"]]
  gn  <- p[["gene"]]
  ct  <- p[["cell_type"]]

  # Look up the hits row for model stats
  row <- hits |> dplyr::filter(Compound.ID == cid, gene == gn, cell_type == ct)

  if (nrow(row) == 0) {
    warning(sprintf("\n[SKIP] No row found in hits CSV for: %s / %s / %s\n", cid, gn, ct))
    next
  }
  if (nrow(row) > 1) {
    warning(sprintf("[WARN] Multiple rows for %s / %s / %s — using first.\n", cid, gn, ct))
    row <- row[1L, ]
  }

  # 1. Print model summary to console
  print_model_info(row)

  # 2. Build and display scatter
  pl <- make_pair_plot(cid, gn, ct, row$Name)
  if (is.null(pl)) next
  print(pl)

  # 3. Optionally save PNG + SVG
  if (SAVE_OUTPUT) {
    base_name <- gsub("[^A-Za-z0-9]+", "-", paste(cid, gn, ct, sep = "_"))
    save_dual_format(pl, scatter_dir, base_name, width = 16, height = 10)
    cat(sprintf("Saved → %s\n\n", file.path(scatter_dir, base_name)))
  }
}

# ---- AI assistance disclosure ------------------------------------------------
# Code in this script was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.

# ---- session info ------------------------------------------------------------
cat("\n\n---- Session Info ----\n")
print(sessionInfo())
