#!/usr/bin/env Rscript
# mofa_umap_exploration.R
# ---------------------------------------------------------------------------
# Quick exploration script: MOFA factor values + key gene expression on UMAP.
# Generates candidate panels for Figures 5 & 6 upgrades.
#
# Outputs to results/mofa/umap-exploration/ (PNG + SVG via save_dual_format).
#
# Run interactively in RStudio or:
#   Rscript mofa_umap_exploration.R
#
# Code developed by author with assistance from Claude (Opus 4.6).
# Output checked and validated by author.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(MOFA2)
  library(patchwork)
  library(viridis)
  library(svglite)
})

# ---- Paths ------------------------------------------------------------------
# Detect project root via .git marker
find_project_root <- function() {
  path <- getwd()
  while (path != dirname(path)) {
    if (file.exists(file.path(path, ".git"))) return(path)
    path <- dirname(path)
  }
  stop("Could not find project root (.git marker)")
}

PROJECT_ROOT <- find_project_root()
cache_dir    <- file.path(PROJECT_ROOT, "data", "cache")
out_dir      <- file.path(PROJECT_ROOT, "results", "mofa", "umap-exploration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT_ROOT, "src", "gene-metabolite", "pseudobulk_functions.R"))

# ---- TODO: Set your MOFA model path ----------------------------------------
# Update this to wherever mofa_integration.qmd saves the trained model.
# Likely something like:
#   results/mofa/mofa_model.rds
#   results/mofa/mofa_trained_model.hdf5
mofa_path <- file.path(PROJECT_ROOT, "results", "mofa", "mofa_model.hdf5")
# If saved as HDF5:
# mofa_path <- file.path(PROJECT_ROOT, "results", "mofa", "mofa_model.hdf5")

# ---- Load objects -----------------------------------------------------------
cat("Loading Seurat meta+UMAP object...\n")
ssd <- readRDS(file.path(cache_dir, "SSD_meta_umap.rds"))
cat("  Cells:", ncol(ssd), "\n")

cat("Loading MOFA model...\n")
if (grepl("\\.hdf5$", mofa_path)) {
  mofa <- load_model(mofa_path)
} else {
  mofa <- readRDS(mofa_path)
}
cat("  Factors:", MOFA2::get_dimensions(mofa)$K, "\n")
cat("  Samples:", MOFA2::get_dimensions(mofa)$N, "\n")

# ---- Extract MOFA factor values per sample ----------------------------------
# Factor matrix: samples x factors
Z <- MOFA2::get_factors(mofa)[[1]]  # first group (should be only group)
factor_df <- as.data.frame(Z) %>%
  tibble::rownames_to_column("sample_id")

cat("\nMOFA factor matrix dimensions:", nrow(factor_df), "x", ncol(factor_df) - 1, "\n")

# ---- Extract UMAP coordinates -----------------------------------------------
umap_embed <- as.data.frame(Embeddings(ssd, reduction = "umap.rna"))
colnames(umap_embed) <- c("UMAP_1", "UMAP_2")

cell_df <- umap_embed %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    Sample   = ssd@meta.data[cell_id, "Sample"],
    celltype = ssd@meta.data[cell_id, "celltype"],
    Sex      = ssd@meta.data[cell_id, "genotype"]
  )

# ---- Map MOFA factors to cells via Sample -----------------------------------
# MOFA sample IDs are pseudobulk-level; need to match to Seurat Sample column.
# Check naming convention — adjust if MOFA uses a different ID scheme.
cat("\nMOFA sample IDs (first 5):\n")
print(head(factor_df$sample_id, 5))
cat("Seurat Sample values (first 5 unique):\n")
print(head(unique(cell_df$Sample), 5))

# Try direct join first; if no matches, print diagnostic
cell_factors <- cell_df %>%
  left_join(factor_df, by = c("Sample" = "sample_id"))

n_matched <- sum(!is.na(cell_factors$Factor1))
cat(sprintf("\nCells with MOFA factor values: %d / %d (%.1f%%)\n",
            n_matched, nrow(cell_factors), 100 * n_matched / nrow(cell_factors)))

if (n_matched == 0) {
  cat("\n*** NO MATCHES — sample ID mismatch between MOFA and Seurat. ***\n")
  cat("Check naming: MOFA uses '", factor_df$sample_id[1],
      "' vs Seurat uses '", unique(cell_df$Sample)[1], "'\n")
  cat("You may need to adjust the join key above.\n")
  cat("Stopping here — fix the join and re-run.\n")
  quit(status = 1)
}

# ============================================================================
# PART 1: MOFA FACTOR VALUES ON UMAP
# ============================================================================
cat("\n========== MOFA Factor UMAPs ==========\n")

n_factors <- MOFA2::get_dimensions(mofa)$K
factor_cols <- paste0("Factor", seq_len(n_factors))

# ---- Shared color scale for the sex-associated factors ----------------------
# Factor1 and Factor5 are plotted on a COMMON symmetric scale so the panels are
# directly comparable side by side. Latent factor signs are arbitrary in MOFA
# (only the factor x weight product is identified), so signs are retained here
# rather than flipped, and a diverging palette centered at zero is used: warm =
# positive, cool = negative, magnitude comparable across both panels.
# Limits are computed on the pseudobulk-sample factor values (not per cell),
# so replicated cells from the same sample cannot inflate the range.
shared_factors <- intersect(c("Factor1", "Factor5"), factor_cols)

shared_lim <- max(abs(as.matrix(factor_df[, shared_factors, drop = FALSE])),
                  na.rm = TRUE)
shared_lim <- ceiling(shared_lim * 100) / 100  # round up for a clean legend
cat(sprintf("\nShared symmetric limits for %s: [%.2f, %.2f]\n",
            paste(shared_factors, collapse = " & "), -shared_lim, shared_lim))

# Returns the appropriate color scale for a given factor:
#   Factor1 / Factor5 -> shared diverging scale, fixed limits, centered at 0
#   all others        -> per-factor magma scale (free limits, as before)
factor_scale <- function(f) {
  if (f %in% shared_factors) {
    ggplot2::scale_color_gradient2(
      low      = "#2166AC",
      mid      = "grey95",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-shared_lim, shared_lim),
      oob      = scales::squish,
      na.value = "grey60",
      name     = f
    )
  } else {
    ggplot2::scale_color_viridis_c(option = "magma", na.value = "grey85",
                                   name = f)
  }
}

# Individual factor plots
factor_plots <- list()
for (f in factor_cols) {
  p <- ggplot(cell_factors %>% arrange(!is.na(.data[[f]]), .data[[f]]),
              aes(x = UMAP_1, y = UMAP_2, color = .data[[f]])) +
    geom_point(size = 0.1, alpha = 0.6) +
    factor_scale(f) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      plot.title   = element_text(hjust = 0.5, size = 13, face = "bold"),
      legend.position = "right"
    ) +
    labs(title = f, x = "UMAP 1", y = "UMAP 2")

  factor_plots[[f]] <- p

  save_dual_format(p, out_dir, paste0("umap_", tolower(f)),
                   width = 7, height = 5.5)
}

# Combined panel: all factors
combined_factors <- wrap_plots(factor_plots, ncol = 3) +
  plot_annotation(title = "MOFA Factor Values on UMAP",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
save_dual_format(combined_factors, out_dir, "umap_all_factors",
                 width = 18, height = ifelse(n_factors <= 3, 5.5, 11))

# ============================================================================
# PART 2: FACTOR VALUES SPLIT BY SEX (side-by-side)
# ============================================================================
cat("\n========== Factor UMAPs split by sex ==========\n")

# Focus on the interesting factors: Factor2 (GW), Factor1 & Factor5 (sex)
key_factors <- intersect(c("Factor1", "Factor2", "Factor5"), factor_cols)

for (f in key_factors) {
  p_sex <- ggplot(cell_factors %>% arrange(!is.na(.data[[f]]), .data[[f]]),
                  aes(x = UMAP_1, y = UMAP_2, color = .data[[f]])) +
    geom_point(size = 0.1, alpha = 0.6) +
    factor_scale(f) +
    facet_wrap(~ Sex, ncol = 2) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      strip.text   = element_text(size = 12, face = "bold"),
      plot.title   = element_text(hjust = 0.5, size = 13, face = "bold"),
      legend.position = "right"
    ) +
    labs(title = paste(f, "— by genotype"), x = "UMAP 1", y = "UMAP 2")

  save_dual_format(p_sex, out_dir, paste0("umap_", tolower(f), "_by_sex"),
                   width = 12, height = 5.5)
}

# ============================================================================
# PART 3: GENE EXPRESSION FEATURE PLOTS (Figure 5 candidates)
# ============================================================================
cat("\n========== Gene Expression UMAPs ==========\n")

# Key genes from Figure 5 gene-metabolite associations:
#   PLCG2  — Phospholipase C, OPC lipid substrate association
#   NRG1   — Oligodendrocyte differentiation signal
#   BCAS1  — Active myelination marker
#   TSHZ2  — Zinc-finger TF, indirect metabolite association
#   RPL37A — Ribosomal protein, pyrrole-2-carboxylic acid association
fig5_genes <- c("PLCG2", "NRG1", "BCAS1", "TSHZ2", "RPL37A")

# Check if expression data is available in the loaded object
has_expr <- tryCatch({
  test <- GetAssayData(ssd, assay = "RNA", layer = "data")
  nrow(test) > 1
}, error = function(e) FALSE)

if (!has_expr) {
  cat("\nSSD_meta_umap has no expression data (it's a minimal meta+UMAP object).\n")
  cat("Loading the full filtered Seurat object for gene expression...\n")

  seurat_full_path <- file.path(PROJECT_ROOT, "data", "SSD_RNAonly_filtered.rds")
  if (!file.exists(seurat_full_path)) {
    cat("  Full Seurat not found at:", seurat_full_path, "\n")
    cat("  Skipping gene expression UMAPs.\n")
    cat("  Set seurat_full_path to the correct location and re-run.\n")
    has_expr <- FALSE
  } else {
    cat("  This may take a few minutes...\n")
    ssd_full <- readRDS(seurat_full_path)
    DefaultAssay(ssd_full) <- "RNA"

    # Check which genes are present
    genes_present <- intersect(fig5_genes, rownames(ssd_full))
    genes_missing <- setdiff(fig5_genes, genes_present)
    if (length(genes_missing) > 0)
      cat("  Genes not found in Seurat:", paste(genes_missing, collapse = ", "), "\n")
    cat("  Genes available:", paste(genes_present, collapse = ", "), "\n")

    has_expr <- length(genes_present) > 0
  }
} else {
  ssd_full <- ssd
  genes_present <- intersect(fig5_genes, rownames(ssd_full))
  genes_missing <- setdiff(fig5_genes, genes_present)
  if (length(genes_missing) > 0)
    cat("Genes not found:", paste(genes_missing, collapse = ", "), "\n")
  cat("Genes available:", paste(genes_present, collapse = ", "), "\n")
}

if (has_expr && length(genes_present) > 0) {

  # Extract expression for target genes
  expr_mat <- GetAssayData(ssd_full, assay = "RNA", layer = "data")[genes_present, , drop = FALSE]

  # Build plot df using the same UMAP coordinates
  # (ssd_full should have the same umap.rna reduction)
  umap_full <- as.data.frame(Embeddings(ssd_full, reduction = "umap.rna"))
  colnames(umap_full) <- c("UMAP_1", "UMAP_2")

  gene_plots <- list()

  for (gene in genes_present) {
    expr_vec <- expr_mat[gene, ]

    plot_df <- umap_full %>%
      tibble::rownames_to_column("cell_id") %>%
      mutate(
        expression = expr_vec[cell_id],
        celltype   = ssd_full@meta.data[cell_id, "celltype"]
      )

    # Order: zero-expression cells plotted first (grey), expressing on top
    plot_df <- plot_df %>% arrange(expression)

    p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = expression)) +
      geom_point(size = 0.1, alpha = 0.6) +
      scale_color_viridis_c(
        option  = "viridis",
        name    = "log-norm\nexpr",
        na.value = "grey85",
        limits  = c(0, quantile(expr_vec[expr_vec > 0], 0.99, na.rm = TRUE))
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid   = element_blank(),
        axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        plot.title   = element_text(hjust = 0.5, size = 13, face = "bold.italic"),
        legend.position = "right"
      ) +
      labs(title = gene, x = "UMAP 1", y = "UMAP 2")

    gene_plots[[gene]] <- p

    save_dual_format(p, out_dir, paste0("umap_gene_", gene),
                     width = 7, height = 5.5)
  }

  # Combined panel
  combined_genes <- wrap_plots(gene_plots, ncol = 3) +
    plot_annotation(
      title = "Key Gene Expression — Figure 5 Candidates",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  save_dual_format(combined_genes, out_dir, "umap_genes_combined",
                   width = 18, height = ifelse(length(genes_present) <= 3, 5.5, 11))
}

# ============================================================================
# PART 4: MOFA VARIANCE EXPLAINED HEATMAP (cleaner than weight lollipops)
# ============================================================================
cat("\n========== MOFA Variance Explained ==========\n")

p_var <- MOFA2::plot_variance_explained(mofa, max_r2 = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
save_dual_format(p_var, out_dir, "mofa_variance_explained",
                 width = 10, height = 5)

# ============================================================================
# Summary
# ============================================================================
cat("\n========================================\n")
cat("All plots saved to:", out_dir, "\n")
cat("Files:\n")
list.files(out_dir, pattern = "\\.png$") %>% cat(sep = "\n  ")
cat("\n\nDone.\n")
