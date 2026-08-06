#!/usr/bin/env Rscript
# figure5_panels.R
# ---------------------------------------------------------------------------
# Standalone script to generate manuscript Figure 5 panels while the limma
# rerun is writing parquets. Reads the hits CSV directly (no parquet needed)
# and the cached data objects for scatter plots.
#
# Run from src/gene-metabolite/:
#   Rscript figure5_panels.R
#
# After the rerun finishes, fold these chunks back into
# metabolite_gene_associations.qmd (the code is identical — just copy-paste
# the labeled sections).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(svglite)
  library(edgeR)
  library(patchwork)
})

source("pseudobulk_functions.R")
source("metabolite_functions.R")
set.seed(123)

# ---- Config ----------------------------------------------------------------
METHOD       <- "log2_na"
cache_dir    <- "../../data/cache"
metab_path   <- "../../results/untargeted/batch2_peak_area_clean.csv"
hormone_xlsx <- "../../doc/targeted/targeted_hormones.xlsx"
hits_path    <- "../../results/gene-metabolite/csv-log2_na/metabolite_gene_hits.csv"
out_dir      <- "../../results/gene-metabolite"
fig_dir      <- file.path(out_dir, "figure5-panels")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data objects (same as QMD setup, no parquet) ---------------------
pb_base     <- readRDS(file.path(cache_dir, "pb_base.rds"))
batch_props <- readRDS(file.path(cache_dir, "batch_props_sample.rds"))

targeted_comparison_raw <- readxl::read_xlsx(hormone_xlsx)
targeted_comparison_raw$E2 <- targeted_comparison_raw$`17B-E2 (ng/mL)`
targeted_comparison_raw$TT <- targeted_comparison_raw$`TT (ng/mL)`
targeted_comparison_raw$P4 <- targeted_comparison_raw$`P4 (ng/mL)`

pb       <- build_metab_covariates(pb_base, targeted_comparison_raw, batch_props,
                                   min_cells = 30)
metab_mat <- read_metabolite_matrix(metab_path)
metab_tx  <- transform_metabolite_matrix(metab_mat, method = METHOD)

rn        <- rownames(pb$counts)
row_ct    <- sub("^[^.]+\\.(.*)$", "\\1", rn)
metab_col <- metab_sample_to_pb(colnames(metab_tx), sub("^([^.]+)\\..*$", "\\1", rn))

# ---- Load hits directly from CSV ------------------------------------------
hits <- read_csv(hits_path, show_col_types = FALSE)
cat(sprintf("Loaded %s hits from %s\n", format(nrow(hits), big.mark = ","), hits_path))

# ---- Cell-type recurrence (needed for bar chart) ---------------------------
ct_recur <- hits |>
  group_by(cell_type) |>
  summarise(n_associations = n(),
            n_metabolites  = n_distinct(Compound.ID),
            n_genes        = n_distinct(gene), .groups = "drop") |>
  arrange(desc(n_associations))


# ============================================================================
# FIGURE 5A — Associations by cell type (grouped bar)
# ============================================================================
cat("\n--- Figure 5A: Associations by cell type ---\n")

ct_long <- ct_recur |>
  pivot_longer(cols = c(n_associations, n_metabolites, n_genes),
               names_to  = "metric",
               values_to = "count") |>
  mutate(metric = factor(metric,
    levels = c("n_associations", "n_metabolites", "n_genes"),
    labels = c("Associations", "Distinct metabolites", "Distinct genes")))

p_assoc <- ggplot(ct_long,
                  aes(x = reorder(cell_type, count, .fun = max),
                      y = count, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Associations"         = "#4E79A7",
                                "Distinct metabolites" = "#F28E2B",
                                "Distinct genes"       = "#59A14F"),
                    name = NULL) +
  labs(x = NULL, y = "Count",
       title = "Gene\u2013metabolite associations by cell type") +
  theme_cowplot(12) +
  theme(legend.position = "top")

save_dual_format(p_assoc, fig_dir, "fig5a_associations_by_celltype",
                 width = 10, height = 6)
print(p_assoc)
cat("Saved fig5a\n")


# ============================================================================
# Spotlight scatter helper
# ============================================================================

make_spotlight_scatter <- function(cid, gene, ct, meta_name,
                                   label = NULL, subtitle = NULL) {
  idx      <- which(row_ct == ct)
  pos      <- match(metab_col[idx], colnames(metab_tx))
  present  <- !is.na(pos) & gene %in% colnames(pb$counts)
  rows_use <- idx[present]; pos_p <- pos[present]
  if (length(rows_use) < 3) {
    message("Too few samples for ", gene, " ~ ", meta_name, " in ", ct)
    return(NULL)
  }

  # Gene expression: log1p(CPM) computed per cell type
  expr <- as.numeric(log1p(
    edgeR::cpm(t(pb$counts[rows_use, , drop = FALSE]), log = FALSE)[gene, ]
  ))
  md <- pb$metadata[rn[rows_use], ]

  df <- tibble(
    metab       = as.numeric(metab_tx[cid, pos_p]),
    expr        = expr,
    Sex         = md$Sex,
    GW          = md$GW,
    batch1_frac = md$batch1_frac
  ) |> filter(!is.na(metab))

  if (nrow(df) < 4) {
    message("Too few non-NA for ", gene, " ~ ", meta_name, " in ", ct)
    return(NULL)
  }

  # Partial residual: remove covariate effects, keep metabolite slope
  covs   <- c("GW", "Sex", "batch1_frac")
  keepc  <- covs[vapply(covs, function(cv) length(unique(df[[cv]])) >= 2, logical(1))]
  fit    <- lm(reformulate(c("metab", keepc), response = "expr"), data = df)
  df$partial <- resid(fit) + coef(fit)[["metab"]] * df$metab

  # Stats annotation
  hit_row <- hits |> filter(Compound.ID == cid, gene == !!gene, cell_type == ct)
  lfc_txt <- if (nrow(hit_row) > 0) sprintf("logFC = %.2f", hit_row$logFC[1]) else ""
  fdr_txt <- if (nrow(hit_row) > 0) sprintf("FDR = %.3f",   hit_row$FDR_global[1]) else ""
  stat_label <- paste(lfc_txt, fdr_txt, sep = "  |  ")

  if (is.null(label)) label <- paste0(gene, " ~ ", meta_name, "  (", ct, ")")

  p <- ggplot(df, aes(metab, partial, color = GW)) +
    geom_point(aes(shape = Sex), size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9,
                color = "grey30") +
    scale_color_viridis_c(name = "GW") +
    scale_shape_manual(values = c("M" = 16, "F" = 17)) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = stat_label, size = 3.5, color = "grey40") +
    labs(x = paste0(meta_name, " [", METHOD, " transformed]"),
         y = "Gene expression [log1p(CPM), adjusted]",
         title    = label,
         subtitle = subtitle) +
    theme_cowplot(13) +
    theme(plot.subtitle = element_text(color = "grey40", size = 10))

  p
}


# ============================================================================
# FIGURE 5C — OPC lipid spotlight
# ============================================================================
cat("\n--- Figure 5C: OPC lipid spotlight ---\n")

# BCAS1 ~ phospholipid (active myelination marker)
bcas1_row <- hits |> filter(gene == "BCAS1", cell_type == "OPC",
                             grepl("Glycerophospholipid", Class))
if (nrow(bcas1_row) > 0) {
  p_bcas1 <- make_spotlight_scatter(
    cid = bcas1_row$Compound.ID[1], gene = "BCAS1", ct = "OPC",
    meta_name = bcas1_row$Name[1],
    subtitle  = "Active-myelination marker with phospholipid"
  )
  if (!is.null(p_bcas1)) {
    save_dual_format(p_bcas1, fig_dir, "fig5c_OPC_BCAS1_phospholipid", width = 8, height = 6)
    cat("Saved BCAS1 spotlight\n")
  }
}

# PLCG2 ~ Palmitic Acid (phospholipase + fatty acid substrate)
plcg2_row <- hits |> filter(gene == "PLCG2", Name == "Palmitic Acid", cell_type == "OPC")
if (nrow(plcg2_row) > 0) {
  p_plcg2 <- make_spotlight_scatter(
    cid = plcg2_row$Compound.ID[1], gene = "PLCG2", ct = "OPC",
    meta_name = "Palmitic Acid",
    subtitle  = "Phospholipase C with fatty acid substrate"
  )
  if (!is.null(p_plcg2)) {
    save_dual_format(p_plcg2, fig_dir, "fig5c_OPC_PLCG2_palmitic", width = 8, height = 6)
    cat("Saved PLCG2 spotlight\n")
  }
}

# NRG1 ~ phospholipid (OPC differentiation signal)
nrg1_row <- hits |> filter(gene == "NRG1", cell_type == "OPC")
if (nrow(nrg1_row) > 0) {
  p_nrg1 <- make_spotlight_scatter(
    cid = nrg1_row$Compound.ID[1], gene = "NRG1", ct = "OPC",
    meta_name = nrg1_row$Name[1],
    subtitle  = "Oligodendrocyte differentiation signal with phospholipid"
  )
  if (!is.null(p_nrg1)) {
    save_dual_format(p_nrg1, fig_dir, "fig5c_OPC_NRG1_phospholipid", width = 8, height = 6)
    cat("Saved NRG1 spotlight\n")
  }
}


# ============================================================================
# INTERPRETATION SPOTLIGHTS
# ============================================================================
cat("\n--- Interpretation spotlights ---\n")

# ---- Direct: PLCG2 ~ Cer(D15:2_7:0) in OPC ----
# Phospholipase C directly cleaves lipid substrates
direct_row <- hits |> filter(gene == "PLCG2", cell_type == "OPC", grepl("Cer", Name))
if (nrow(direct_row) > 0) {
  p_direct <- make_spotlight_scatter(
    cid = direct_row$Compound.ID[1], gene = "PLCG2", ct = "OPC",
    meta_name = direct_row$Name[1],
    label    = "Direct: PLCG2 ~ ceramide (OPC)",
    subtitle = "Phospholipase C with lipid substrate \u2014 expected association"
  )
  if (!is.null(p_direct)) {
    save_dual_format(p_direct, fig_dir, "interp_direct_PLCG2_ceramide", width = 8, height = 6)
    cat("Saved direct example\n")
  }
}

# ---- Ribosomal: RPL37A ~ Pyrrole-2-Carboxylic Acid in EN-Non-IT-Immature ----
# 73 RPL/RPS genes associate with this metabolite; reflects coordinated
# translational machinery and metabolic activity
ribo_row <- hits |> filter(gene == "RPL37A", Name == "Pyrrole-2-Carboxylic Acid")
if (nrow(ribo_row) > 0) {
  p_ribo <- make_spotlight_scatter(
    cid = ribo_row$Compound.ID[1], gene = "RPL37A", ct = "EN-Non-IT-Immature",
    meta_name = "Pyrrole-2-Carboxylic Acid",
    label    = "Ribosomal: RPL37A ~ Pyrrole-2-Carboxylic Acid (EN-Non-IT-Immature)",
    subtitle = "1 of 73 ribosomal genes associated with this metabolite"
  )
  if (!is.null(p_ribo)) {
    save_dual_format(p_ribo, fig_dir, "interp_ribosomal_RPL37A_pyrrole", width = 8, height = 6)
    cat("Saved ribosomal example\n")
  }
}

# ---- Indirect: TSHZ2 ~ Palmitic Acid in IN-MGE-SST ----
# Zinc-finger TF for cortical neuron specification; no enzymatic link to
# fatty acid metabolism
tf_row <- hits |> filter(gene == "TSHZ2", Name == "Palmitic Acid", cell_type == "IN-MGE-SST")
if (nrow(tf_row) > 0) {
  p_tf <- make_spotlight_scatter(
    cid = tf_row$Compound.ID[1], gene = "TSHZ2", ct = "IN-MGE-SST",
    meta_name = "Palmitic Acid",
    label    = "Indirect: TSHZ2 ~ Palmitic Acid (IN-MGE-SST)",
    subtitle = "Zinc-finger TF \u2014 no enzymatic metabolite link"
  )
  if (!is.null(p_tf)) {
    save_dual_format(p_tf, fig_dir, "interp_indirect_TSHZ2_palmitic", width = 8, height = 6)
    cat("Saved TF/indirect example\n")
  }
}

cat("\n--- All Figure 5 panels saved to", fig_dir, "---\n")
