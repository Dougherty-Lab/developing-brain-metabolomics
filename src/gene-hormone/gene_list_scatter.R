# =============================================================================
# gene_list_scatter.R
# =============================================================================
#
# PURPOSE
#   Plot a manually specified gene list against current hormone + expression
#   data to check whether hormone–expression structure is preserved across
#   dataset versions.
#
#   For each hormone × gene: scatter of hormone concentration (x) vs
#   log1p(CPM) expression (y), one panel per cell type.
#   Annotated with Pearson r and raw p-value (exploratory only — no FDR).
#
# OUTPUT
#   results/gene-hormone/gene_list_scatter/
#     TT_gene_list_scatter.pdf
#     E2_gene_list_scatter.pdf
#     P4_gene_list_scatter.pdf
#
# USAGE
#   Rscript gene_list_scatter.R
#   (run from src/gene-hormone/ — same working directory as the primary QMDs)
#
# =============================================================================

setwd("src/gene-hormone/")
# ---- 0. Manual gene list (edit here) ----------------------------------------

GENE_LIST <- c(
  "CADM2",
  "LSAMP",
  "MARCH1",
  "SEMA3E",
  "SULF1",
  "DIAPH3",
  "LINC02232",
  "SYT9",
  "PDZRN4",
  "AL162493.1",
  "AC091078.1", 
  "RBFOX1", 
  "HAS2",
  "AC02305.1",
  "KCNH5",
  "MGAT4C",
  "FAM106A",
  "PBX1",
  "DLG2",
  "SOX5",
  "AC013265.1",
  "LINC02503",
  "LINC01446",
  "LINC01242",
  "DDX11",
  "AC107419.1"
)


# ---- 1. Packages & source ---------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(edgeR)
  library(patchwork)
  library(svglite)
})

source("pseudobulk_functions.R")
set.seed(123)

# ---- 2. Load data -----------------------------------------------------------

pb_base      <- readRDS("../../data/cache/pb_base.rds")
SSD_data     <- readRDS("../../data/cache/SSD_meta_umap.rds")

targeted_raw <- read_xlsx("../../doc/targeted/targeted_hormones.xlsx") |>
  dplyr::mutate(
    E2 = `17B-E2 (ng/mL)`,
    TT = `TT (ng/mL)`,
    P4 = `P4 (ng/mL)`
  )

# ---- 3. Outlier removal (identical to primary QMDs) -------------------------

detect_outliers <- function(data, variable) {
  limit <- quantile(data[[variable]], 0.75, na.rm = TRUE) +
           1.5 * IQR(data[[variable]], na.rm = TRUE)
  dplyr::filter(data, .data[[variable]] <= limit)
}

targeted_TT <- detect_outliers(targeted_raw, "TT")
targeted_E2 <- detect_outliers(targeted_raw, "E2")
targeted_P4 <- detect_outliers(targeted_raw, "P4")

# ---- 4. Build pb objects (identical to primary QMDs) ------------------------

pb_sample_pattern <- "tissue\\d+"

build_pb <- function(hormone_df) {
  pb <- merge_hormone_metadata(
    pb_base,
    hormone_df,
    pb_sample_pattern = pb_sample_pattern,
    min_cells         = 30
  )
  add_batch_fraction(pb, SSD_data)
}

pb_TT <- build_pb(targeted_TT)
pb_E2 <- build_pb(targeted_E2)
pb_P4 <- build_pb(targeted_P4)

# ---- 5. Core plot function --------------------------------------------------

#' One-page scatter for a single gene across all cell types.
#'
#' @param pb         pb object from merge_hormone_metadata() + add_batch_fraction()
#' @param gene       Gene name. Must be in colnames(pb$counts).
#' @param hormone    Column name of x-axis hormone: "TT", "E2", or "P4".
#' @param male_only  Logical. TRUE restricts to Sex == "M" (TT analysis).
#' @param min_obs    Minimum non-NA observations per cell type to include a panel.
#'
#' @return A patchwork plot (one page), or NULL if no valid cell types.
plot_gene_scatter <- function(pb, gene, hormone, male_only = FALSE,
                              min_obs = 5) {

  meta <- pb$metadata

  if (male_only) {
    meta <- meta[meta$Sex == "M", , drop = FALSE]
  }

  # pb$counts is samples x genes; cpm() expects genes x samples → transpose twice
  counts_sub <- pb$counts[rownames(meta), , drop = FALSE]    # samples x genes
  cpm_mat    <- t(log1p(edgeR::cpm(t(counts_sub), log = FALSE)))  # samples x genes

  if (!gene %in% colnames(cpm_mat)) {
    message("  Gene not found in counts: ", gene)
    return(NULL)
  }

  cell_types <- sort(unique(meta$cell_type))

  panel_list <- purrr::map(cell_types, function(ct) {

    idx <- which(meta$cell_type == ct)
    if (length(idx) < min_obs) return(NULL)

    df <- data.frame(
      hormone_val = as.numeric(meta[[hormone]][idx]),
      expr        = cpm_mat[idx, gene],
      cell_count  = meta$cell_count[idx],
      Sex         = meta$Sex[idx],
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$hormone_val) & !is.na(df$expr), , drop = FALSE]
    if (nrow(df) < min_obs) return(NULL)

    # Pearson r — appropriate since primary model is linear (limma-voom)
    ct_test <- cor.test(df$hormone_val, df$expr, method = "pearson")
    ann <- sprintf("r = %.2f  p = %.3f  n = %d",
                   ct_test$estimate, ct_test$p.value, nrow(df))

    # For male-only (TT), Sex is constant so drop the shape aesthetic
    if (male_only) {
      aes_call <- aes(x = hormone_val, y = expr, color = cell_count)
    } else {
      aes_call <- aes(x = hormone_val, y = expr,
                      color = cell_count, shape = Sex)
    }

    p <- ggplot(df, aes_call) +
      geom_point(size = 2, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE,
                  color = "#E74C3C", fill = "#E74C3C",
                  alpha = 0.15, linewidth = 0.8,
                  show.legend = FALSE) +
      scale_color_viridis_c(name = "cells") +
      labs(
        title    = ct,
        subtitle = ann,
        x        = paste0(hormone, " (ng/mL)"),
        y        = "log1p(CPM)"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        legend.position  = "none",
        plot.title       = element_text(size = 9,  face = "bold"),
        plot.subtitle    = element_text(size = 7.5, color = "grey40"),
        axis.title       = element_text(size = 8),
        panel.border     = element_rect(color = "grey80", fill = NA,
                                        linewidth = 0.4)
      )

    if (!male_only) {
      p <- p + scale_shape_manual(values = c("M" = 16, "F" = 17))
    }

    p
  }) |> purrr::compact()

  if (length(panel_list) == 0) return(NULL)

  n_panels <- length(panel_list)
  ncols    <- min(4L, n_panels)

  patchwork::wrap_plots(panel_list, ncol = ncols) +
    patchwork::plot_annotation(
      title    = paste0(gene, "  \u2014  ", hormone),
      subtitle = sprintf(
        "%d cell type%s with \u2265 %d valid observations%s",
        n_panels,
        if (n_panels == 1) "" else "s",
        min_obs,
        if (male_only) "  |  Male samples only" else ""
      ),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "grey40")
      )
    )
}

# ---- 6. Generate PDFs per hormone -------------------------------------------

output_dir <- "../../results/gene-hormone/gene_list_scatter"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

hormone_config <- list(
  TT = list(pb = pb_TT, male_only = TRUE),
  E2 = list(pb = pb_E2, male_only = FALSE),
  P4 = list(pb = pb_P4, male_only = FALSE)
)

for (h in names(hormone_config)) {

  pb        <- hormone_config[[h]]$pb
  male_only <- hormone_config[[h]]$male_only

  valid_genes <- intersect(GENE_LIST, colnames(pb$counts))
  missing     <- setdiff(GENE_LIST, colnames(pb$counts))

  if (length(missing) > 0)
    message(h, ": genes not found in counts — ",
            paste(missing, collapse = ", "))

  if (length(valid_genes) == 0) {
    message(h, ": no valid genes to plot — skipping\n")
    next
  }

  pdf_path <- file.path(output_dir, paste0(h, "_gene_list_scatter.pdf"))
  pdf(pdf_path, width = 14, height = 10)

  n_plotted <- 0L
  for (gene in valid_genes) {
    message(h, " / ", gene, " ...")
    p <- plot_gene_scatter(pb, gene, hormone = h, male_only = male_only)
    if (!is.null(p)) {
      print(p)
      n_plotted <- n_plotted + 1L
    } else {
      message("  No valid cell types — page skipped")
    }
  }

  dev.off()
  message(h, ": ", n_plotted, " / ", length(valid_genes),
          " genes plotted  →  ", pdf_path, "\n")
}

message("Done.")
