# pseudobulk_functions.R
#
# Shared functions for pseudobulk + limma-voom DE analysis of snRNA-seq data.
# Sourced by hormone_gene_limma_sexseparate.Rmd and sex-combined scripts.
#
# Design:
#   build_pseudobulk()         -> pure aggregation; no hormone/phenotype knowledge.
#                                 Returns counts (samples x genes), logcounts, metadata.
#   merge_hormone_metadata()   -> joins hormone phenotype data onto pb$metadata and
#                                 re-aligns counts/logcounts to the matched samples.
#   run_voom_lm*()             -> three model variants (predictor only, + cell_count,
#                                 + cell_count + GW). Sex-separate analysis.
#   run_voom_lm*_sex()         -> three sex-combined variants; same covariates as above
#                                 plus Sex as an additive covariate. Designed for use
#                                 with analyze_hormone_combined() so that sex differences
#                                 are accounted for when estimating the hormone effect
#                                 across all samples together.
#   analyze_hormone()          -> per (cell_type, Sex) split, filter, run model.
#                                 Use for sex-separate analyses.
#   analyze_hormone_combined() -> per cell_type split only; Sex covariate enters
#                                 the model rather than the grouping. Use for
#                                 sex-combined analyses.
#   plot_top_genes_from_subsets(), summarize_sig_genes(), plot_sig_gene_summary()
#                              -> plotting / summary helpers shared by both modes.
#
# Note on the gene filter inside analyze_hormone() and analyze_hormone_combined():
#   Uses edgeR::filterByExpr() with default min.count = 10 and min.prop = 0.7.
#   Filter is library-size-aware (CPM scale), runs per cell_type (x Sex) group
#   using only that group's samples. See Chen, Lun & Smyth 2016 (F1000Research)
#   and Law et al. 2018 (RNAseq123) for the canonical limma-edgeR filter pattern.


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
#'                          pseudobulk side is `str_extract(sample_id, pb_sample_pattern)`
#'                          rather than the full `sample_id`. Use this when the
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

# Helper: backtick-quote variable names that aren't syntactically valid.
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
  formula_vars <- .safe_vars(c(predictor, "GW"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_cc <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  formula_vars <- .safe_vars(c(predictor, "cell_count"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_cc_gw <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"GW" %in% colnames(metadata)) metadata$GW <- NA
  formula_vars <- .safe_vars(c(predictor, "cell_count", "GW"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}


# ---- 3b. limma-voom model wrappers (sex-combined) ---------------------------
#
# These are identical to the sex-separate wrappers above, except that Sex is
# included as an additive covariate. Including Sex controls for sex-related
# expression differences when estimating the hormone coefficient across the
# pooled (male + female) sample. See Law et al. 2018 (F1000Research) for the
# general limma-edgeR workflow and the rationale for blocking or covarying
# on known batch/biological factors before testing.

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

run_voom_lm_sex_gw <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  formula_vars <- .safe_vars(c(predictor, "Sex", "GW"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_cc_sex <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  formula_vars <- .safe_vars(c(predictor, "cell_count", "Sex"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}

run_voom_lm_cc_gw_sex <- function(counts, metadata, predictor) {
  stopifnot(ncol(counts) == nrow(metadata))
  if (!"GW" %in% colnames(metadata)) metadata$GW <- NA
  formula_vars <- .safe_vars(c(predictor, "cell_count", "GW", "Sex"))
  design <- model.matrix(reformulate(formula_vars), data = metadata)
  v   <- limma::voom(counts, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  predictor_col <- colnames(design)[2]
  limma::topTable(fit, coef = predictor_col, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene")
}


# ---- 4. Per-group analysis (cell_type x Sex) — sex-separate -----------------

# NOTE: the gene-filter line below (colSums + [, keep_genes]) is preserved
# exactly as in the previous version. It is currently under review — see the
# Rmd / project memory before changing it.

analyze_hormone <- function(pb_counts, metadata_merged, hormone,
                            overall_min_samples = 8,
                            min_count = 10,
                            min_prop = 0.7,
                            model = run_voom_lm) {

  split_ids <- split(rownames(metadata_merged),
                     interaction(metadata_merged$cell_type, metadata_merged$Sex,
                                 drop = TRUE))

  results_by_group <- purrr::map(names(split_ids), function(group_name) {
    ids          <- split_ids[[group_name]]
    counts_sub   <- t(pb_counts)[, ids, drop = FALSE]
    metadata_sub <- metadata_merged[ids, , drop = FALSE]

    cat("\nGroup:", group_name, "\n")
    cat("Initial number of samples:", nrow(metadata_sub), "\n")

    if (nrow(metadata_sub) < overall_min_samples) {
      cat("Group removed: insufficient samples", nrow(metadata_sub), "\n")
      return(NULL)
    }

    # Gene-level filtering via edgeR::filterByExpr.
    # Operates on genes x samples; library-size-aware (CPM scale) so the
    # effective threshold scales with each group's pseudobulk depth.
    keep_genes <- edgeR::filterByExpr(counts_sub,
                                      min.count = min_count,
                                      min.prop  = min_prop)
    counts_sub <- counts_sub[keep_genes, , drop = FALSE]
    cat("Genes passing filterByExpr:", sum(keep_genes), "\n")

    if (ncol(counts_sub) == 0) {
      cat("Group removed: no genes passed the filtering threshold\n")
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


# ---- 4b. Per-group analysis (cell_type only) — sex-combined -----------------
#
# Mirrors analyze_hormone() with two key differences:
#   1. Splits only by cell_type (not cell_type x Sex); Sex enters the model
#      as a covariate via the run_voom_lm_*_sex() wrappers above. This is the
#      standard approach for controlling a known biological factor when pooling
#      groups — see Ritchie et al. 2015 (limma paper) and Law et al. 2018.
#   2. result$sex is set to "combined" rather than a single-sex value, to
#      distinguish sex-combined output rows in downstream summary tables.
#
# Note: nrow() check after filterByExpr (not ncol()) corrects the direction of
# the subsetting check — filterByExpr removes genes (rows), not samples (cols).

analyze_hormone_combined <- function(pb_counts, metadata_merged, hormone,
                                     overall_min_samples = 8,
                                     min_count = 10,
                                     min_prop = 0.7,
                                     model = run_voom_lm_sex) {

  split_ids <- split(rownames(metadata_merged),
                     metadata_merged$cell_type)   # cell_type only; no Sex split

  results_by_group <- purrr::map(names(split_ids), function(group_name) {
    ids          <- split_ids[[group_name]]
    counts_sub   <- t(pb_counts)[, ids, drop = FALSE]
    metadata_sub <- metadata_merged[ids, , drop = FALSE]

    cat("\nGroup:", group_name, "\n")
    cat("Initial number of samples:", nrow(metadata_sub), "\n")

    if (nrow(metadata_sub) < overall_min_samples) {
      cat("Group removed: insufficient samples", nrow(metadata_sub), "\n")
      return(NULL)
    }

    keep_genes <- edgeR::filterByExpr(counts_sub,
                                      min.count = min_count,
                                      min.prop  = min_prop)
    counts_sub <- counts_sub[keep_genes, , drop = FALSE]
    cat("Genes passing filterByExpr:", sum(keep_genes), "\n")

    if (nrow(counts_sub) == 0) {
      cat("Group removed: no genes passed the filtering threshold\n")
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
      cat("Group removed: only one sex present after hormone filtering (",
          unique(final_metadata_sub$Sex), ")\n")
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

plot_top_genes_from_subsets <- function(by_group_data, output_dir, hormone,
                                        fdr_cutoff = 0.1,
                                        logfc_cutoff = 0.25,
                                        top_n = 15) {
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

    counts_cpm <- edgeR::cpm(counts_sub, log = FALSE)
    counts_log <- t(log1p(counts_cpm))           # samples x genes

    # Apply the same FDR + logFC filter as summarize_sig_genes(), so the number
    # of trajectories plotted matches the bar-chart count for this group.
    # `top_n` is a safety cap (defaults to Inf = no cap).
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
        message("    Gene ", gene, " not in counts — skipping")
        return(NULL)
      }

      expr_vals <- counts_log[, gene]
      if (all(is.na(expr_vals)) || all(expr_vals == 0)) {
        message("    Gene ", gene, " has all NA or zero expression — skipping")
        return(NULL)
      }

      message("    Plotting gene: ", gene,
              " (adj. p-value: ", format(adj_pval, digits = 3), ")")
      flush.console()

      df <- metadata_sub %>% dplyr::mutate(expr = expr_vals)

      ggplot(df, aes(x = !!sym(hormone), y = expr,
                     color = cell_count, shape = Sex)) +
        geom_point(size = 2, alpha = 0.8) +
        geom_smooth(method = "lm", se = TRUE) +
        scale_color_viridis_c() +
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

    ggsave(file.path(output_dir, paste0(clean_name, "_top_genes.png")),
           plot = combined_plot,
           width = 10,
           height = 6 + ceiling(length(plot_list) / 2),
           dpi = 300)
    message("  ✅ Saved plot for ", clean_name)
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
