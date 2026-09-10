# geneset_enrichment.R
# ---------------------------------------------------------------------------
# Gene-set enrichment analysis (fgsea) on metabolite-associated gene rankings
# from limma-voom results.
#
# Inputs:  Association CSVs, gene sets from build_gene_sets.R
# Outputs: Enrichment CSVs and figures in results/gene-metabolite/geneset/
#
# Upstream:  metabolite_gene_associations.Rmd, build_gene_sets.R
# Downstream: None
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(cowplot)
  library(svglite)
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
source(file.path(root, "src/gene-hormone/pseudobulk_functions.R"))

set.seed(123)

# ---- Config ----------------------------------------------------------------
METHOD      <- "log2_na"          # "int" | "zscore" | "zscore_trim" | "log2_na"
suffix      <- if (METHOD == "int") "" else paste0("-", METHOD)
parquet_dir <- if (METHOD == "int") {
  file.path(root, "results/gene-metabolite/parquet")
} else {
  sprintf(file.path(root, "results/gene-metabolite/parquet-%s"), METHOD)
}

hits_path     <- sprintf(file.path(root, "results/gene-metabolite/csv%s/metabolite_gene_hits.csv"), suffix)
gene_set_path <- file.path(root, "doc/gene_lists/disease_gene_sets.csv")
cache_dir     <- file.path(root, "data/cache")
out_dir       <- file.path(root, "results/gene-metabolite/geneset-enrichment")
csv_dir       <- out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Minimum distinct hit genes for a cell type to be tested.
MIN_HIT_GENES_PER_CT <- 20

FDR_ALPHA        <- 0.05   # significance is AFTER BH correction
REBUILD_UNIVERSE <- FALSE  # TRUE to force a re-scan of the parquet archive

# ---- Load hits and gene sets -----------------------------------------------
hits <- read_csv(hits_path, show_col_types = FALSE)

gene_sets <- read_csv(gene_set_path, show_col_types = FALSE)

set_meta <- gene_sets |>
  distinct(set_id, set_label, condition, source, level, set_type)

# One reported set per condition in plot
# Where a condition draws on more than one publication (SCZ, Epilepsy) the union set is reported
# where it draws on one (NDD, T2D, IBD) the paper set already IS the condition set. 
# Per-paper tests written to the CSVs 
reported_set_ids <- set_meta |>
  dplyr::group_by(condition) |>
  dplyr::filter(if (any(level == "condition")) level == "condition"
                else level == "paper") |>
  dplyr::ungroup() |>
  dplyr::pull(set_id)

adjust_split <- function(df) {
  df |>
    dplyr::mutate(reported = set_id %in% reported_set_ids) |>
    dplyr::group_by(reported) |>
    dplyr::mutate(FDR = p.adjust(p, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::mutate(significant = !is.na(FDR) & FDR < FDR_ALPHA)
}

cat(sprintf("Hits: %d associations, %d unique genes, %d cell types\n",
            nrow(hits), n_distinct(hits$gene), n_distinct(hits$cell_type)))
cat(sprintf("Gene sets: %d (%d paper-level, %d condition-level)\n",
            nrow(set_meta), sum(set_meta$level == "paper"),
            sum(set_meta$level == "condition")))

# ---- Universes -------------------------------------------------------------
# Global universe is shared with sfari_enrichment.R 
universe_cache <- file.path(cache_dir,
                            sprintf("tested_gene_universe_%s.rds", METHOD))
# Per-cell-type universe is specific to this analysis.
ct_universe_cache <- file.path(cache_dir,
                               sprintf("tested_gene_universe_by_celltype_%s.rds", METHOD))

if (!REBUILD_UNIVERSE && file.exists(universe_cache)) {
  background_genes <- readRDS(universe_cache)
} else {
  cat("Scanning parquet archive for global tested-gene universe ...\n")
  background_genes <- open_dataset(parquet_dir) |>
    dplyr::distinct(gene) |> collect() |> dplyr::pull(gene) |> unique() |> sort()
  saveRDS(background_genes, universe_cache)
}
cat(sprintf("Global tested-gene universe: %d\n", length(background_genes)))

if (!REBUILD_UNIVERSE && file.exists(ct_universe_cache)) {
  ct_universe <- readRDS(ct_universe_cache)
} else {
  cat("Scanning parquet archive for per-cell-type tested-gene universe ...\n")
  ct_universe <- open_dataset(parquet_dir) |>
    dplyr::distinct(cell_type, gene) |> collect()
  saveRDS(ct_universe, ct_universe_cache)
}
cat(sprintf("Per-cell-type universe: %d cell types, %s gene-by-cell-type pairs\n",
            n_distinct(ct_universe$cell_type),
            format(nrow(ct_universe), big.mark = ",")))

# ---- Gene set coverage -----------------------------------------------------
coverage <- gene_sets |>
  group_by(set_id) |>
  dplyr::summarise(n_genes   = n_distinct(gene_symbol),
            n_tested  = n_distinct(gene_symbol[gene_symbol %in% background_genes]),
            .groups   = "drop") |>
  dplyr::mutate(pct_tested = 100 * n_tested / n_genes) |>
  left_join(set_meta, by = "set_id") |>
  dplyr::arrange(set_type, condition, level)

cat("\nGene set coverage in the tested universe:\n")
print(as.data.frame(coverage |>
        dplyr::select(set_label, level, set_type, n_genes, n_tested, pct_tested)),
      right = FALSE, digits = 3)

# Symbol-nomenclature warning. 
# Several published sets use older GENCODE symbols and clone-based lincRNA names (RP11-*). 
low_cov <- coverage |> dplyr::filter(pct_tested < 60)
if (nrow(low_cov)) {
  warning("Sets with <60% of genes in the tested universe -- check for ",
          "deprecated gene symbols before interpreting: ",
          paste(low_cov$set_id, collapse = ", "))
}

set_list <- split(gene_sets$gene_symbol, gene_sets$set_id) |> purrr::map(unique)

# ---- Fisher helper ---------------------------------------------------------
# One-sided (greater): the hypothesis is enrichment, not depletion. 
# The CI is taken from the two-sided test 
fisher_one <- function(fg, bg, set_genes) {
  fg <- unique(fg)
  bg <- unique(bg)
  fg <- intersect(fg, bg)             # foreground must be inside background
  sg <- intersect(set_genes, bg)      # set restricted to what was testable

  fg_in  <- sum(fg %in% sg)
  fg_out <- length(fg) - fg_in
  bg_only <- setdiff(bg, fg)
  bg_in  <- sum(bg_only %in% sg)
  bg_out <- length(bg_only) - bg_in

  m <- matrix(c(fg_in, fg_out, bg_in, bg_out), nrow = 2, byrow = TRUE)

  if (fg_in == 0 || length(sg) == 0) {
    return(tibble(n_fg = length(fg), n_set_tested = length(sg),
                  fg_in_set = fg_in, bg_in_set = bg_in,
                  OR = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                  p = NA_real_))
  }

  ft  <- fisher.test(m, alternative = "greater")
  ft2 <- fisher.test(m)

  tibble(n_fg = length(fg), n_set_tested = length(sg),
         fg_in_set = fg_in, bg_in_set = bg_in,
         OR = unname(ft$estimate),
         ci_lo = ft2$conf.int[1], ci_hi = ft2$conf.int[2],
         p = ft$p.value)
}

# ---- Analysis 1: pooled ----------------------------------------------------
hit_genes_all <- unique(hits$gene)

pooled <- imap_dfr(set_list, function(sg, sid) {
  fisher_one(hit_genes_all, background_genes, sg) |> dplyr::mutate(set_id = sid)
}) |>
  left_join(set_meta, by = "set_id") |>
  dplyr::mutate(pct_fg = 100 * fg_in_set / n_fg,
         pct_bg = 100 * bg_in_set / (length(background_genes) - n_fg))

# BH across all sets tested in this analysis.
pooled <- pooled |>
  adjust_split() |>
  dplyr::arrange(set_type, dplyr::desc(reported), p)

cat("\n=== Analysis 1: pooled (hit genes vs tested universe) ===\n")
print(as.data.frame(pooled |>
        dplyr::select(set_label, level, set_type, n_set_tested, fg_in_set,
               OR, ci_lo, ci_hi, p, FDR, significant)),
      right = FALSE, digits = 3)

write_csv(pooled, file.path(csv_dir, "geneset_enrichment_pooled.csv"))

# ---- Panel B: pooled forest, one row per condition -------------------------

condition_level <- pooled |> dplyr::filter(reported)

forest_pooled <- condition_level |>
  dplyr::filter(is.finite(OR), OR > 0, is.finite(ci_hi)) |>
  dplyr::mutate(
    panel = ifelse(set_type == "control", "Negative controls",
                   "Neurodevelopmental"),
    sig   = ifelse(significant, sprintf("FDR < %.2f", FDR_ALPHA), "n.s."),
    lab   = sprintf("%d/%d", fg_in_set, n_set_tested),
    condition = fct_reorder(condition, OR))

p_forest_pooled <- ggplot(forest_pooled,
                          aes(x = OR, y = condition, colour = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18,
                 linewidth = 0.6) +
  geom_point(size = 3) +
  geom_text(aes(label = lab), vjust = -1.1, size = 2.9, show.legend = FALSE) +
  scale_x_log10() +
  scale_colour_manual(
    values = setNames(c("#B2182B", "grey55"),
                      c(sprintf("FDR < %.2f", FDR_ALPHA), "n.s.")),
    name = NULL) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Odds ratio (log scale)", y = NULL,
       title = "Disease gene set enrichment among metabolite-associated genes",
       subtitle = "Pooled across cell types; background = all tested genes",
       caption = "Points labelled set genes among hits / set genes tested. Bars are 95% CI.") +
  theme_cowplot(12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"))

save_dual_format(p_forest_pooled, out_dir, "geneset_forest_pooled",
                 width = 9, height = 6)
print(p_forest_pooled)

cat("\nCondition-level pooled results (plotted):\n")
print(as.data.frame(condition_level |>
        dplyr::select(condition, set_label, level, n_set_tested, fg_in_set,
                      OR, ci_lo, ci_hi, p, FDR, significant)),
      right = FALSE, digits = 3)

# ---- Analysis 1b: covariate-matched permutation null -----------------------
# Fisher's test compares set genes among hits against ALL tested genes. 

# For each set, draw N random gene sets of the same size, sampling WITHIN
# strata defined by expression decile x number of cell types 
N_PERM_MATCH <- 10000

pb_base   <- readRDS(file.path(cache_dir, "pb_base.rds"))
cpm       <- pb_base$counts / rowSums(pb_base$counts) * 1e6
mean_expr <- log2(colMeans(cpm) + 1)

n_ct_tested <- ct_universe |> dplyr::count(gene, name = "n_ct_tested")

match_df <- tibble(gene = background_genes) |>
  dplyr::left_join(n_ct_tested, by = "gene") |>
  dplyr::mutate(mean_expr   = unname(mean_expr[gene]),
                n_ct_tested = tidyr::replace_na(n_ct_tested, 0L),
                is_hit      = gene %in% hit_genes_all) |>
  dplyr::filter(!is.na(mean_expr)) |>
  dplyr::mutate(expr_decile = dplyr::ntile(mean_expr, 10),
                stratum     = paste(expr_decile, n_ct_tested, sep = "_"))

pool <- split(match_df$is_hit, match_df$stratum)

cat(sprintf("\n=== Analysis 1b: detectability-matched null (%s permutations) ===\n",
            format(N_PERM_MATCH, big.mark = ",")))
cat(sprintf("Matching strata: %d (expression decile x cell types tested)\n",
            length(pool)))

set.seed(123)
matched <- imap_dfr(set_list, function(sg, sid) {
  in_set <- match_df |> dplyr::filter(gene %in% sg)
  if (nrow(in_set) == 0)
    return(tibble(set_id = sid, n_set_matched = 0L, observed = 0L,
                  null_mean = NA_real_, null_sd = NA_real_,
                  fold = NA_real_, p_matched = NA_real_))

  obs      <- sum(in_set$is_hit)
  strata_n <- in_set |> dplyr::count(stratum, name = "k")

  null <- replicate(N_PERM_MATCH, {
    sum(map2_int(strata_n$stratum, strata_n$k, function(s, k) {
      v <- pool[[s]]
      if (is.null(v) || length(v) == 0) return(0L)
      sum(sample(v, min(k, length(v))))
    }))
  })

  tibble(set_id = sid, n_set_matched = nrow(in_set), observed = obs,
         null_mean = mean(null), null_sd = sd(null),
         fold = obs / mean(null),
         p_matched = (sum(null >= obs) + 1) / (N_PERM_MATCH + 1))
}) |>
  dplyr::left_join(set_meta, by = "set_id") |>
  dplyr::mutate(reported = set_id %in% reported_set_ids) |>
  dplyr::group_by(reported) |>
  dplyr::mutate(FDR_matched = p.adjust(p_matched, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::mutate(significant_matched = !is.na(FDR_matched) &
                  FDR_matched < FDR_ALPHA) |>
  dplyr::arrange(set_type, dplyr::desc(reported), p_matched)

print(as.data.frame(matched |>
        dplyr::filter(reported) |>
        dplyr::select(condition, set_label, set_type, n_set_matched, observed,
                      null_mean, fold, p_matched, FDR_matched,
                      significant_matched)),
      right = FALSE, digits = 3)

write_csv(matched, file.path(csv_dir, "geneset_matched_permutation.csv"))

# Side-by-side: unmatched odds ratio against matched fold enrichment. 
compare <- pooled |>
  dplyr::select(set_id, condition, set_label, set_type, reported,
                OR, p_unmatched = p, FDR_unmatched = FDR) |>
  dplyr::left_join(matched |> dplyr::select(set_id, fold, p_matched,
                                            FDR_matched, observed, null_mean),
                   by = "set_id")

write_csv(compare, file.path(csv_dir, "geneset_matched_vs_unmatched.csv"))

comp_plot <- compare |>
  dplyr::filter(reported, is.finite(OR), is.finite(fold)) |>
  tidyr::pivot_longer(c(OR, fold), names_to = "test", values_to = "estimate") |>
  dplyr::mutate(
    test  = factor(test, levels = c("OR", "fold"),
                   labels = c("Unmatched (odds ratio)",
                              "Detectability-matched (fold)")),
    panel = ifelse(set_type == "control", "Negative controls",
                   "Neurodevelopmental"),
    sig   = dplyr::case_when(
      test == "Unmatched (odds ratio)"        & FDR_unmatched < FDR_ALPHA ~ TRUE,
      test == "Detectability-matched (fold)"  & FDR_matched   < FDR_ALPHA ~ TRUE,
      TRUE ~ FALSE))

p_compare <- ggplot(comp_plot, aes(x = estimate, y = condition,
                                   colour = test, shape = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.55)) +
  scale_x_log10() +
  scale_colour_manual(values = c("Unmatched (odds ratio)" = "#B2182B",
                                 "Detectability-matched (fold)" = "#2166AC"),
                      name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     name = paste0("FDR < ", FDR_ALPHA)) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Enrichment (log scale)", y = NULL,
       title = "Enrichment before and after matching on detectability",
       subtitle = sprintf("Matched null: %s permutations within expression-decile x cell-types-tested strata",
                          format(N_PERM_MATCH, big.mark = ",")),
       caption = "A large gap between the two estimates indicates enrichment driven by gene detectability rather than gene set membership.") +
  theme_cowplot(12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"))

save_dual_format(p_compare, out_dir, "geneset_matched_vs_unmatched",
                 width = 9, height = 6)
print(p_compare)

# ---- Analysis 2: per cell type ---------------------------------------------
hits_by_ct <- hits |> distinct(cell_type, gene)

ct_counts <- hits_by_ct |> dplyr::count(cell_type, name = "n_hit_genes") |>
  dplyr::mutate(tested = n_hit_genes >= MIN_HIT_GENES_PER_CT) |>
  dplyr::arrange(desc(n_hit_genes))

cat(sprintf("\n=== Analysis 2: per cell type (floor = %d hit genes) ===\n",
            MIN_HIT_GENES_PER_CT))
print(as.data.frame(ct_counts), right = FALSE)

ct_keep <- ct_counts |> dplyr::filter(tested) |> pull(cell_type)

per_ct <- expand_grid(cell_type = ct_keep, set_id = names(set_list)) |>
  pmap_dfr(function(cell_type, set_id) {
    fg <- hits_by_ct$gene[hits_by_ct$cell_type == cell_type]
    bg <- ct_universe$gene[ct_universe$cell_type == cell_type]
    fisher_one(fg, bg, set_list[[set_id]]) |>
      dplyr::mutate(cell_type = cell_type, set_id = set_id, n_bg = length(unique(bg)))
  }) |>
  left_join(set_meta, by = "set_id")

# BH across the whole set x cell-type grid -- these are all tests of the same
# question, so correcting within each set (or each cell type) would understate
# the burden.
per_ct <- per_ct |>
  adjust_split() |>
  dplyr::arrange(set_type, condition, cell_type)

cat(sprintf("\nGrid: %d cell types x %d sets = %d tests (%d reported, %d supplementary per-paper)\n",
            length(ct_keep), length(set_list), nrow(per_ct),
            sum(per_ct$reported), sum(!per_ct$reported)))
cat(sprintf("Significant at FDR < %.2f: %d reported, %d supplementary\n",
            FDR_ALPHA,
            sum(per_ct$significant & per_ct$reported),
            sum(per_ct$significant & !per_ct$reported)))
print(as.data.frame(per_ct |> dplyr::filter(significant, reported) |>
        dplyr::select(set_label, cell_type, n_set_tested, fg_in_set, OR, p, FDR)),
      right = FALSE, digits = 3)

write_csv(per_ct,    file.path(csv_dir, "geneset_enrichment_by_celltype.csv"))
write_csv(ct_counts, file.path(csv_dir, "geneset_celltype_inclusion.csv"))

# ---- Panel A: set x cell type dot matrix -----------------------------------
# Colour = log2 OR, size = -log10 FDR. Paper-level sets only; condition-level
# unions are their own panel so a reader never counts the same gene twice.
dot_df <- per_ct |>
  dplyr::filter(reported, !is.na(OR), is.finite(OR)) |>
  dplyr::mutate(log2_OR   = log2(OR),
         neglog10  = -log10(pmax(FDR, 1e-300)),
         set_label = fct_reorder(set_label, as.integer(factor(condition))),
         panel     = ifelse(set_type == "control",
                            "Negative controls", "Neurodevelopmental"))

p_dot <- ggplot(dot_df, aes(x = cell_type, y = set_label)) +
  geom_point(aes(size = neglog10, fill = log2_OR, colour = significant),
             shape = 21, stroke = 0.7) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = expression(log[2]~OR)) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"),
                      name = paste0("FDR < ", FDR_ALPHA)) +
  scale_size_continuous(range = c(1.5, 8),
                        name = expression(-log[10]~FDR)) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title = "Disease gene set enrichment among metabolite-associated genes",
       subtitle = sprintf("Background = genes tested in that cell type; cell types with >= %d hit genes",
                          MIN_HIT_GENES_PER_CT)) +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey92"),
        panel.grid.major = element_line(colour = "grey93"))

save_dual_format(p_dot, out_dir, "geneset_dotmatrix_by_celltype",
                 width = 11, height = 7)

# ---- Panel C: OR vs hit-gene count -----------------------------------------
power_df <- per_ct |>
  dplyr::filter(reported, !is.na(OR), is.finite(OR)) |>
  left_join(ct_counts, by = "cell_type")

p_power <- ggplot(power_df, aes(x = n_hit_genes, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = set_type, shape = significant), size = 2.4,
             alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(neuro = "#B2182B", control = "grey45"),
                      name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     name = paste0("FDR < ", FDR_ALPHA)) +
  labs(x = "Distinct hit genes in cell type (log scale)",
       y = "Odds ratio (log scale)",
       title = "Enrichment vs. detection power") +
  theme_cowplot(12)

save_dual_format(p_power, out_dir, "geneset_or_vs_power", width = 8, height = 6)

# ---- Panel D: gene set overlap ---------------------------------------------
sids <- reported_set_ids
jac <- expand_grid(a = sids, b = sids) |>
  pmap_dfr(function(a, b) {
    ga <- set_list[[a]]; gb <- set_list[[b]]
    tibble(a = a, b = b,
           jaccard  = length(intersect(ga, gb)) / length(union(ga, gb)),
           n_shared = length(intersect(ga, gb)))
  }) |>
  left_join(set_meta |> dplyr::select(a = set_id, lab_a = set_label), by = "a") |>
  left_join(set_meta |> dplyr::select(b = set_id, lab_b = set_label), by = "b")

p_jac <- ggplot(jac, aes(x = lab_a, y = lab_b, fill = jaccard)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(a == b, "", n_shared)), size = 3) +
  scale_fill_gradient(low = "white", high = "#4A1486", name = "Jaccard") +
  labs(x = NULL, y = NULL, title = "Gene set overlap (cell labels = shared genes)") +
  theme_cowplot(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_dual_format(p_jac, out_dir, "geneset_overlap_jaccard", width = 9, height = 7)
write_csv(jac |> dplyr::filter(a != b), file.path(csv_dir, "geneset_overlap.csv"))

# ---- Which genes and which metabolites? ------------------------------------
# Two levels, each computed for the full dataset and then for the highlighted
# cell type x set pairs:
#   gene level       -- which disease-set genes carry associations
#   metabolite level -- what metabolites they associate with, and whether it differs from
#                       the metabolites associated with non-set genes
#
# Highlights are DERIVED from the results by default rather than hardcoded, so
# the figure tracks the analysis if thresholds or data change. 
# Set HIGHLIGHT_OVERRIDE to pin specific pairs instead.
HIGHLIGHT_OVERRIDE <- NULL
# e.g. tibble::tribble(~cell_type,     ~set_id,
#                      "IN-MGE-SST",   "scz_union",
#                      "IN-MGE-SST",   "ndd_fu_dd477",
#                      "OPC",          "epilepsy_union",
#                      "EN-L6b",       "ndd_fu_dd477")

# Long table of every hit that falls in a disease set, one row per
# gene x metabolite x cell type x set.
hit_set_long <- hits |>
  dplyr::select(dplyr::any_of(c("gene", "cell_type", "Compound.ID", "Name",
                                "Class", "Super.Class", "logFC"))) |>
  dplyr::inner_join(gene_sets |> dplyr::select(set_id, gene_symbol),
                    by = c("gene" = "gene_symbol"),
                    relationship = "many-to-many") |>
  dplyr::left_join(set_meta, by = "set_id")

write_csv(hit_set_long, file.path(csv_dir, "geneset_hit_associations_long.csv"))

cat(sprintf("\nHit associations falling in a disease gene set: %s rows, %d genes\n",
            format(nrow(hit_set_long), big.mark = ","),
            dplyr::n_distinct(hit_set_long$gene)))

# ---- Panel F: gene x cell type, full dataset -------------------------------
# Paper-level neuro sets only. 
gene_ct <- hit_set_long |>
  dplyr::filter(set_type == "neuro", set_id %in% reported_set_ids) |>
  dplyr::group_by(condition, gene, cell_type) |>
  dplyr::summarise(n_metabolites = dplyr::n_distinct(Compound.ID),
                   .groups = "drop")

if (nrow(gene_ct) > 0) {
  gene_order <- gene_ct |>
    dplyr::group_by(gene) |>
    dplyr::summarise(total = sum(n_metabolites), .groups = "drop") |>
    dplyr::arrange(total) |> dplyr::pull(gene)

  p_gene_ct <- gene_ct |>
    dplyr::mutate(gene = factor(gene, levels = gene_order)) |>
    ggplot(aes(x = cell_type, y = gene, fill = n_metabolites)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = n_metabolites), size = 2.4, colour = "grey15") +
    facet_grid(condition ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "#DEEBF7", high = "#08519C", name = "Metabolites") +
    labs(x = NULL, y = NULL,
         title = "Disease gene set members with metabolite associations",
         subtitle = "All cell types; genes ordered by total distinct metabolites") +
    theme_cowplot(11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 6.5),
          strip.background = element_rect(fill = "grey92"))

  save_dual_format(p_gene_ct, out_dir, "geneset_gene_by_celltype_heatmap",
                   width = 9,
                   height = max(6, 0.13 * dplyr::n_distinct(gene_ct$gene)))
  print(p_gene_ct)
  write_csv(gene_ct, file.path(csv_dir, "geneset_gene_by_celltype.csv"))
}

# ---- Panel G: metabolite class composition, full dataset -------------------

class_col <- intersect(c("Class", "Super.Class", "Subclass"), names(hits))

if (length(class_col) > 0) {
  cc <- class_col[1]
  neuro_genes <- gene_sets |>
    dplyr::filter(set_type == "neuro") |> dplyr::pull(gene_symbol) |> unique()

  class_df <- hits |>
    dplyr::mutate(group = ifelse(gene %in% neuro_genes,
                                 "Disease-set gene", "Other hit gene"),
                  metab_class = tidyr::replace_na(as.character(.data[[cc]]),
                                                  "Unannotated")) |>
    dplyr::distinct(gene, Compound.ID, group, metab_class) |>
    dplyr::count(group, metab_class, name = "n") |>
    dplyr::group_by(group) |>
    dplyr::mutate(pct = 100 * n / sum(n)) |>
    dplyr::ungroup()

  keep <- class_df |> dplyr::group_by(metab_class) |>
    dplyr::summarise(mx = max(pct), .groups = "drop") |>
    dplyr::filter(mx >= 2) |> dplyr::pull(metab_class)

  class_plot_df <- class_df |>
    dplyr::mutate(metab_class = ifelse(metab_class %in% keep, metab_class,
                                       "Other (<2%)")) |>
    dplyr::group_by(group, metab_class) |>
    dplyr::summarise(n = sum(n), pct = sum(pct), .groups = "drop")

  tab <- class_df |>
    dplyr::select(group, metab_class, n) |>
    tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0) |>
    tibble::column_to_rownames("metab_class") |> as.matrix()

  chisq_res <- if (all(dim(tab) >= 2) && sum(tab) > 0)
    suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000)) else NULL

  if (!is.null(chisq_res))
    cat(sprintf("Metabolite class composition, disease-set vs other hit genes: chi-square p = %.4g (simulated)\n",
                chisq_res$p.value))

  p_class <- ggplot(class_plot_df,
                    aes(x = pct, y = fct_reorder(metab_class, pct), fill = group)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.9) +
    scale_fill_manual(values = c("Disease-set gene" = "#08519C",
                                 "Other hit gene" = "#BAB0AC"), name = NULL) +
    labs(x = "Percent of gene-metabolite pairs", y = NULL,
         title = "Metabolite classes associated with disease-set vs other hit genes",
         subtitle = if (!is.null(chisq_res))
           sprintf("Chi-square p = %.3g (simulated); classes below 2%% collapsed",
                   chisq_res$p.value) else "Classes below 2% collapsed") +
    theme_cowplot(11) + theme(legend.position = "bottom")

  save_dual_format(p_class, out_dir, "geneset_metabolite_class_composition",
                   width = 9, height = 7)
  print(p_class)
  write_csv(class_df, file.path(csv_dir, "geneset_metabolite_class_composition.csv"))
} else {
  cat("\nMetabolite class panel skipped: no Class/Super.Class column in the hits CSV.\n")
}

# ---- Panels H-I: highlighted cell type x set pairs -------------------------

highlights <- if (!is.null(HIGHLIGHT_OVERRIDE)) {
  HIGHLIGHT_OVERRIDE
} else {
  per_ct |>
    dplyr::filter(significant, reported, set_type == "neuro") |>
    dplyr::select(cell_type, set_id)
}

if (nrow(highlights) == 0) {
  cat("\nNo significant neuro cell-type x set pairs to highlight.\n")
} else {
  cat(sprintf("\nHighlighting %d significant cell type x set pairs:\n",
              nrow(highlights)))
  print(as.data.frame(highlights |>
          dplyr::left_join(set_meta |> dplyr::select(set_id, set_label),
                           by = "set_id")), right = FALSE)

  metab_label <- if ("Name" %in% names(hits)) "Name" else "Compound.ID"

  hl_long <- highlights |>
    dplyr::inner_join(hit_set_long, by = c("cell_type", "set_id")) |>
    dplyr::mutate(metabolite = .data[[metab_label]],
                  metabolite = ifelse(is.na(metabolite), Compound.ID,
                                      as.character(metabolite)))

  write_csv(hl_long, file.path(csv_dir, "geneset_highlight_associations.csv"))

  for (i in seq_len(nrow(highlights))) {
    ct  <- highlights$cell_type[i]
    sid <- highlights$set_id[i]
    lab <- set_meta$set_label[match(sid, set_meta$set_id)]

    d <- hl_long |> dplyr::filter(cell_type == ct, set_id == sid)
    if (nrow(d) == 0) next

    # Truncate long IUPAC-style metabolite names, which are otherwise wide
    # enough to squeeze the panel to nothing.
    d <- d |>
      dplyr::mutate(metabolite = ifelse(nchar(metabolite) > 42,
                                        paste0(substr(metabolite, 1, 39), "..."),
                                        metabolite))

    has_lfc <- "logFC" %in% names(d)

    p_hl <- ggplot(d, aes(x = metabolite, y = gene)) +
      {
        if (has_lfc) geom_tile(aes(fill = logFC), colour = "white", linewidth = 0.4)
        else geom_tile(fill = "#08519C", colour = "white", linewidth = 0.4)
      } +
      {
        if (has_lfc) scale_fill_gradient2(low = "#2166AC", mid = "white",
                                          high = "#B2182B", midpoint = 0,
                                          name = expression(log[2]~FC))
      } +
      labs(x = NULL, y = NULL,
           title = sprintf("%s \u2014 %s", ct, lab),
           subtitle = sprintf("%d genes, %d metabolites, %d associations",
                              dplyr::n_distinct(d$gene),
                              dplyr::n_distinct(d$Compound.ID), nrow(d))) +
      theme_cowplot(11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            axis.text.y = element_text(size = 8))

    fn <- sprintf("geneset_highlight_%s_%s",
                  gsub("[^A-Za-z0-9]+", "-", ct), sid)
    save_dual_format(p_hl, out_dir, fn,
                     width  = max(7, 0.30 * dplyr::n_distinct(d$metabolite) + 3),
                     height = max(4, 0.22 * dplyr::n_distinct(d$gene) + 2))
    print(p_hl)
  }
}

cat("\nSession info:\n"); print(sessionInfo())

# ---- AI assistance disclosure ------------------------------------------------
# Code in this script was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.

# ---- session info ------------------------------------------------------------
cat("\n\n---- Session Info ----\n")
print(sessionInfo())
