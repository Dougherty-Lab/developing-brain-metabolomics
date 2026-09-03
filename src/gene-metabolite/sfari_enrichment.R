#!/usr/bin/env Rscript
# sfari_enrichment.R
# ---------------------------------------------------------------------------
# Two clean enrichment tests for SFARI autism risk genes:
#
#   Test 1: Are SFARI genes enriched among genes with ANY significant
#           gene-metabolite association?
#           Foreground = hit genes, Background = TESTED gene universe
#
#   Test 2: Among hit genes, are SFARI genes more likely to be associated
#           with >= N distinct metabolites?
#           Foreground = genes associated with >= N metabolites
#           Background = hit genes
#
# Each test: Fisher's exact + Wilcoxon rank-sum + SFARI score stratification.
#
# Background note (Test 1): the background is the set of genes that actually
# entered a limma model -- i.e. every gene appearing in the parquet archive --
# NOT colnames(pb_base$counts). Genes dropped by the per-cell-type expression
# filter had zero opportunity to become a hit; including them inflates the
# non-hit/non-SFARI cell and biases the OR upward. 
# Universe = union over cell types (a gene tested in >=1 cell type), which
# matches the foreground definition ("a hit in >=1 cell type").
#
# Run from src/gene-metabolite/:
#   Rscript sfari_enrichment.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(cowplot)
  library(svglite)
})

source("pseudobulk_functions.R")

# ---- Config ----------------------------------------------------------------
# METHOD keys the parquet archive, the hits CSV, and the cached gene universe
# together so a run can never mix transforms.
METHOD      <- "log2_na"                       # "int" | "zscore" | "zscore_trim" | "log2_na"
suffix      <- if (METHOD == "int") "" else paste0("-", METHOD)
parquet_dir <- if (METHOD == "int") {
  "../../results/gene-metabolite/parquet"
} else {
  sprintf("../../results/gene-metabolite/parquet-%s", METHOD)
}

hits_path  <- sprintf("../../results/gene-metabolite/csv%s/metabolite_gene_hits.csv", suffix)
sfari_path <- "../../doc/gene_lists/SFARI-Gene_genes_07-12-2026release_08-13-2026export.csv"
cache_dir  <- "../../data/cache"
out_dir    <- "../../results/gene-metabolite/sfari-enrichment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_METAB_RECURRENCE <- 5
MIN_HIT_GENES_PER_CT <- 20      # cell-type floor, matches geneset_enrichment.R
FDR_ALPHA            <- 0.05    # significance judged after BH, not on raw p
REBUILD_UNIVERSE     <- FALSE   # TRUE to force a re-scan of the parquet archive

# ---- Load data -------------------------------------------------------------
hits <- read_csv(hits_path, show_col_types = FALSE)

# ---- Tested-gene universe (Test 1 background) ------------------------------
# Union of genes that entered >=1 limma model, scanned once from the parquet
# archive and cached. distinct() on a single column is a cheap Arrow query --
# only the `gene` column is read -- but it still touches every file, so cache it.
universe_cache <- file.path(cache_dir,
                            sprintf("tested_gene_universe_%s.rds", METHOD))

if (!REBUILD_UNIVERSE && file.exists(universe_cache)) {
  background_genes <- readRDS(universe_cache)
  cat(sprintf("Tested-gene universe (cached): %d\n", length(background_genes)))
} else {
  cat("Scanning parquet archive for tested-gene universe ...\n")
  background_genes <- open_dataset(parquet_dir) |>
    dplyr::distinct(gene) |>
    collect() |>
    dplyr::pull(gene) |>
    unique() |>
    sort()
  saveRDS(background_genes, universe_cache)
  cat(sprintf("Tested-gene universe (scanned, cached to %s): %d\n",
              universe_cache, length(background_genes)))
}

# Sanity check against the unfiltered pseudobulk matrix, for reporting only.
pb_base <- readRDS(file.path(cache_dir, "pb_base.rds"))
cat(sprintf("  (pseudobulk matrix has %d genes; %d dropped by expression filtering)\n",
            ncol(pb_base$counts), ncol(pb_base$counts) - length(background_genes)))
cat(sprintf("Gene-metabolite hits: %d associations, %d unique genes\n",
            nrow(hits), n_distinct(hits$gene)))

sfari <- read_csv(sfari_path, show_col_types = FALSE) |>
  rename(symbol = `gene-symbol`, ensembl = `ensembl-id`,
         score = `gene-score`, gene_name = `gene-name`)

sfari_symbols <- sfari$symbol
sfari_ensembl <- sfari$ensembl

# ---- SFARI lookup -----------------------------------------------------------
sfari_by_sym <- sfari |> dplyr::filter(!is.na(symbol),  symbol  != "") |>
  dplyr::select(symbol, score)  |> deframe()
sfari_by_ens <- sfari |> dplyr::filter(!is.na(ensembl), ensembl != "") |>
  dplyr::select(ensembl, score) |> deframe()

is_sfari_gene <- function(g) g %in% sfari_symbols | g %in% sfari_ensembl

get_sfari_score <- function(g) {
  case_when(
    g %in% names(sfari_by_sym) ~ sfari_by_sym[g],
    g %in% names(sfari_by_ens) ~ sfari_by_ens[g],
    TRUE ~ NA_real_
  )
}

# ---- Gene recurrence from hits ---------------------------------------------
hit_recur <- hits |>
  group_by(gene) |>
  summarise(n_metabolites  = n_distinct(Compound.ID),
            n_cell_types   = n_distinct(cell_type),
            n_associations = n(), .groups = "drop") |>
  mutate(is_sfari    = is_sfari_gene(gene),
         sfari_score = get_sfari_score(gene),
         is_recurrent      = n_metabolites >= MIN_METAB_RECURRENCE) |>
  arrange(desc(n_metabolites))

hit_genes <- hit_recur$gene

# ---- Helper: Fisher + Wilcoxon + stratified --------------------------------
run_fisher <- function(fg, bg, label_fg, label_bg) {
  fg_sfari    <- sum(is_sfari_gene(fg))
  fg_nonsfari <- length(fg) - fg_sfari
  bg_only     <- setdiff(bg, fg)
  bg_sfari    <- sum(is_sfari_gene(bg_only))
  bg_nonsfari <- length(bg_only) - bg_sfari

  contingency <- matrix(
    c(fg_sfari, fg_nonsfari, bg_sfari, bg_nonsfari),
    nrow = 2, byrow = TRUE,
    dimnames = list(c(label_fg, paste0("Non-", label_fg)),
                    c("SFARI", "Non-SFARI"))
  )
  cat("\nContingency table:\n")
  print(contingency)

  ft <- fisher.test(contingency, alternative = "greater")
  cat(sprintf("\nFisher's exact (one-sided, greater):\n  OR = %.2f, 95%% CI [%.2f, Inf], p = %.2e\n",
              ft$estimate, ft$conf.int[1], ft$p.value))
  cat(sprintf("  %s: %d/%d SFARI (%.1f%%)\n  %s: %d/%d SFARI (%.1f%%)\n",
              label_fg, fg_sfari, length(fg), 100 * fg_sfari / length(fg),
              label_bg, sum(is_sfari_gene(bg)), length(bg),
              100 * sum(is_sfari_gene(bg)) / length(bg)))
  ft
}

# ---- Helper: SFARI score-stratified Fisher ---------------------------------
# For each SFARI score s, ask whether score-s genes are over-represented in the
# foreground relative to NON-SFARI genes, with both restricted to `bg`. Using
# non-SFARI (rather than all-other) as the comparator keeps every score's test
# against the same reference group, so ORs are comparable across scores.
#
# `bg` MUST be the tested/eligible universe: a SFARI gene that never entered a
# model is not a failed foreground, it is an unobserved one, and counting it as
# a non-foreground SFARI gene deflates the OR.
#
# p is one-sided (greater), matching run_fisher(); the CI comes from the
# two-sided test because the one-sided interval is [lo, Inf) and cannot be drawn.
run_score_strata <- function(fg, bg, test_label, fg_label) {
  bg     <- unique(bg)
  in_fg  <- bg %in% unique(fg)
  is_s   <- is_sfari_gene(bg)
  scores <- get_sfari_score(bg)

  ns_fg    <- sum(!is_s &  in_fg)
  ns_nonfg <- sum(!is_s & !in_fg)

  one_row <- function(sel, label, n_pool) {
    s_fg    <- sum(sel &  in_fg)
    s_nonfg <- sum(sel & !in_fg)
    mat <- matrix(c(s_fg, ns_fg, s_nonfg, ns_nonfg), nrow = 2)
    ft1 <- fisher.test(mat, alternative = "greater")
    ft2 <- fisher.test(mat)                       # two-sided, for the CI only
    tibble(test = test_label, stratum = label,
           n_tested = n_pool, n_fg = s_fg,
           pct_fg = 100 * s_fg / max(n_pool, 1),
           OR = unname(ft1$estimate),
           ci_lo = ft2$conf.int[1], ci_hi = ft2$conf.int[2],
           p = ft1$p.value)
  }

  rows <- purrr::map_dfr(
    sort(unique(scores[!is.na(scores)])),
    ~ one_row(!is.na(scores) & scores == .x, paste("Score", .x),
              sum(!is.na(scores) & scores == .x))
  )

  # Unscored (syndromic-only) SFARI genes and an all-SFARI reference row.
  n_unscored <- sum(is_s & is.na(scores))
  if (n_unscored > 0)
    rows <- bind_rows(rows, one_row(is_s & is.na(scores), "Unscored", n_unscored))
  rows <- bind_rows(rows, one_row(is_s, "All SFARI", sum(is_s)))

  cat(sprintf("\n--- SFARI score stratification (%s) ---\n", test_label))
  cat(sprintf("  Comparator: %d non-SFARI genes (%d %s)\n",
              ns_fg + ns_nonfg, ns_fg, fg_label))
  for (i in seq_len(nrow(rows)))
    cat(sprintf("  %-10s %3d/%4d %s (%.1f%%), OR=%.2f [%.2f, %.2f], p=%.2e\n",
                rows$stratum[i], rows$n_fg[i], rows$n_tested[i], fg_label,
                rows$pct_fg[i], rows$OR[i], rows$ci_lo[i], rows$ci_hi[i],
                rows$p[i]))
  rows
}

# ---- TEST 1: SFARI enrichment among hit genes ------------------------------
cat("\n##########################################################\n")
cat("TEST 1: Are SFARI genes enriched among significant genes?\n")
cat("  Foreground: genes with any hit (", length(hit_genes), ")\n")
cat("  Background: tested-gene universe (", length(background_genes), ")\n")
cat("##########################################################\n")

# Every hit gene must be in the tested universe by construction; a non-empty
# setdiff means the hits CSV and the parquet archive came from different runs.
orphan_hits <- setdiff(hit_genes, background_genes)
if (length(orphan_hits) > 0)
  warning(sprintf("%d hit genes absent from the tested universe (METHOD/run mismatch?): %s",
                  length(orphan_hits),
                  paste(head(orphan_hits, 5), collapse = ", ")))

fisher_t1 <- run_fisher(hit_genes, background_genes, "Hit", "Background")
strata_t1 <- run_score_strata(hit_genes, background_genes, "Test 1", "hit")

# Wilcoxon: do SFARI genes have more associations than non-SFARI (among all genes)?
full_df <- tibble(gene = background_genes) |>
  left_join(hit_recur |> dplyr::select(gene, n_metabolites), by = "gene") |>
  mutate(n_metabolites = replace_na(n_metabolites, 0),
         is_sfari = is_sfari_gene(gene),
         has_hit  = n_metabolites > 0)

wilcox_t1 <- wilcox.test(has_hit ~ is_sfari, data = full_df |>
                           mutate(has_hit = as.numeric(has_hit)),
                         alternative = "greater")
cat(sprintf("\nWilcoxon (SFARI more likely to have hits): p = %.4f\n",
            wilcox_t1$p.value))

# ---- TEST 2: SFARI genes and multi-metabolite association ------------------
cat("\n##########################################################\n")
cat("TEST 2: Among hit genes, are SFARI genes associated with more metabolites?\n")
cat("  Foreground: genes associated with >=", MIN_METAB_RECURRENCE, "metabolites (",
    sum(hit_recur$is_recurrent), ")\n")
cat("  Background: all hit genes (", nrow(hit_recur), ")\n")
cat("##########################################################\n")

recurrent_genes <- hit_recur |> filter(is_recurrent) |> pull(gene)
fisher_t2 <- run_fisher(recurrent_genes, hit_genes, "Multi-metabolite", "Hit")

# Wilcoxon: among hit genes, do SFARI genes have higher recurrence?
wilcox_t2 <- wilcox.test(n_metabolites ~ is_sfari, data = hit_recur,
                         alternative = "greater")
cat(sprintf("\nWilcoxon (SFARI higher recurrence among hits): W = %.0f, p = %.4f\n",
            wilcox_t2$statistic, wilcox_t2$p.value))

hit_recur |>
  group_by(is_sfari) |>
  summarise(n = n(), median_metab = median(n_metabolites),
            mean_metab = mean(n_metabolites),
            n_recurrent = sum(is_recurrent), .groups = "drop") |>
  print()

# ---- SFARI score stratification (Test 2) -----------------------------------
cat("\n--- SFARI score stratification (among hit genes) ---\n")
sfari_hits <- hit_recur |> filter(is_sfari, !is.na(sfari_score))

sfari_hits |>
  group_by(sfari_score) |>
  summarise(n = n(), n_recurrent = sum(is_recurrent),
            pct_recurrent = 100 * mean(is_recurrent),
            median_metab = median(n_metabolites), .groups = "drop") |>
  print()

strata_t2 <- run_score_strata(recurrent_genes, hit_genes, "Test 2",
                              sprintf(">=%d metabolites", MIN_METAB_RECURRENCE))

# Both tests share the helper, so the score axis is defined identically:
# Test 1 asks "does this score reach ANY association?", Test 2 asks "given an
# association, does this score reach >= N metabolites?".
strata_all <- bind_rows(strata_t1, strata_t2)
write_csv(strata_all, file.path(out_dir, "sfari_score_stratified_enrichment.csv"))

# ---- Visualizations --------------------------------------------------------

# A0) Score-stratified OR forest, both tests. Log x-axis so OR and 1/OR are
# symmetric about the null. Zero-count strata give OR = 0 or Inf with an
# unbounded CI; those are dropped from the panel and reported in the CSV only.
forest_df <- strata_all |>
  mutate(stratum = factor(stratum,
                          levels = rev(c("Score 1", "Score 2", "Score 3",
                                         "Unscored", "All SFARI"))),
         test = factor(test, levels = c("Test 1", "Test 2")),
         sig  = ifelse(p < 0.05, "p < 0.05", "n.s."),
         lab  = sprintf("%d/%d", n_fg, n_tested)) |>
  filter(is.finite(OR), OR > 0, is.finite(ci_hi))

if (nrow(forest_df) < nrow(strata_all))
  cat(sprintf("\nNote: %d stratum/strata with unbounded OR omitted from the forest plot.\n",
              nrow(strata_all) - nrow(forest_df)))

p_forest <- ggplot(forest_df, aes(x = OR, y = stratum, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.6) +
  geom_point(size = 3) +
  geom_text(aes(label = lab), vjust = -1.1, size = 3, show.legend = FALSE) +
  facet_wrap(~ test, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(
               "Test 1" = "Test 1: any association (vs tested universe)",
               "Test 2" = sprintf("Test 2: >=%d metabolites (vs all hit genes)",
                                  MIN_METAB_RECURRENCE)))) +
  scale_x_log10() +
  scale_color_manual(values = c("p < 0.05" = "#E15759", "n.s." = "grey55"),
                     name = NULL) +
  labs(x = "Odds ratio vs non-SFARI genes (log scale)", y = NULL,
       title = "SFARI enrichment by gene score",
       caption = "Points labelled foreground/tested. p one-sided (greater); CI two-sided.") +
  theme_cowplot(12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92", color = NA))

save_dual_format(p_forest, out_dir, "sfari_score_forest", width = 8, height = 7)
print(p_forest)

# A) Recurrence distribution among hit genes: SFARI vs non-SFARI
p_dist <- ggplot(hit_recur |> mutate(group = ifelse(is_sfari, "SFARI", "Non-SFARI")),
                 aes(x = n_metabolites, fill = group)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.8) +
  geom_vline(xintercept = MIN_METAB_RECURRENCE - 0.5,
             linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("SFARI" = "#E15759", "Non-SFARI" = "#76B7B2"),
                    name = NULL) +
  annotate("text", x = MIN_METAB_RECURRENCE + 0.5, y = Inf, vjust = 1.5,
           hjust = 0, size = 3.5, color = "grey40",
           label = sprintf("%d-metabolite threshold", MIN_METAB_RECURRENCE)) +
  labs(x = "Distinct metabolites associated", y = "Number of genes",
       title = "Metabolite recurrence: SFARI vs non-SFARI genes (Test 2)") +
  theme_cowplot(12) +
  theme(legend.position = c(0.75, 0.85))

save_dual_format(p_dist, out_dir, "sfari_recurrence_distribution", width = 9, height = 6)
print(p_dist)

# B) Lollipop: genes associated with >= N metabolites, coloured by SFARI score
recurrent_df <- hit_recur |>
  filter(is_recurrent) |>
  mutate(sfari_label = case_when(
    !is_sfari              ~ "Non-SFARI",
    sfari_score == 1       ~ "SFARI (score 1)",
    sfari_score == 2       ~ "SFARI (score 2)",
    sfari_score == 3       ~ "SFARI (score 3)",
    TRUE                   ~ "SFARI (unscored)"
  ))

p_recurrent <- ggplot(recurrent_df,
                aes(x = reorder(gene, n_metabolites), y = n_metabolites,
                    color = sfari_label)) +
  geom_point(size = 3) +
  geom_segment(aes(xend = gene, yend = 0), linewidth = 0.5) +
  coord_flip() +
  scale_color_manual(values = c("Non-SFARI"         = "grey60",
                                 "SFARI (score 1)"   = "#E15759",
                                 "SFARI (score 2)"   = "#F28E2B",
                                 "SFARI (score 3)"   = "#EDC948",
                                 "SFARI (unscored)"  = "#BAB0AC"),
                     name = NULL) +
  labs(x = NULL, y = "Distinct metabolites associated",
       title = sprintf("Genes associated with \u2265%d metabolites \u2014 SFARI status",
                       MIN_METAB_RECURRENCE)) +
  theme_cowplot(12) +
  theme(legend.position = "right")

save_dual_format(p_recurrent, out_dir, "sfari_multimetabolite_genes_lollipop",
                 width = 10, height = 8)
print(p_recurrent)

# ---- TEST 3: SFARI enrichment per cell type --------------------------------
# Tests 1 and 2 pool cell types ("a hit in >=1 cell type"), which cannot tell a
# signal carried by one well-powered population from one present across the
# tissue. Test 3 repeats Test 1 within each cell type.
#
# BACKGROUND. Foreground = SFARI status of genes hit IN THAT CELL TYPE;
# background = genes TESTED in that cell type. 
cat("\n##########################################################\n")
cat("TEST 3: Is SFARI enrichment cell-type specific?\n")
cat("##########################################################\n")

ct_universe_cache <- file.path(cache_dir,
                               sprintf("tested_gene_universe_by_celltype_%s.rds", METHOD))

if (!REBUILD_UNIVERSE && file.exists(ct_universe_cache)) {
  ct_universe <- readRDS(ct_universe_cache)
} else {
  cat("Scanning parquet archive for per-cell-type tested-gene universe ...\n")
  ct_universe <- open_dataset(parquet_dir) |>
    dplyr::distinct(cell_type, gene) |> collect()
  saveRDS(ct_universe, ct_universe_cache)
}

hits_by_ct <- hits |> distinct(cell_type, gene)

ct_counts <- hits_by_ct |>
  count(cell_type, name = "n_hit_genes") |>
  mutate(tested = n_hit_genes >= MIN_HIT_GENES_PER_CT) |>
  arrange(desc(n_hit_genes))

cat(sprintf("Cell-type floor: %d distinct hit genes\n", MIN_HIT_GENES_PER_CT))
print(as.data.frame(ct_counts), right = FALSE)

# Silent, tabular version of run_fisher() -- Test 3 runs it once per cell type,
# so printing a contingency table each time would bury the result.
fisher_quiet <- function(fg, bg) {
  fg <- intersect(unique(fg), unique(bg))
  bg <- unique(bg)
  fg_s  <- sum(is_sfari_gene(fg))
  fg_ns <- length(fg) - fg_s
  bg_only <- setdiff(bg, fg)
  bg_s  <- sum(is_sfari_gene(bg_only))
  bg_ns <- length(bg_only) - bg_s

  if (fg_s == 0)
    return(tibble(n_hit = length(fg), n_bg = length(bg), n_sfari_hit = 0L,
                  pct_sfari_hit = 0, pct_sfari_bg = 100 * sum(is_sfari_gene(bg)) / length(bg),
                  OR = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_))

  m   <- matrix(c(fg_s, fg_ns, bg_s, bg_ns), nrow = 2, byrow = TRUE)
  ft  <- fisher.test(m, alternative = "greater")
  ft2 <- fisher.test(m)
  tibble(n_hit = length(fg), n_bg = length(bg), n_sfari_hit = fg_s,
         pct_sfari_hit = 100 * fg_s / length(fg),
         pct_sfari_bg  = 100 * sum(is_sfari_gene(bg)) / length(bg),
         OR = unname(ft$estimate),
         ci_lo = ft2$conf.int[1], ci_hi = ft2$conf.int[2], p = ft$p.value)
}

ct_keep <- ct_counts |> filter(tested) |> pull(cell_type)

sfari_by_ct <- map_dfr(ct_keep, function(ct) {
  fg <- hits_by_ct$gene[hits_by_ct$cell_type == ct]
  bg <- ct_universe$gene[ct_universe$cell_type == ct]
  fisher_quiet(fg, bg) |> mutate(cell_type = ct, .before = 1)
}) |>
  mutate(FDR = p.adjust(p, method = "BH"),
         significant = !is.na(FDR) & FDR < FDR_ALPHA) |>
  arrange(p)

print(as.data.frame(sfari_by_ct), right = FALSE, digits = 3)
cat(sprintf("%d of %d cell types enriched for SFARI genes at FDR < %.2f\n",
            sum(sfari_by_ct$significant), nrow(sfari_by_ct), FDR_ALPHA))

write_csv(sfari_by_ct, file.path(out_dir, "sfari_enrichment_by_celltype.csv"))
write_csv(ct_counts,   file.path(out_dir, "sfari_celltype_inclusion.csv"))

# ---- Additional visualizations ---------------------------------------------

# C) Per-cell-type forest. Ordered by OR; unbounded intervals (no SFARI hits)
# are dropped from the panel but retained in the CSV.
forest_ct <- sfari_by_ct |>
  filter(is.finite(OR), OR > 0, is.finite(ci_hi)) |>
  mutate(cell_type = fct_reorder(cell_type, OR),
         sig = ifelse(significant, sprintf("FDR < %.2f", FDR_ALPHA), "n.s."),
         lab = sprintf("%d/%d", n_sfari_hit, n_hit))

p_ct_forest <- ggplot(forest_ct, aes(x = OR, y = cell_type, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.6) +
  geom_point(size = 3) +
  geom_text(aes(label = lab), vjust = -1.1, size = 3, show.legend = FALSE) +
  scale_x_log10() +
  scale_color_manual(values = setNames(c("#E15759", "grey55"),
                                       c(sprintf("FDR < %.2f", FDR_ALPHA), "n.s.")),
                     name = NULL) +
  labs(x = "Odds ratio (log scale)", y = NULL,
       title = "SFARI enrichment by cell type",
       caption = "Points labelled SFARI hits / hit genes. Background = genes tested in that cell type.") +
  theme_cowplot(12) +
  theme(legend.position = "bottom")

save_dual_format(p_ct_forest, out_dir, "sfari_forest_by_celltype", width = 8, height = 7)
print(p_ct_forest)

# D) Score-response: proportion of TESTED SFARI genes at each score that reach
# at least one association. The forest gives odds ratios against a comparator;
# this gives the raw rate on an absolute scale, which is what shows whether the
# effect is graded across scores 1 -> 3 rather than driven by one stratum.
# Wilson intervals, which stay inside [0, 1] at small n where Wald does not.
wilson_ci <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2); ph <- k / n
  d <- 1 + z^2 / n
  ctr <- (ph + z^2 / (2 * n)) / d
  hw  <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / d
  c(max(0, ctr - hw), min(1, ctr + hw))
}

score_rate <- tibble(gene = background_genes) |>
  mutate(is_hit   = gene %in% hit_genes,
         is_sfari = is_sfari_gene(gene),
         score    = get_sfari_score(gene),
         stratum  = case_when(
           !is_sfari            ~ "Non-SFARI",
           is.na(score)         ~ "Unscored",
           TRUE                 ~ paste("Score", score))) |>
  group_by(stratum) |>
  summarise(n_tested = n(), n_hit = sum(is_hit), .groups = "drop") |>
  mutate(rate = n_hit / n_tested,
         ci   = map2(n_hit, n_tested, wilson_ci),
         lo   = map_dbl(ci, 1), hi = map_dbl(ci, 2),
         stratum = factor(stratum,
                          levels = c("Non-SFARI", "Unscored",
                                     "Score 3", "Score 2", "Score 1"))) |>
  dplyr::select(-ci) |>
  arrange(stratum)

print(as.data.frame(score_rate), right = FALSE, digits = 3)
write_csv(score_rate, file.path(out_dir, "sfari_hit_rate_by_score.csv"))

p_score_rate <- ggplot(score_rate, aes(x = stratum, y = 100 * rate,
                                       fill = stratum)) +
  geom_col(width = 0.68, alpha = 0.9) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.16,
                linewidth = 0.5) +
  geom_text(aes(label = sprintf("%d/%d", n_hit, n_tested)),
            vjust = -0.6, hjust = -0.15, size = 3) +
  geom_hline(yintercept = 100 * with(score_rate,
                                     sum(n_hit) / sum(n_tested)),
             linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Score 1" = "#E15759", "Score 2" = "#F28E2B",
                               "Score 3" = "#EDC948", "Unscored" = "#BAB0AC",
                               "Non-SFARI" = "#76B7B2"), guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Tested genes with >=1 metabolite association (%)",
       title = "Association rate across SFARI confidence scores",
       subtitle = "Wilson 95% intervals; dashed line = overall rate across the tested universe") +
  theme_cowplot(12)

save_dual_format(p_score_rate, out_dir, "sfari_hit_rate_by_score",
                 width = 8, height = 5)
print(p_score_rate)

# E) Score x cell type. Sparse by construction -- score 1 has few genes and a
# single cell type contributes few hits -- so this is plotted as raw counts
# with a rate fill rather than as a grid of odds ratios, which would be mostly
# unstable estimates. Guarded: skipped entirely if too thin to read.
score_ct <- ct_universe |>
  filter(cell_type %in% ct_keep) |>
  mutate(is_sfari = is_sfari_gene(gene),
         score    = get_sfari_score(gene)) |>
  filter(is_sfari) |>
  mutate(stratum = ifelse(is.na(score), "Unscored", paste("Score", score))) |>
  left_join(hits_by_ct |> mutate(is_hit = TRUE), by = c("cell_type", "gene")) |>
  mutate(is_hit = replace_na(is_hit, FALSE)) |>
  group_by(cell_type, stratum) |>
  summarise(n_tested = n(), n_hit = sum(is_hit),
            rate = mean(is_hit), .groups = "drop")

write_csv(score_ct, file.path(out_dir, "sfari_score_by_celltype.csv"))

if (nrow(score_ct) >= 6 && sum(score_ct$n_hit) >= 10) {
  p_score_ct <- ggplot(score_ct, aes(x = cell_type, y = stratum, fill = 100 * rate)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%d/%d", n_hit, n_tested)), size = 2.8) +
    scale_fill_gradient(low = "white", high = "#E15759",
                        name = "Hit rate (%)") +
    labs(x = NULL, y = NULL,
         title = "SFARI genes reaching an association, by score and cell type",
         subtitle = "Cell labels = hit genes / tested SFARI genes in that cell type") +
    theme_cowplot(11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_dual_format(p_score_ct, out_dir, "sfari_score_by_celltype",
                   width = 10, height = 5)
  print(p_score_ct)
} else {
  cat("\nScore x cell type panel skipped: too few SFARI hits to plot meaningfully.\n")
}

# ---- Threshold-free enrichment (no FDR cut anywhere) -----------------------
# Every test above depends on the FDR < 0.10 / |logFC| >= 0.25 cut. A reviewer
# can reasonably ask whether SFARI enrichment is a property of the data or of
# where that line was drawn -- which is exactly what the co-author's
# "what about FDR < 0.05?" question was getting at.
#
# So: rank EVERY tested gene by its strongest association (smallest p-value
# across all metabolite x cell-type tests) and walk down the ranking, adding
# when a gene is SFARI and subtracting when it is not. A curve that climbs
# smoothly to an early peak means SFARI genes sit toward the top of the whole
# ranking, independent of any threshold. The statistic is the classic
# unweighted Kolmogorov-Smirnov running sum used by GSEA; significance comes
# from permuting gene labels, which preserves the ranking and the set size.
rank_cache <- file.path(cache_dir, sprintf("gene_min_pvalue_%s.rds", METHOD))

if (!REBUILD_UNIVERSE && file.exists(rank_cache)) {
  gene_rank <- readRDS(rank_cache)
} else {
  cat("\nScanning parquet archive for per-gene minimum p-value ...\n")
  gene_rank <- open_dataset(parquet_dir) |>
    group_by(gene) |>
    summarise(min_p = min(P.Value, na.rm = TRUE)) |>
    collect()
  saveRDS(gene_rank, rank_cache)
}

gene_rank <- gene_rank |>
  filter(is.finite(min_p)) |>
  arrange(min_p) |>
  mutate(is_sfari = is_sfari_gene(gene))

cat(sprintf("Ranked genes: %s (%d SFARI)\n",
            format(nrow(gene_rank), big.mark = ","), sum(gene_rank$is_sfari)))

# Running enrichment score over a 0/1 membership vector.
running_es <- function(hit_vec) {
  n <- length(hit_vec); nh <- sum(hit_vec)
  if (nh == 0 || nh == n) return(rep(0, n))
  cumsum(ifelse(hit_vec, 1 / nh, -1 / (n - nh)))
}

es_curve <- running_es(gene_rank$is_sfari)
es_obs   <- es_curve[which.max(abs(es_curve))]
peak_i   <- which.max(abs(es_curve))

set.seed(123)
N_PERM_ES <- 10000
es_null <- replicate(N_PERM_ES, {
  v <- sample(gene_rank$is_sfari)
  e <- running_es(v)
  e[which.max(abs(e))]
})
p_es <- (sum(es_null >= es_obs) + 1) / (N_PERM_ES + 1)

cat(sprintf("Threshold-free enrichment: ES = %+.4f at rank %s of %s, permutation p = %.4g\n",
            es_obs, format(peak_i, big.mark = ","),
            format(nrow(gene_rank), big.mark = ","), p_es))

# F) Running enrichment curve with the SFARI rug beneath it.
curve_df <- tibble(rank = seq_along(es_curve), es = es_curve)
rug_df   <- tibble(rank = which(gene_rank$is_sfari))

p_curve <- plot_grid(
  ggplot(curve_df, aes(x = rank, y = es)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_line(colour = "#E15759", linewidth = 0.8) +
    geom_vline(xintercept = peak_i, linetype = "dotted", colour = "grey40") +
    annotate("text", x = peak_i, y = es_obs, hjust = -0.1, vjust = -0.5,
             size = 3.4,
             label = sprintf("ES = %+.3f\np = %.3g", es_obs, p_es)) +
    labs(x = NULL, y = "Running enrichment",
         title = "SFARI enrichment without a significance threshold",
         subtitle = "All tested genes ranked by strongest metabolite association") +
    theme_cowplot(12) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  ggplot(rug_df, aes(x = rank)) +
    geom_segment(aes(xend = rank, y = 0, yend = 1), linewidth = 0.2,
                 colour = "grey25") +
    scale_x_continuous(limits = c(1, nrow(gene_rank)),
                       labels = scales::comma) +
    labs(x = "Gene rank (1 = strongest association)", y = NULL) +
    theme_cowplot(12) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
  ncol = 1, align = "v", rel_heights = c(4, 1))

save_dual_format(p_curve, out_dir, "sfari_enrichment_curve", width = 9, height = 6)
print(p_curve)

write_csv(tibble(es = es_obs, peak_rank = peak_i, n_ranked = nrow(gene_rank),
                 n_sfari = sum(gene_rank$is_sfari), n_perm = N_PERM_ES, p = p_es),
          file.path(out_dir, "sfari_enrichment_curve_stats.csv"))

# ---- Covariate-matched permutation null ------------------------------------
# The confound both this script and constraint_enrichment.R document: hit genes
# are longer and better expressed than background, so ANY gene set biased
# toward large, well-expressed genes will appear enriched. Fisher's test cannot
# see that. This can.
#
# Draw N random gene sets the same size as the tested SFARI set, sampling
# WITHIN strata defined by expression decile x number of cell types in which
# the gene was testable -- the two things that most directly determine whether
# a gene could become a hit at all. If the observed overlap still sits in the
# tail of that null, the enrichment is not a detectability artefact.
N_PERM_MATCH <- 10000

cpm       <- pb_base$counts / rowSums(pb_base$counts) * 1e6
mean_expr <- log2(colMeans(cpm) + 1)

n_ct_tested <- ct_universe |> count(gene, name = "n_ct_tested")

match_df <- tibble(gene = background_genes) |>
  left_join(n_ct_tested, by = "gene") |>
  mutate(mean_expr   = unname(mean_expr[gene]),
         n_ct_tested = replace_na(n_ct_tested, 0L),
         is_hit      = gene %in% hit_genes,
         is_sfari    = is_sfari_gene(gene)) |>
  filter(!is.na(mean_expr)) |>
  mutate(expr_decile = dplyr::ntile(mean_expr, 10),
         stratum     = paste(expr_decile, n_ct_tested, sep = "_"))

obs_overlap <- sum(match_df$is_sfari & match_df$is_hit)

# Sample the same number of genes from each stratum as SFARI occupies there.
strata_n <- match_df |> filter(is_sfari) |> count(stratum, name = "k")
pool     <- split(match_df$is_hit, match_df$stratum)

set.seed(123)
null_overlap <- replicate(N_PERM_MATCH, {
  sum(map2_int(strata_n$stratum, strata_n$k, function(s, k) {
    v <- pool[[s]]
    if (is.null(v) || length(v) == 0) return(0L)
    sum(sample(v, min(k, length(v))))
  }))
})

p_match <- (sum(null_overlap >= obs_overlap) + 1) / (N_PERM_MATCH + 1)
fold    <- obs_overlap / mean(null_overlap)

cat(sprintf("\nMatched permutation: observed %d SFARI hit genes vs null mean %.1f (%.2fx), p = %.4g\n",
            obs_overlap, mean(null_overlap), fold, p_match))
cat(sprintf("  Strata: %d (expression decile x cell types tested)\n", nrow(strata_n)))

# G) Null distribution with the observed value marked.
p_null <- ggplot(tibble(x = null_overlap), aes(x = x)) +
  geom_histogram(bins = 40, fill = "#BAB0AC", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs_overlap, colour = "#E15759", linewidth = 1) +
  annotate("text", x = obs_overlap, y = Inf, hjust = -0.08, vjust = 1.6,
           size = 3.6, colour = "#E15759",
           label = sprintf("observed = %d\n%.2fx, p = %.3g",
                           obs_overlap, fold, p_match)) +
  labs(x = "SFARI genes among hits in matched random sets",
       y = "Permutations",
       title = "Enrichment against a detectability-matched null",
       subtitle = sprintf("%s permutations, sampled within expression-decile x cell-types-tested strata",
                          format(N_PERM_MATCH, big.mark = ","))) +
  theme_cowplot(12)

save_dual_format(p_null, out_dir, "sfari_matched_permutation_null",
                 width = 8, height = 5)
print(p_null)

write_csv(tibble(observed = obs_overlap, null_mean = mean(null_overlap),
                 null_sd = sd(null_overlap), fold = fold,
                 n_strata = nrow(strata_n), n_perm = N_PERM_MATCH, p = p_match),
          file.path(out_dir, "sfari_matched_permutation.csv"))

# ---- H) Which SFARI genes, in which cell types? ----------------------------
# The odds ratio says enrichment exists; this says what it is made of. Fill is
# the number of distinct metabolites, so a reader can see whether one cell type
# carries the result and which genes drive it.
sfari_ct <- hits |>
  filter(is_sfari_gene(gene)) |>
  group_by(gene, cell_type) |>
  summarise(n_metabolites = n_distinct(Compound.ID), .groups = "drop")

if (nrow(sfari_ct) > 0) {
  gene_order <- sfari_ct |>
    group_by(gene) |>
    summarise(total = sum(n_metabolites), .groups = "drop") |>
    arrange(total) |>
    pull(gene)

  p_sfari_ct <- sfari_ct |>
    mutate(gene = factor(gene, levels = gene_order),
           score = get_sfari_score(as.character(gene))) |>
    ggplot(aes(x = cell_type, y = gene, fill = n_metabolites)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = n_metabolites), size = 2.6, colour = "grey15") +
    scale_fill_gradient(low = "#FDE0DD", high = "#B2182B",
                        name = "Metabolites") +
    labs(x = NULL, y = NULL,
         title = "SFARI genes with metabolite associations",
         subtitle = "Genes ordered by total distinct metabolites across cell types") +
    theme_cowplot(11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7))

  save_dual_format(p_sfari_ct, out_dir, "sfari_gene_by_celltype_heatmap",
                   width = 9, height = max(6, 0.16 * n_distinct(sfari_ct$gene)))
  print(p_sfari_ct)
  write_csv(sfari_ct, file.path(out_dir, "sfari_gene_by_celltype.csv"))
}

# ---- I) Metabolite class composition ---------------------------------------
# A question none of the odds-ratio panels ask: do SFARI genes associate with a
# DIFFERENT KIND of metabolite, not merely more of them? A flat comparison is a
# clean negative worth one sentence; a skew toward one class is a result.
# Guarded on the annotation column, which may not be carried in the hits CSV.
class_col <- intersect(c("Class", "Super.Class", "Subclass"), names(hits))

if (length(class_col) > 0) {
  cc <- class_col[1]

  class_df <- hits |>
    mutate(group = ifelse(is_sfari_gene(gene), "SFARI", "Non-SFARI"),
           metab_class = replace_na(as.character(.data[[cc]]), "Unannotated")) |>
    distinct(gene, Compound.ID, group, metab_class) |>
    count(group, metab_class, name = "n") |>
    group_by(group) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup()

  # Keep the classes that matter to either group; collapse the rest so the
  # panel is readable rather than a 40-row strip of near-zero bars.
  keep <- class_df |> group_by(metab_class) |>
    summarise(mx = max(pct), .groups = "drop") |>
    filter(mx >= 2) |> pull(metab_class)

  class_plot_df <- class_df |>
    mutate(metab_class = ifelse(metab_class %in% keep, metab_class,
                                "Other (<2%)")) |>
    group_by(group, metab_class) |>
    summarise(n = sum(n), pct = sum(pct), .groups = "drop")

  # Chi-square on the un-collapsed counts, if the table is large enough.
  tab <- class_df |>
    dplyr::select(group, metab_class, n) |>
    pivot_wider(names_from = group, values_from = n, values_fill = 0) |>
    column_to_rownames("metab_class") |>
    as.matrix()

  chisq_res <- if (all(dim(tab) >= 2) && sum(tab) > 0) {
    suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000))
  } else NULL

  if (!is.null(chisq_res)) {
    cat(sprintf("\nMetabolite class composition, SFARI vs non-SFARI: chi-square p = %.4g (simulated)\n",
                chisq_res$p.value))

    # Overall test statistics. The p-value is Monte Carlo (B = 10,000) because
    # most metabolite classes have small expected counts, where the asymptotic
    # chi-square approximation is unreliable; the asymptotic version is carried
    # alongside so the two can be compared. `min_expected` is the diagnostic
    # that justifies the simulation -- the usual rule of thumb is that the
    # asymptotic test needs expected counts of at least 5.
    chisq_asym <- suppressWarnings(chisq.test(tab))

    chisq_stats <- tibble(
      comparison      = "SFARI vs non-SFARI hit genes",
      unit            = "distinct gene-metabolite pairs",
      class_column    = cc,
      n_classes       = nrow(tab),
      n_pairs_total   = sum(tab),
      n_pairs_sfari   = sum(tab[, intersect("SFARI", colnames(tab)), drop = FALSE]),
      n_pairs_nonsfari = sum(tab[, intersect("Non-SFARI", colnames(tab)), drop = FALSE]),
      statistic       = unname(chisq_res$statistic),
      df_asymptotic   = unname(chisq_asym$parameter),
      p_simulated     = chisq_res$p.value,
      p_asymptotic    = chisq_asym$p.value,
      B               = 10000,
      min_expected    = min(chisq_asym$expected),
      n_cells_exp_lt5 = sum(chisq_asym$expected < 5),
      method          = chisq_res$method)

    print(as.data.frame(chisq_stats), right = FALSE, digits = 4)
    write_csv(chisq_stats,
              file.path(out_dir, "sfari_metabolite_class_chisq.csv"))

    # Per-class contributions. The overall p-value says only that composition
    # differs somewhere; standardised residuals say WHICH classes drive it and
    # in which direction (positive = over-represented in that group). These are
    # the numbers to quote if the test is significant.
    resid_df <- as.data.frame.table(chisq_asym$stdres,
                                    responseName = "std_residual") |>
      rename(metab_class = Var1, group = Var2) |>
      left_join(as.data.frame.table(chisq_asym$expected,
                                    responseName = "expected") |>
                  rename(metab_class = Var1, group = Var2),
                by = c("metab_class", "group")) |>
      left_join(as.data.frame.table(tab, responseName = "observed") |>
                  rename(metab_class = Var1, group = Var2),
                by = c("metab_class", "group")) |>
      mutate(across(c(metab_class, group), as.character),
             # Two-sided normal approximation on the standardised residual,
             # BH-corrected across classes. Descriptive follow-up to a
             # significant omnibus test, not an independent set of tests.
             p_residual   = 2 * pnorm(-abs(std_residual)),
             FDR_residual = p.adjust(p_residual, method = "BH")) |>
      arrange(desc(abs(std_residual)))

    write_csv(resid_df,
              file.path(out_dir, "sfari_metabolite_class_residuals.csv"))

    cat("\nLargest class-level deviations (standardised residuals):\n")
    print(as.data.frame(head(resid_df, 10)), right = FALSE, digits = 3)
  }

  p_class <- ggplot(class_plot_df,
                    aes(x = pct, y = fct_reorder(metab_class, pct),
                        fill = group)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             alpha = 0.9) +
    scale_fill_manual(values = c("SFARI" = "#E15759",
                                 "Non-SFARI" = "#76B7B2"), name = NULL) +
    labs(x = "Percent of gene-metabolite pairs", y = NULL,
         title = "Metabolite classes associated with SFARI vs non-SFARI genes",
         subtitle = if (!is.null(chisq_res))
           sprintf("Chi-square p = %.3g (simulated); classes below 2%% collapsed",
                   chisq_res$p.value) else "Classes below 2% collapsed") +
    theme_cowplot(11) +
    theme(legend.position = "bottom")

  save_dual_format(p_class, out_dir, "sfari_metabolite_class_composition",
                   width = 9, height = 7)
  print(p_class)
  write_csv(class_df, file.path(out_dir, "sfari_metabolite_class_composition.csv"))
} else {
  cat("\nMetabolite class panel skipped: no Class/Super.Class column in the hits CSV.\n")
}

# ---- Summary ---------------------------------------------------------------
cat("\n========== SUMMARY ==========\n")
cat(sprintf("Method: %s | background = %d tested genes (union across cell types)\n",
            METHOD, length(background_genes)))
cat(sprintf("Test 1 (SFARI in hits vs tested universe): Fisher OR=%.2f, p=%.2e\n",
            fisher_t1$estimate, fisher_t1$p.value))
cat(sprintf("Test 2 (SFARI among >=%d-metabolite genes): Fisher OR=%.2f, p=%.2e; Wilcoxon p=%.2e\n",
            MIN_METAB_RECURRENCE,
            fisher_t2$estimate, fisher_t2$p.value, wilcox_t2$p.value))
cat(sprintf("Test 3 (per cell type): %d/%d cell types enriched at FDR < %.2f (floor = %d hit genes)\n",
            sum(sfari_by_ct$significant), nrow(sfari_by_ct), FDR_ALPHA,
            MIN_HIT_GENES_PER_CT))
cat(sprintf("Threshold-free (no FDR cut):  ES = %+.4f, permutation p = %.4g\n",
            es_obs, p_es))
cat(sprintf("Detectability-matched null:   %d observed vs %.1f expected (%.2fx), p = %.4g\n",
            obs_overlap, mean(null_overlap), fold, p_match))

# ---- Save results ----------------------------------------------------------
sfari_overlap <- hit_recur |>
  filter(is_sfari) |>
  left_join(sfari |> dplyr::select(symbol, gene_name, score, syndromic,
                                    `genetic-category`, `number-of-reports`),
            by = c("gene" = "symbol")) |>
  arrange(desc(n_metabolites))

write_csv(sfari_overlap, file.path(out_dir, "sfari_gene_metabolite_overlap.csv"))
write_csv(hit_recur,     file.path(out_dir, "gene_recurrence_with_sfari.csv"))

cat(sprintf("\nResults saved to %s\n", out_dir))
cat(sprintf("SFARI genes in hits: %d\n", nrow(sfari_overlap)))
cat(sprintf("SFARI genes associated with >=%d metabolites: %d\n",
            MIN_METAB_RECURRENCE,
            sum(sfari_overlap$n_metabolites >= MIN_METAB_RECURRENCE)))
