#!/usr/bin/env Rscript
# sfari_enrichment.R
# ---------------------------------------------------------------------------
# Two clean enrichment tests for SFARI autism risk genes:
#
#   Test 1: Are SFARI genes enriched among genes with ANY significant
#           gene-metabolite association?
#           Foreground = hit genes, Background = TESTED gene universe
#
#   Test 2: Among hit genes, are SFARI genes more likely to be hubs
#           (>= N distinct metabolite associations)?
#           Foreground = hub genes, Background = hit genes
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
         is_hub      = n_metabolites >= MIN_METAB_RECURRENCE) |>
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

# ---- TEST 2: Among hit genes, are SFARI genes more likely hubs? ------------
cat("\n##########################################################\n")
cat("TEST 2: Among hit genes, are SFARI genes more likely hubs?\n")
cat("  Foreground: hub genes >= ", MIN_METAB_RECURRENCE, " metabolites (",
    sum(hit_recur$is_hub), ")\n")
cat("  Background: all hit genes (", nrow(hit_recur), ")\n")
cat("##########################################################\n")

hub_genes <- hit_recur |> filter(is_hub) |> pull(gene)
fisher_t2 <- run_fisher(hub_genes, hit_genes, "Hub", "Hit")

# Wilcoxon: among hit genes, do SFARI genes have higher recurrence?
wilcox_t2 <- wilcox.test(n_metabolites ~ is_sfari, data = hit_recur,
                         alternative = "greater")
cat(sprintf("\nWilcoxon (SFARI higher recurrence among hits): W = %.0f, p = %.4f\n",
            wilcox_t2$statistic, wilcox_t2$p.value))

hit_recur |>
  group_by(is_sfari) |>
  summarise(n = n(), median_metab = median(n_metabolites),
            mean_metab = mean(n_metabolites),
            n_hub = sum(is_hub), .groups = "drop") |>
  print()

# ---- SFARI score stratification (Test 2) -----------------------------------
cat("\n--- SFARI score stratification (among hit genes) ---\n")
sfari_hits <- hit_recur |> filter(is_sfari, !is.na(sfari_score))

sfari_hits |>
  group_by(sfari_score) |>
  summarise(n = n(), n_hub = sum(is_hub),
            pct_hub = 100 * mean(is_hub),
            median_metab = median(n_metabolites), .groups = "drop") |>
  print()

strata_t2 <- run_score_strata(hub_genes, hit_genes, "Test 2", "hub")

# Both tests share the helper, so the score axis is defined identically:
# Test 1 asks "does this score reach ANY association?", Test 2 asks "given an
# association, does this score reach hub status?".
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
               "Test 2" = "Test 2: hub status (vs all hit genes)"))) +
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
           label = sprintf("Hub threshold (%d)", MIN_METAB_RECURRENCE)) +
  labs(x = "Distinct metabolites associated", y = "Number of genes",
       title = "Metabolite recurrence: SFARI vs non-SFARI genes (Test 2)") +
  theme_cowplot(12) +
  theme(legend.position = c(0.75, 0.85))

save_dual_format(p_dist, out_dir, "sfari_recurrence_distribution", width = 9, height = 6)
print(p_dist)

# B) Lollipop: hub genes colored by SFARI status/score
hub_df <- hit_recur |>
  filter(is_hub) |>
  mutate(sfari_label = case_when(
    !is_sfari              ~ "Non-SFARI",
    sfari_score == 1       ~ "SFARI (score 1)",
    sfari_score == 2       ~ "SFARI (score 2)",
    sfari_score == 3       ~ "SFARI (score 3)",
    TRUE                   ~ "SFARI (unscored)"
  ))

p_hub <- ggplot(hub_df,
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
       title = sprintf("Hub genes (\u2265%d metabolites) \u2014 SFARI status",
                       MIN_METAB_RECURRENCE)) +
  theme_cowplot(12) +
  theme(legend.position = "right")

save_dual_format(p_hub, out_dir, "sfari_hub_genes_lollipop", width = 10, height = 8)
print(p_hub)

# ---- Summary ---------------------------------------------------------------
cat("\n========== SUMMARY ==========\n")
cat(sprintf("Method: %s | background = %d tested genes (union across cell types)\n",
            METHOD, length(background_genes)))
cat(sprintf("Test 1 (SFARI in hits vs tested universe): Fisher OR=%.2f, p=%.2e\n",
            fisher_t1$estimate, fisher_t1$p.value))
cat(sprintf("Test 2 (SFARI hubs vs hits):          Fisher OR=%.2f, p=%.2e; Wilcoxon p=%.2e\n",
            fisher_t2$estimate, fisher_t2$p.value, wilcox_t2$p.value))

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
cat(sprintf("SFARI hub genes: %d\n", sum(sfari_overlap$n_metabolites >= MIN_METAB_RECURRENCE)))
