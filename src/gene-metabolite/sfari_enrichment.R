#!/usr/bin/env Rscript
# sfari_enrichment.R
# ---------------------------------------------------------------------------
# Two clean enrichment tests for SFARI autism risk genes:
#
#   Test 1: Are SFARI genes enriched among genes with ANY significant
#           gene-metabolite association?
#           Foreground = hit genes (542), Background = all pseudobulk genes (~20k)
#
#   Test 2: Among hit genes, are SFARI genes more likely to be hubs
#           (>= N distinct metabolite associations)?
#           Foreground = hub genes, Background = hit genes (542)
#
# Each test: Fisher's exact + Wilcoxon rank-sum + SFARI score stratification.
#
# Run from src/gene-metabolite/:
#   Rscript sfari_enrichment.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(svglite)
})

source("pseudobulk_functions.R")

# ---- Config ----------------------------------------------------------------
hits_path  <- "../../results/gene-metabolite/csv-log2_na/metabolite_gene_hits.csv"
sfari_path <- "../../doc/gene_lists/SFARI-Gene_genes_07-12-2026release_08-13-2026export.csv"
cache_dir  <- "../../data/cache"
out_dir    <- "../../results/gene-metabolite/sfari-enrichment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_METAB_RECURRENCE <- 5

# ---- Load data -------------------------------------------------------------
hits <- read_csv(hits_path, show_col_types = FALSE)

pb_base <- readRDS(file.path(cache_dir, "pb_base.rds"))
background_genes <- colnames(pb_base$counts)
cat(sprintf("Background genes (from pseudobulk): %d\n", length(background_genes)))
cat(sprintf("Gene-metabolite hits: %d associations, %d unique genes\n",
            nrow(hits), n_distinct(hits$gene)))

sfari <- read_csv(sfari_path, show_col_types = FALSE) |>
  rename(symbol = `gene-symbol`, ensembl = `ensembl-id`,
         score = `gene-score`, gene_name = `gene-name`)

sfari_symbols <- sfari$symbol
sfari_ensembl <- sfari$ensembl

# ---- SFARI lookup -----------------------------------------------------------
sfari_by_sym <- sfari |> dplyr::select(symbol, score) |> deframe()
sfari_by_ens <- sfari |> dplyr::select(ensembl, score) |> deframe()

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

# ---- TEST 1: SFARI enrichment among hit genes ------------------------------
cat("\n##########################################################\n")
cat("TEST 1: Are SFARI genes enriched among significant genes?\n")
cat("  Foreground: genes with any hit (", length(hit_genes), ")\n")
cat("  Background: all pseudobulk genes (", length(background_genes), ")\n")
cat("##########################################################\n")

fisher_t1 <- run_fisher(hit_genes, background_genes, "Hit", "Background")

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

non_sfari_hits <- hit_recur |> filter(!is_sfari)
ns_hub    <- sum(non_sfari_hits$is_hub)
ns_nonhub <- sum(!non_sfari_hits$is_hub)

for (s in sort(unique(sfari_hits$sfari_score))) {
  s_hub    <- sum(sfari_hits$sfari_score == s & sfari_hits$is_hub)
  s_nonhub <- sum(sfari_hits$sfari_score == s & !sfari_hits$is_hub)
  mat <- matrix(c(s_hub, ns_hub, s_nonhub, ns_nonhub), nrow = 2)
  ft  <- fisher.test(mat, alternative = "greater")
  cat(sprintf("  Score %s: %d/%d hub (%.0f%%), OR=%.2f, p=%.2e\n",
              s, s_hub, s_hub + s_nonhub,
              100 * s_hub / max(s_hub + s_nonhub, 1),
              ft$estimate, ft$p.value))
}

# ---- Visualizations --------------------------------------------------------

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
cat(sprintf("Test 1 (SFARI in hits vs background): Fisher OR=%.2f, p=%.2e\n",
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
