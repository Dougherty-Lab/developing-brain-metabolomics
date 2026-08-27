#!/usr/bin/env Rscript
# constraint_enrichment.R
# ---------------------------------------------------------------------------
# Is pLoF constraint (gnomAD LOEUF) elevated among genes with significant
# gene-metabolite associations? Mirrors sfari_enrichment.R.
#
#   Test 1: hit genes vs the TESTED gene universe
#   Test 2: hub genes (>= N metabolites) vs all hit genes
#   Test 3: LOEUF vs metabolite recurrence among hit genes (Spearman)
#
# LOEUF is continuous and gnomAD explicitly recommends using it that way, so
# the PRIMARY test for 1 and 2 is a Wilcoxon rank-sum on the LOEUF
# distribution. The percentile ladder (0.15 / 0.27 / 0.36 / 0.45 / 0.60) is a
# SECONDARY sensitivity analysis only -- reported as a Fisher OR at every
# threshold rather than at one hand-picked cut, so the reader can see whether
# the answer depends on where the line is drawn.
#
# THE CONFOUND. LOEUF is not independent of statistical power. It is bounded
# by the expected pLoF count, so short genes get wide CIs and inflated LOEUF
# regardless of true constraint (Karczewski et al. 2020 Nature). Constrained
# genes are also longer and more highly expressed, which raises power to
# detect a gene-metabolite association in limma-voom. "Hits are more
# constrained" can therefore be manufactured entirely by gene size and
# expression. Three defences, all reported:
#   (a) background = tested genes only  (see sfari_enrichment.R rationale)
#   (b) logistic regression: hit ~ LOEUF + mean expression + log10(expected LoF)
#   (c) expression-decile-matched resampling of the non-hit background
# If (a)-(c) agree with the raw Wilcoxon, report the Wilcoxon and cite the
# rest as robustness. If they disagree, the raw result was power, not biology.
#
# Input: gnomAD v4.1 constraint metrics TSV (MANE Select), or v2.1.1 -- the
# LOEUF and expected-LoF columns are auto-detected. Download from
# https://gnomad.broadinstitute.org/downloads#v4-constraint
#
# Run from src/gene-metabolite/:
#   Rscript constraint_enrichment.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(cowplot)
  library(svglite)
})

source("pseudobulk_functions.R")

set.seed(123)

# ---- Config ----------------------------------------------------------------
METHOD      <- "log2_na"                       # "int" | "zscore" | "zscore_trim" | "log2_na"
suffix      <- if (METHOD == "int") "" else paste0("-", METHOD)
parquet_dir <- if (METHOD == "int") {
  "../../results/gene-metabolite/parquet"
} else {
  sprintf("../../results/gene-metabolite/parquet-%s", METHOD)
}

hits_path      <- sprintf("../../results/gene-metabolite/csv%s/metabolite_gene_hits.csv", suffix)
constraint_path <- "../../doc/gene_lists/gnomad.v4.1.1.constraint_metrics.tsv"
sfari_path     <- "../../doc/gene_lists/SFARI-Gene_genes_07-12-2026release_08-13-2026export.csv"
cache_dir      <- "../../data/cache"
out_dir        <- "../../results/gene-metabolite/constraint-enrichment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_METAB_RECURRENCE <- 5      # hub definition, matches sfari_enrichment.R
MIN_EXP_LOF          <- 10     # LOEUF unreliable below ~10 expected pLoF variants
N_MATCH_REPS         <- 1000   # expression-matched resampling iterations
REBUILD_UNIVERSE     <- FALSE

# gnomAD's published LOEUF -> percentile ladder (17,063 MANE Select transcripts)
LOEUF_LADDER <- tibble(
  threshold  = c(0.15, 0.27, 0.36, 0.45, 0.60),
  percentile = c(1,    5,    10,   15,   25)
)
PRIMARY_THRESHOLD <- 0.45      # gnomAD's suggested Mendelian cut, for highlighting

# ---- Load hits + tested universe -------------------------------------------
hits <- read_csv(hits_path, show_col_types = FALSE)

universe_cache <- file.path(cache_dir, sprintf("tested_gene_universe_%s.rds", METHOD))
if (!REBUILD_UNIVERSE && file.exists(universe_cache)) {
  background_genes <- readRDS(universe_cache)
  cat(sprintf("Tested-gene universe (cached): %d\n", length(background_genes)))
} else {
  cat("Scanning parquet archive for tested-gene universe ...\n")
  background_genes <- open_dataset(parquet_dir) |>
    dplyr::distinct(gene) |> collect() |> dplyr::pull(gene) |> unique() |> sort()
  saveRDS(background_genes, universe_cache)
  cat(sprintf("Tested-gene universe (scanned): %d\n", length(background_genes)))
}

hit_recur <- hits |>
  group_by(gene) |>
  summarise(n_metabolites  = n_distinct(Compound.ID),
            n_cell_types   = n_distinct(cell_type),
            n_associations = n(), .groups = "drop")

cat(sprintf("Hits: %d associations, %d unique genes\n",
            nrow(hits), nrow(hit_recur)))

# ---- Load gnomAD constraint ------------------------------------------------
# gnomAD ships the constraint TSV bgzipped (.bgz). BGZF is a valid multi-member
# gzip stream, so gzfile() reads it directly -- no need to decompress on disk.
# readr only auto-detects .gz/.bz2/.xz/.zip by extension, so a bare .bgz path
# would fail with an opaque parse error; resolve and wrap it explicitly.
resolve_constraint_path <- function(p) {
  cands <- c(p, paste0(p, ".bgz"), paste0(p, ".gz"))
  hit   <- cands[file.exists(cands)]
  if (length(hit) == 0)
    stop("gnomAD constraint file not found. Looked for:\n  ",
         paste(cands, collapse = "\n  "),
         "\nDownload the constraint metrics TSV from ",
         "https://gnomad.broadinstitute.org/downloads#v4-constraint")
  hit[1]
}

constraint_path <- resolve_constraint_path(constraint_path)
cat(sprintf("Constraint file: %s\n", basename(constraint_path)))

constraint_raw <- if (grepl("\\.(bgz|gz)$", constraint_path)) {
  read_tsv(gzfile(constraint_path), show_col_types = FALSE, guess_max = 100000)
} else {
  read_tsv(constraint_path, show_col_types = FALSE, guess_max = 100000)
}

# Column names differ between releases; detect rather than hard-code.
#   v4.1 : lof.oe_ci.upper / lof.exp / gene / gene_id / mane_select
#   v2.1 : oe_lof_upper    / exp_lof / gene / gene_id / canonical
pick_col <- function(df, candidates, what) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0)
    stop(sprintf("No %s column found. Looked for: %s\nAvailable: %s", what,
                 paste(candidates, collapse = ", "),
                 paste(head(names(df), 40), collapse = ", ")))
  hit[1]
}

col_loeuf  <- pick_col(constraint_raw, c("lof.oe_ci.upper", "oe_lof_upper"), "LOEUF")
col_exp    <- pick_col(constraint_raw, c("lof.exp", "exp_lof"), "expected-pLoF")
col_symbol <- pick_col(constraint_raw, c("gene", "gene_symbol"), "gene-symbol")
col_ensg   <- pick_col(constraint_raw, c("gene_id", "ensembl_gene_id"), "gene-ID")
cat(sprintf("gnomAD columns: LOEUF=%s | expected pLoF=%s\n", col_loeuf, col_exp))

# One row per gene. Prefer MANE Select / canonical; if neither flag exists,
# keep the transcript with the most expected pLoF variants (best-powered).
flag_col <- intersect(c("mane_select", "canonical"), names(constraint_raw))
constraint <- constraint_raw |>
  dplyr::transmute(symbol = .data[[col_symbol]],
                   ensembl = sub("\\..*$", "", .data[[col_ensg]]),
                   loeuf  = as.numeric(.data[[col_loeuf]]),
                   exp_lof = as.numeric(.data[[col_exp]]),
                   flag = if (length(flag_col) > 0)
                     as.logical(constraint_raw[[flag_col[1]]]) else NA) |>
  dplyr::filter(!is.na(loeuf))

if (length(flag_col) > 0 && any(constraint$flag, na.rm = TRUE)) {
  constraint <- constraint |> dplyr::filter(flag %in% TRUE)
  cat(sprintf("Filtered to %s transcripts: %d rows\n", flag_col[1], nrow(constraint)))
}
constraint <- constraint |>
  arrange(symbol, desc(exp_lof)) |>
  distinct(symbol, .keep_all = TRUE)

# ---- Map constraint onto the tested universe -------------------------------
# Hits carry gene symbols; match on symbol, fall back to Ensembl ID.
by_sym <- constraint |> dplyr::select(symbol, loeuf, exp_lof) |>
  dplyr::rename(key = symbol)
by_ens <- constraint |> dplyr::filter(!is.na(ensembl), ensembl != "") |>
  dplyr::select(ensembl, loeuf, exp_lof) |> dplyr::rename(key = ensembl)
lookup <- bind_rows(by_sym, by_ens) |> distinct(key, .keep_all = TRUE)

# Mean expression across pseudobulk samples -- the power covariate.
pb_base  <- readRDS(file.path(cache_dir, "pb_base.rds"))
cpm      <- pb_base$counts / rowSums(pb_base$counts) * 1e6
mean_expr <- log2(colMeans(cpm) + 1)

gene_df <- tibble(gene = background_genes) |>
  left_join(lookup, by = c("gene" = "key")) |>
  left_join(hit_recur |> dplyr::select(gene, n_metabolites), by = "gene") |>
  mutate(n_metabolites = replace_na(n_metabolites, 0L),
         is_hit    = n_metabolites > 0,
         is_hub    = n_metabolites >= MIN_METAB_RECURRENCE,
         mean_expr = unname(mean_expr[gene]))

cat(sprintf("\nLOEUF matched: %d / %d tested genes (%.1f%%)\n",
            sum(!is.na(gene_df$loeuf)), nrow(gene_df),
            100 * mean(!is.na(gene_df$loeuf))))
cat(sprintf("  Match rate among hit genes: %.1f%% | non-hit genes: %.1f%%\n",
            100 * mean(!is.na(gene_df$loeuf[gene_df$is_hit])),
            100 * mean(!is.na(gene_df$loeuf[!gene_df$is_hit]))))

# Unmatched genes are non-random (novel/short/non-coding loci). If the match
# rate differs between hits and non-hits, dropping them is itself a filter.
mr_test <- fisher.test(table(gene_df$is_hit, !is.na(gene_df$loeuf)))
cat(sprintf("  Differential match rate (hit vs non-hit): p = %.3g%s\n",
            mr_test$p.value,
            if (mr_test$p.value < 0.05) "  <-- INTERPRET WITH CARE" else ""))

analysis_df <- gene_df |>
  dplyr::filter(!is.na(loeuf), !is.na(mean_expr))
low_power <- sum(analysis_df$exp_lof < MIN_EXP_LOF, na.rm = TRUE)
cat(sprintf("  Genes with < %d expected pLoF (LOEUF unstable): %d (%.1f%%)\n",
            MIN_EXP_LOF, low_power, 100 * low_power / nrow(analysis_df)))

# ---- Helper: continuous LOEUF comparison -----------------------------------
# Wilcoxon rank-sum, one-sided: is LOEUF LOWER (more constrained) in fg?
compare_loeuf <- function(df, fg_col, label_fg, label_bg, restrict_exp = FALSE) {
  d <- if (restrict_exp) dplyr::filter(df, exp_lof >= MIN_EXP_LOF) else df
  fg <- d$loeuf[d[[fg_col]]]
  bg <- d$loeuf[!d[[fg_col]]]
  if (length(fg) < 3 || length(bg) < 3) return(NULL)

  wt <- wilcox.test(fg, bg, alternative = "less", conf.int = TRUE)
  cat(sprintf("\n  %s (n=%d): median LOEUF %.3f | %s (n=%d): median %.3f\n",
              label_fg, length(fg), median(fg), label_bg, length(bg), median(bg)))
  cat(sprintf("  Wilcoxon (%s more constrained): W = %.0f, p = %.3g\n",
              label_fg, wt$statistic, wt$p.value))
  cat(sprintf("  Hodges-Lehmann shift = %.3f  [%.3f, %.3f]\n",
              wt$estimate, wt$conf.int[1], wt$conf.int[2]))

  tibble(comparison = sprintf("%s vs %s", label_fg, label_bg),
         exp_lof_filtered = restrict_exp,
         n_fg = length(fg), n_bg = length(bg),
         median_fg = median(fg), median_bg = median(bg),
         hl_shift = unname(wt$estimate),
         ci_lo = wt$conf.int[1], ci_hi = wt$conf.int[2],
         p = wt$p.value)
}

# ---- Helper: LOEUF threshold ladder (secondary) ----------------------------
loeuf_ladder <- function(df, fg_col, test_label) {
  purrr::map_dfr(seq_len(nrow(LOEUF_LADDER)), function(i) {
    thr <- LOEUF_LADDER$threshold[i]
    constrained <- df$loeuf < thr
    fg <- df[[fg_col]]
    mat <- matrix(c(sum( fg &  constrained), sum(!fg &  constrained),
                    sum( fg & !constrained), sum(!fg & !constrained)), nrow = 2)
    ft1 <- fisher.test(mat, alternative = "greater")
    ft2 <- fisher.test(mat)
    tibble(test = test_label,
           threshold = thr, percentile = LOEUF_LADDER$percentile[i],
           n_constrained = sum(constrained),
           n_fg_constrained = sum(fg & constrained),
           pct_fg_constrained = 100 * sum(fg & constrained) / max(sum(constrained), 1),
           OR = unname(ft1$estimate),
           ci_lo = ft2$conf.int[1], ci_hi = ft2$conf.int[2],
           p = ft1$p.value)
  })
}

# ---- TEST 1: constraint among hit genes ------------------------------------
cat("\n##########################################################\n")
cat("TEST 1: Are hit genes more constrained than tested genes?\n")
cat("##########################################################\n")

cont_t1     <- compare_loeuf(analysis_df, "is_hit", "Hit", "Non-hit")
cont_t1_exp <- compare_loeuf(analysis_df, "is_hit", "Hit", "Non-hit", restrict_exp = TRUE)
ladder_t1   <- loeuf_ladder(analysis_df, "is_hit", "Test 1")

cat("\n  LOEUF threshold ladder (secondary):\n")
for (i in seq_len(nrow(ladder_t1)))
  cat(sprintf("    LOEUF < %.2f (top %2d%%): OR = %.2f [%.2f, %.2f], p = %.3g\n",
              ladder_t1$threshold[i], ladder_t1$percentile[i],
              ladder_t1$OR[i], ladder_t1$ci_lo[i], ladder_t1$ci_hi[i], ladder_t1$p[i]))

# --- (b) adjusted model: does LOEUF survive expression + gene size? ---------
cat("\n  Confound-adjusted logistic regression:\n")
model_df <- analysis_df |>
  dplyr::filter(exp_lof > 0) |>
  mutate(log_exp_lof = log10(exp_lof))

fit_raw <- glm(is_hit ~ loeuf, data = model_df, family = binomial)
fit_adj <- glm(is_hit ~ loeuf + mean_expr + log_exp_lof,
               data = model_df, family = binomial)

report_coef <- function(fit, label) {
  s  <- summary(fit)$coefficients
  b  <- s["loeuf", "Estimate"]; se <- s["loeuf", "Std. Error"]
  cat(sprintf("    %-28s LOEUF beta = %+.3f (SE %.3f), OR/unit = %.2f, p = %.3g\n",
              label, b, se, exp(b), s["loeuf", "Pr(>|z|)"]))
  tibble(model = label, beta = b, se = se, or_per_unit = exp(b),
         p = s["loeuf", "Pr(>|z|)"])
}
glm_res <- bind_rows(report_coef(fit_raw, "Unadjusted"),
                     report_coef(fit_adj, "+ expression + gene size"))
cat("    (OR/unit < 1 = higher LOEUF lowers hit odds, i.e. constraint enriched)\n")

# --- (c) expression-decile-matched background ------------------------------
# Draw a non-hit set with the same mean-expression profile as the hits, and
# ask how often its median LOEUF is at or below the observed hit median.
cat("\n  Expression-matched resampling:\n")
matched_df <- analysis_df |>
  mutate(expr_decile = dplyr::ntile(mean_expr, 10))
hit_med   <- median(matched_df$loeuf[matched_df$is_hit])
need      <- matched_df |> dplyr::filter(is_hit) |> count(expr_decile, name = "k")
pool      <- matched_df |> dplyr::filter(!is_hit)

null_meds <- replicate(N_MATCH_REPS, {
  idx <- purrr::map(seq_len(nrow(need)), function(i) {
    p <- which(pool$expr_decile == need$expr_decile[i])
    if (length(p) == 0) return(integer(0))
    sample(p, min(need$k[i], length(p)), replace = FALSE)
  }) |> unlist()
  median(pool$loeuf[idx])
})
p_match <- (sum(null_meds <= hit_med) + 1) / (N_MATCH_REPS + 1)
cat(sprintf("    Hit median LOEUF = %.3f | matched-null median = %.3f [%.3f, %.3f]\n",
            hit_med, median(null_meds),
            quantile(null_meds, 0.025), quantile(null_meds, 0.975)))
cat(sprintf("    Empirical one-sided p = %.4f (%d reps)\n", p_match, N_MATCH_REPS))

# ---- TEST 2: constraint among hub genes ------------------------------------
cat("\n##########################################################\n")
cat("TEST 2: Among hit genes, are hubs more constrained?\n")
cat("##########################################################\n")

hit_only    <- analysis_df |> dplyr::filter(is_hit)
cont_t2     <- compare_loeuf(hit_only, "is_hub", "Hub", "Non-hub hit")
ladder_t2   <- loeuf_ladder(hit_only, "is_hub", "Test 2")

cat("\n  LOEUF threshold ladder (secondary):\n")
for (i in seq_len(nrow(ladder_t2)))
  cat(sprintf("    LOEUF < %.2f (top %2d%%): OR = %.2f [%.2f, %.2f], p = %.3g\n",
              ladder_t2$threshold[i], ladder_t2$percentile[i],
              ladder_t2$OR[i], ladder_t2$ci_lo[i], ladder_t2$ci_hi[i], ladder_t2$p[i]))

# ---- TEST 3: LOEUF vs metabolite recurrence --------------------------------
cat("\n##########################################################\n")
cat("TEST 3: Does LOEUF track metabolite recurrence among hits?\n")
cat("##########################################################\n")

sp <- cor.test(hit_only$loeuf, hit_only$n_metabolites,
               method = "spearman", alternative = "less", exact = FALSE)
cat(sprintf("  Spearman rho = %+.3f, p = %.3g (n = %d hit genes)\n",
            sp$estimate, sp$p.value, nrow(hit_only)))
cat("  (negative rho = more constrained genes track more metabolites)\n")

# ---- SFARI overlap check ---------------------------------------------------
# SFARI genes are strongly constrained as a class, so a constraint signal could
# just be the SFARI signal restated. Re-run Test 1 with SFARI genes removed.
if (file.exists(sfari_path)) {
  sfari <- read_csv(sfari_path, show_col_types = FALSE)
  sfari_ids <- c(sfari[["gene-symbol"]], sfari[["ensembl-id"]]) |> na.omit() |> unique()
  cat("\n--- Test 1 excluding SFARI genes (is constraint a separate signal?) ---\n")
  cont_t1_nosfari <- compare_loeuf(
    analysis_df |> dplyr::filter(!gene %in% sfari_ids),
    "is_hit", "Hit (non-SFARI)", "Non-hit (non-SFARI)")
} else {
  cat("\nSFARI file not found; skipping the SFARI-excluded sensitivity check.\n")
  cont_t1_nosfari <- NULL
}

# ---- Visualizations --------------------------------------------------------

# A) LOEUF ECDF: the whole distribution, no threshold. Separation anywhere
# along the curve is the result; the ladder below just samples this picture.
ecdf_df <- analysis_df |>
  mutate(group = ifelse(is_hit, "Hit genes", "Non-hit tested genes"))

p_ecdf <- ggplot(ecdf_df, aes(x = loeuf, color = group)) +
  stat_ecdf(linewidth = 1) +
  geom_vline(xintercept = PRIMARY_THRESHOLD, linetype = "dashed", color = "grey50") +
  annotate("text", x = PRIMARY_THRESHOLD, y = 0.05, hjust = -0.08, size = 3,
           color = "grey40", label = sprintf("LOEUF %.2f", PRIMARY_THRESHOLD)) +
  scale_color_manual(values = c("Hit genes" = "#E15759",
                                "Non-hit tested genes" = "#76B7B2"), name = NULL) +
  labs(x = "LOEUF (lower = more pLoF constrained)",
       y = "Cumulative fraction of genes",
       title = "pLoF constraint: hit genes vs tested background",
       subtitle = sprintf("Wilcoxon p = %.3g; median %.3f vs %.3f",
                          cont_t1$p, cont_t1$median_fg, cont_t1$median_bg)) +
  theme_cowplot(12) +
  theme(legend.position = c(0.55, 0.25))

save_dual_format(p_ecdf, out_dir, "loeuf_ecdf_hit_vs_background", width = 8, height = 6)
print(p_ecdf)

# B) OR across the whole threshold ladder, both tests. Shown as a ladder
# precisely so no single cutoff carries the claim.
ladder_all <- bind_rows(ladder_t1, ladder_t2)

p_ladder <- ggplot(ladder_all,
                   aes(x = OR, y = factor(sprintf("< %.2f (top %d%%)", threshold, percentile),
                                          levels = rev(sprintf("< %.2f (top %d%%)",
                                                               LOEUF_LADDER$threshold,
                                                               LOEUF_LADDER$percentile))),
                       color = p < 0.05)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.6) +
  geom_point(size = 3) +
  facet_wrap(~ test, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(
               "Test 1" = "Test 1: hit genes vs tested universe",
               "Test 2" = "Test 2: hub genes vs all hit genes"))) +
  scale_x_log10() +
  scale_color_manual(values = c(`TRUE` = "#E15759", `FALSE` = "grey55"),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."), name = NULL) +
  labs(x = "Odds ratio (log scale)", y = "LOEUF threshold",
       title = "Constraint enrichment across the LOEUF ladder",
       caption = "Secondary analysis. Primary result is the continuous Wilcoxon test.") +
  theme_cowplot(12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92", color = NA))

save_dual_format(p_ladder, out_dir, "loeuf_threshold_ladder", width = 8, height = 7)
print(p_ladder)

# C) Hit rate by LOEUF decile, split by expression tertile. If the constraint
# gradient is real it holds WITHIN expression strata; if the lines flatten
# inside strata, the marginal trend was expression.
decile_df <- analysis_df |>
  mutate(loeuf_decile = dplyr::ntile(loeuf, 10),
         expr_tertile = factor(dplyr::ntile(mean_expr, 3),
                               labels = c("Low expr", "Mid expr", "High expr"))) |>
  group_by(expr_tertile, loeuf_decile) |>
  summarise(hit_rate = 100 * mean(is_hit), n = n(), .groups = "drop")

p_decile <- ggplot(decile_df, aes(x = loeuf_decile, y = hit_rate, color = expr_tertile)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c("Low expr" = "#BAB0AC", "Mid expr" = "#F28E2B",
                                "High expr" = "#4E79A7"), name = NULL) +
  labs(x = "LOEUF decile (1 = most constrained)", y = "Hit rate (%)",
       title = "Hit rate by constraint, stratified by expression",
       subtitle = "A real constraint effect persists within each expression tertile") +
  theme_cowplot(12) +
  theme(legend.position = "bottom")

save_dual_format(p_decile, out_dir, "loeuf_hitrate_by_expression", width = 8, height = 6)
print(p_decile)

# ---- Save results ----------------------------------------------------------
continuous_all <- bind_rows(cont_t1, cont_t1_exp, cont_t2, cont_t1_nosfari)
write_csv(continuous_all, file.path(out_dir, "loeuf_continuous_tests.csv"))
write_csv(ladder_all,     file.path(out_dir, "loeuf_threshold_ladder.csv"))
write_csv(glm_res,        file.path(out_dir, "loeuf_adjusted_models.csv"))
write_csv(gene_df,        file.path(out_dir, "gene_constraint_annotated.csv"))

# ---- Summary ---------------------------------------------------------------
cat("\n========== SUMMARY ==========\n")
cat(sprintf("Method: %s | %d tested genes, %d with LOEUF\n",
            METHOD, nrow(gene_df), nrow(analysis_df)))
cat(sprintf("Test 1 (hits more constrained):  median %.3f vs %.3f, Wilcoxon p = %.3g\n",
            cont_t1$median_fg, cont_t1$median_bg, cont_t1$p))
cat(sprintf("  adjusted for expression + size: OR/unit = %.2f, p = %.3g\n",
            glm_res$or_per_unit[2], glm_res$p[2]))
cat(sprintf("  expression-matched resampling:  empirical p = %.4f\n", p_match))
if (!is.null(cont_t2))
  cat(sprintf("Test 2 (hubs more constrained):  median %.3f vs %.3f, Wilcoxon p = %.3g\n",
              cont_t2$median_fg, cont_t2$median_bg, cont_t2$p))
cat(sprintf("Test 3 (LOEUF vs recurrence):    rho = %+.3f, p = %.3g\n",
            sp$estimate, sp$p.value))
cat(sprintf("\nResults saved to %s\n", out_dir))
