#!/usr/bin/env Rscript
# check_int_smoketest.R
# ---------------------------------------------------------------------------
# Fast validation BEFORE the 5-hour re-run: refit the current top-20 offender
# pairs (from top20_pairs.csv, built under the old z-score transform) using the
# NEW INT transform, and check whether the leave-one-out-unstable hits collapse.
# Assumes pb is already in the session (as in check_hit_nonzero.R).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(tidyverse); library(edgeR) })
source("metabolite_functions.R")     # transform_metabolite_matrix() now = INT

metab_path <- "../../results/untargeted/batch2_peak_area_clean.csv"
pairs_csv  <- "../../results/gene-metabolite/csv/top20_pairs.csv"

# ---- rebuild the predictor under the NEW transform (INT) ---------------------
metab_mat <- read_metabolite_matrix(metab_path)
metab_int <- transform_metabolite_matrix(metab_mat)      # INT (post-swap)

pairs <- readr::read_csv(pairs_csv, show_col_types = FALSE)

rn        <- rownames(pb$counts)
row_ct    <- sub("^[^.]+\\.(.*)$", "\\1", rn)
metab_col <- metab_sample_to_pb(colnames(metab_int), sub("^([^.]+)\\..*$", "\\1", rn))

# ---- same leave-one-out stability check, now on INT values ------------------
loo_unstable_int <- function(cid, gene, ct) {
  if (!gene %in% colnames(pb$counts)) return(list(unstable = NA, max_cooks = NA_real_))
  idx  <- which(row_ct == ct)
  mv   <- metab_int[cid, match(metab_col[idx], colnames(metab_int))]
  keep <- !is.na(mv)
  rows_use <- idx[keep]
  if (length(rows_use) < 4) return(list(unstable = NA, max_cooks = NA_real_))

  csub <- t(pb$counts[rows_use, , drop = FALSE])
  expr <- as.numeric(log1p(edgeR::cpm(csub, log = FALSE)[gene, ]))
  md   <- pb$metadata[rn[rows_use], ]
  df   <- tibble(expr = expr, metab = as.numeric(mv[keep]),
                 GW = md$GW, Sex = md$Sex, batch1_frac = md$batch1_frac)

  covs  <- c("GW", "Sex", "batch1_frac")
  keepc <- covs[vapply(covs, function(cv) length(unique(df[[cv]])) >= 2, logical(1))]
  fml   <- reformulate(c("metab", keepc), response = "expr")
  fit   <- lm(fml, data = df)
  cd    <- cooks.distance(fit)
  b     <- coef(fit)[["metab"]]
  b2    <- coef(lm(fml, data = df[-which.max(cd), ]))[["metab"]]
  list(unstable = sign(b2) != sign(b) || abs(b2) < 0.5 * abs(b),
       max_cooks = max(cd))
}

res <- pmap_dfr(list(pairs$Compound.ID, pairs$gene, pairs$cell_type),
                function(cid, g, ct) {
                  o <- loo_unstable_int(cid, g, ct)
                  tibble(loo_unstable_int = o$unstable, max_cooks_int = o$max_cooks)
                })

cmp <- pairs |>
  dplyr::select(Name, gene, cell_type,
         max_cooks_z = max_cooks, loo_unstable_z = loo_unstable) |>
  bind_cols(res)

# ---- report -----------------------------------------------------------------
was_unstable <- cmp |> filter(loo_unstable_z %in% TRUE)
fixed        <- sum(was_unstable$loo_unstable_int %in% FALSE)

cat(sprintf("Previously loo-unstable (z-score): %d\n", nrow(was_unstable)))
cat(sprintf("Now stable under INT:             %d / %d\n", fixed, nrow(was_unstable)))
cat(sprintf("Still unstable under INT:         %d\n",
            sum(cmp$loo_unstable_int %in% TRUE)))
cat(sprintf("Max Cook's D  z: %.1f  ->  INT: %.1f (max across pairs)\n\n",
            max(cmp$max_cooks_z, na.rm = TRUE),
            max(cmp$max_cooks_int, na.rm = TRUE)))

print(as.data.frame(cmp), row.names = FALSE)
