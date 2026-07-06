#!/usr/bin/env Rscript
# check_hit_nonzero.R
# ---------------------------------------------------------------------------
# Standalone diagnostic (assumes pb and metab_z are already in the session).
# For the top-20 global-FDR metabolite x gene x cell-type hits, recompute the
# per-hit n_nonzero = number of MODELED samples (metabolite present in that
# cell type) with raw gene count >= 2 -- i.e. how many points the fit rested on.
# Also saves the shared top20_pairs.csv the scatter qmd reads (compute-once, so
# the check and the figures cannot drift).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(edgeR) })
source("metabolite_functions.R")            # metab_sample_to_pb

parquet_dir <- "../../results/gene-metabolite/parquet"
out_csv     <- "../../results/gene-metabolite/csv/top20_pairs.csv"
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

FDR_GLOBAL <- 0.10
LFC_MIN    <- 0.25
TOP_N      <- 20

# ---- global BH (memory-light; identical to the qmd) -------------------------
bh_global <- function(p, m) {
  o <- order(p); ro <- order(o)
  adj <- p[o] * m / seq_along(p[o])
  adj <- rev(cummin(rev(adj)))
  pmin(adj[ro], 1)
}

ds      <- open_dataset(parquet_dir)
n_tests <- nrow(ds)

hits <- ds |>
  filter(P.Value < 0.05) |>
  collect() |>
  mutate(FDR_global = bh_global(P.Value, n_tests)) |>
  filter(FDR_global < FDR_GLOBAL, abs(logFC) >= LFC_MIN)

# one plot per metabolite: keep each metabolite's best gene (lowest global FDR,
# ties by larger |logFC|), then take the top 20 distinct metabolites
top20 <- hits |>
  arrange(FDR_global, desc(abs(logFC))) |>
  distinct(Compound.ID, .keep_all = TRUE) |>
  slice_head(n = TOP_N)

# ---- n_nonzero over the modeled sample set ----------------------------------
rn        <- rownames(pb$counts)
row_samp  <- sub("^([^.]+)\\..*$", "\\1", rn)     # sample_id  (before first ".")
row_ct    <- sub("^[^.]+\\.(.*)$",  "\\1", rn)     # cell_type  (after first ".")
metab_col <- metab_sample_to_pb(colnames(metab_z), row_samp)  # metab column per row

nonzero_for <- function(cid, gene, ct) {
  if (!gene %in% colnames(pb$counts)) return(NA_integer_)
  idx      <- which(row_ct == ct)
  mcol     <- metab_col[idx]
  present  <- !is.na(mcol) &
              !is.na(metab_z[cid, match(mcol, colnames(metab_z))])   # modeled samples
  sum(pb$counts[idx[present], gene] >= 2)
}

# ---- influence: is the slope driven by one high-leverage sample? ------------
# Refit the model's lm form on the modeled subset, then (1) Cook's distance
# (flag if any point exceeds the 4/n rule of thumb) and (2) leave-one-out on the
# most influential point: unstable if dropping it flips the metabolite slope's
# sign or shrinks |slope| by >50%. Cook (1977) Technometrics 19:15.
influence_for <- function(cid, gene, ct) {
  na_row <- tibble(max_cooks = NA_real_, cooks_flag = NA,
                   loo_unstable = NA, n_fit = NA_integer_)
  if (!gene %in% colnames(pb$counts)) return(na_row)
  idx  <- which(row_ct == ct)
  mz   <- metab_z[cid, match(metab_col[idx], colnames(metab_z))]
  keep <- !is.na(mz)
  rows_use <- idx[keep]
  if (length(rows_use) < 4) return(na_row)

  csub <- t(pb$counts[rows_use, , drop = FALSE])
  expr <- as.numeric(log1p(edgeR::cpm(csub, log = FALSE)[gene, ]))
  md   <- pb$metadata[rn[rows_use], ]
  df   <- tibble(expr = expr, metab = as.numeric(mz[keep]),
                 GW = md$GW, Sex = md$Sex, batch1_frac = md$batch1_frac)

  covs  <- c("GW", "Sex", "batch1_frac")
  keepc <- covs[vapply(covs, function(cv) length(unique(df[[cv]])) >= 2, logical(1))]
  fml   <- reformulate(c("metab", keepc), response = "expr")
  fit   <- lm(fml, data = df)
  cd    <- cooks.distance(fit); n <- nrow(df)
  b     <- coef(fit)[["metab"]]
  b2    <- coef(lm(fml, data = df[-which.max(cd), ]))[["metab"]]

  tibble(max_cooks    = max(cd),
         cooks_flag   = max(cd) > 4 / n,
         loo_unstable = sign(b2) != sign(b) || abs(b2) < 0.5 * abs(b),
         n_fit        = n)
}

top20 <- top20 |>
  mutate(n_nonzero = pmap_int(list(Compound.ID, gene, cell_type), nonzero_for)) |>
  bind_cols(pmap_dfr(list(top20$Compound.ID, top20$gene, top20$cell_type),
                     influence_for))

write_csv(
  top20 |> dplyr::select(Compound.ID, Name, gene, cell_type, logFC, adj.P.Val,
                  FDR_global, n_samples, n_nonzero,
                  max_cooks, cooks_flag, loo_unstable),
  out_csv
)

# ---- report -----------------------------------------------------------------
cat("Top-20 global-FDR hits -- n_nonzero (modeled samples with count >= 2):\n")
print(summary(top20$n_nonzero))
cat(sprintf("\nResting on exactly 3 non-zero samples: %d / %d\n",
            sum(top20$n_nonzero == 3, na.rm = TRUE), nrow(top20)))
cat(sprintf("Cook's-flagged (influential point, >4/n): %d / %d\n",
            sum(top20$cooks_flag, na.rm = TRUE), nrow(top20)))
cat(sprintf("Slope unstable to dropping top point: %d / %d\n",
            sum(top20$loo_unstable, na.rm = TRUE), nrow(top20)))
cat(sprintf("Saved shared pair list -> %s\n\n", out_csv))

print(as.data.frame(
  top20 |> dplyr::select(Name, gene, cell_type, logFC, FDR_global,
                  n_samples, n_nonzero, max_cooks, cooks_flag, loo_unstable)
), row.names = FALSE)
