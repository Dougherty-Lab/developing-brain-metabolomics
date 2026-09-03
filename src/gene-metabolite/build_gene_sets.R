#!/usr/bin/env Rscript
# build_gene_sets.R
# ---------------------------------------------------------------------------
# Parse doc/gene_lists/disease-associated-genes.xlsx into ONE tidy long CSV
# used as the input to geneset_enrichment.R. Every dedup / stacked-table /
# column-choice decision lives here rather than in a hand-edited file, so the
# gene sets reported in the manuscript are reproducible from the source
# supplementary tables.
#
# Sheets are heterogeneous (one gene per row, comma-delimited gene lists,
# boolean membership columns, two tables stacked in one sheet), so each gets a
# bespoke reader. Counts are asserted against the published figures at the end;
# the script stops if any set has drifted.
#
# OUTPUT: doc/gene_lists/disease_gene_sets.csv, columns
#   set_id       machine key, e.g. "scz_trubetskoy"
#   set_label    display label for figures
#   condition    SCZ | NDD | Epilepsy | T2D | IBD
#   source       first-author-year
#   level        "paper" (one publication) | "condition" (union of papers)
#   set_type     "neuro" (test set) | "control" (non-neural negative control)
#   gene_symbol
#
# Condition-level union rows are emitted ONLY where a condition has >1 paper
# (SCZ, Epilepsy). For single-paper conditions the paper-level set already IS
# the condition set, and duplicating it would add identical tests to the
# multiple-testing grid for no information.
#
# Run from src/gene-metabolite/:
#   Rscript build_gene_sets.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)      # NOT currently in the .sif -- add at container rebuild
})

xlsx_path <- "../../doc/gene_lists/disease_associated_genes.xlsx"
out_path  <- "../../doc/gene_lists/disease_gene_sets.csv"

stopifnot(file.exists(xlsx_path))

# ---- Helpers ---------------------------------------------------------------
# Read one sheet and pull one column BY NAME, failing loudly with the available
# names if it is missing. Sheet layouts change between journal supplements, and
# a silent zero-length vector here is much harder to debug than an error.
read_col <- function(sheet, column, skip = 0, col_index = NULL,
                     col_names = TRUE) {
  df <- suppressMessages(
    read_excel(xlsx_path, sheet = sheet, skip = skip, col_names = col_names)
  )
  if (!is.null(col_index)) {
    if (ncol(df) < col_index)
      stop(sprintf("Sheet '%s': expected at least %d columns, found %d.",
                   sheet, col_index, ncol(df)))
    v <- df[[col_index]]
  } else {
    if (!column %in% names(df))
      stop(sprintf("Sheet '%s': no column '%s'. Available: %s",
                   sheet, column, paste(names(df), collapse = " | ")))
    v <- df[[column]]
  }
  cat(sprintf("  read %-18s %-22s %d rows\n", sheet,
              if (is.null(col_index)) column else paste0("col#", col_index),
              length(v)))
  v
}

# Cells contain non-breaking spaces (\u00a0) and em-spaces (\u2003) from the
# publisher's HTML. Splitting on "," alone leaves those bound to the symbol and
# every downstream join silently misses.
#
# Base subsetting rather than purrr::discard(): the vectorised comparison is
# unambiguous, whereas a predicate lambda applied to an atomic vector is easy
# to get subtly wrong and fails by returning nothing rather than by erroring.
split_genes <- function(x) {
  v <- as.character(x)
  v <- gsub("[\u00a0\u2003\u2002\u2009]", " ", v)
  v <- unlist(strsplit(v, "[,;/]"), use.names = FALSE)
  v <- trimws(v)
  v <- v[!is.na(v) & !(v %in% c("", "NA", "-", "\u2013", "\u2014"))]
  unique(v)
}

tidy_set <- function(genes, set_id, set_label, condition, source, set_type,
                     level = "paper") {
  genes <- sort(unique(genes))
  if (length(genes) == 0)
    stop(sprintf("Set '%s' parsed to zero genes -- check the sheet and column.",
                 set_id))
  tibble(
    set_id      = set_id,
    set_label   = set_label,
    condition   = condition,
    source      = source,
    level       = level,
    set_type    = set_type,
    gene_symbol = genes
  )
}

cat("Reading sheets:\n")

# ---- Trubetskoy et al. 2022 Nature -- SCZ GWAS prioritised genes -----------
# One gene per row; the FINEMAP / SMR / Rare priority columns are 0/1 and are
# NOT used to subset -- the published set is all 120 rows.
trubetskoy <- split_genes(read_col("Trubetskoy et al", "Symbol.ID"))

# ---- Chick et al. 2025 Nat Commun -- 8 additional SCZ genes ----------------
# Two header rows (row 1 = analysis blocks, row 2 = statistic names), so data
# starts at row 3. col_names = FALSE is essential: without it readxl treats
# row 3 as the header and the first gene (STAG1) is silently lost, giving 7.
# Taken by position because there is no header name to match on.
chick <- split_genes(read_col("Chick et al", NULL, skip = 2, col_index = 1,
                              col_names = FALSE))

# ---- Fu et al. 2022 Nat Genet -- developmental disorder genes --------------
# 18,160 rows x boolean membership columns. DD477 = the FDR < 0.05 DD set.
# ASD185 is deliberately NOT used: those genes feed into SFARI, and testing
# both would double-count the same evidence (SFARI is handled entirely in
# sfari_enrichment.R). SCZ10 / SCZ244 are also skipped -- two SCZ sources with
# independent ascertainment are already in hand.
fu_raw <- suppressMessages(read_excel(xlsx_path, sheet = "Fu et al"))
for (need in c("gene", "DD477"))
  if (!need %in% names(fu_raw))
    stop(sprintf("Sheet 'Fu et al': no column '%s'. Available: %s",
                 need, paste(head(names(fu_raw), 30), collapse = " | ")))
# The membership columns may arrive as logical or as the strings "TRUE"/"FALSE"
# depending on how the supplement was exported; accept both.
fu <- split_genes(fu_raw$gene[fu_raw$DD477 %in% c(TRUE, "TRUE")])
cat(sprintf("  read %-18s %-22s %d rows\n", "Fu et al", "gene[DD477]", length(fu)))

# ---- ILAE Consortium 2023 Nat Genet -- epilepsy GWAS -----------------------
ilae <- split_genes(read_col("ILAEC", "Genes"))

# ---- Zhang et al. 2024 Seizure -- clinically curated epilepsy genes --------
# Genes are listed per phenotype across two columns ("Genes by 2017",
# "Genes updated"), comma-delimited, and RECUR across phenotype rows (SCN1A
# appears in 4 rows, GABRG2 in 4, HCN2 in 3). The paper's "168" counts
# gene-phenotype entries; the unique gene set is 133. A gene set is a set, so
# 133 is the correct number to test and to report.
zhang <- split_genes(c(read_col("Zhang et al", "Genes by 2017"),
                       read_col("Zhang et al", "Genes updated")))

# ---- Xue et al. 2018 Nat Commun -- T2D (NEGATIVE CONTROL) ------------------
# This sheet holds TWO stacked tables: the paper's Table 2 and Table 3, with a
# repeated header row partway down. Dropping rows where Gene == "Gene" removes
# that interior header. The two tables share 9 genes; the union is 33.
xue_raw <- read_col("Xue et al", "Gene")
xue <- split_genes(xue_raw[!is.na(xue_raw) & !(xue_raw %in% c("Gene", "probe ID"))])

# ---- Alegbe et al. 2026 Nature -- IBD (NEGATIVE CONTROL) -------------------
# First column header is blank, so readxl repairs it to `...1`; harmless, and
# the message is suppressed in read_col().
alegbe <- split_genes(read_col("Alegbe et al", "Gene symbol"))

# ---- Aragam et al. 2022 Nature Genetics -- CAD (NEGATIVE CONTROL) ----------
# Supplementary Table 31: 279 genome-wide CAD associations, each with a
# prioritised causal gene from an eight-predictor framework. Rows 1-3 are the
# table caption and a merged grouping header, so the real header is row 4
# (skip = 3).
#
# Two quirks specific to this sheet, both absorbed by split_genes():
#   - 23 rows list TIED genes slash-delimited (SORT1/CELSR2, LPA/PLG). All tied
#     genes are kept: the paper prioritises none of them over the others, and
#     dropping either would be an arbitrary choice.
#   - Associations where no gene had more than one supporting predictor carry
#     "-" rather than a symbol; these are dropped, not counted as a gene.
# 279 rows collapse to 220 unique genes, since a gene can be prioritised at
# several nearby associations (PCSK9 appears three times, LPA/PLG four).
#
# This is the LARGE non-neural control. T2D (33) and IBD (26) are small enough
# that a null result could be dismissed as low power; at 220 genes CAD is
# comparable in size to the DD and epilepsy sets, so a null here is
# interpretable rather than merely uninformative.
aragam <- split_genes(read_col("Aragam et al", "most_likely_causal_gene",
                               skip = 3))

cat("\nParsed set sizes before assembly:\n")
print(c(trubetskoy = length(trubetskoy), chick = length(chick),
        fu = length(fu), ilae = length(ilae), zhang = length(zhang),
        xue = length(xue), alegbe = length(alegbe),
        aragam = length(aragam)))

# ---- Assemble --------------------------------------------------------------
paper_sets <- bind_rows(
  tidy_set(trubetskoy, "scz_trubetskoy", "SCZ (Trubetskoy 2022)",
           "SCZ",      "Trubetskoy 2022", "neuro"),
  tidy_set(chick,      "scz_chick",      "SCZ (Chick 2025)",
           "SCZ",      "Chick 2025",      "neuro"),
  tidy_set(fu,         "ndd_fu_dd477",   "DD (Fu 2022, DD477)",
           "NDD",      "Fu 2022",         "neuro"),
  tidy_set(ilae,       "epi_ilae",       "Epilepsy (ILAE 2023)",
           "Epilepsy", "ILAE 2023",       "neuro"),
  tidy_set(zhang,      "epi_zhang",      "Epilepsy (Zhang 2024)",
           "Epilepsy", "Zhang 2024",      "neuro"),
  tidy_set(xue,        "t2d_xue",        "T2D (Xue 2018)",
           "T2D",      "Xue 2018",        "control"),
  tidy_set(alegbe,     "ibd_alegbe",     "IBD (Alegbe 2026)",
           "IBD",      "Alegbe 2026",     "control"),
  tidy_set(aragam,     "cad_aragam",     "CAD (Aragam 2022)",
           "CAD",      "Aragam 2022",     "control")
)

# Condition-level unions, only where a condition draws on >1 paper.
multi_paper <- paper_sets |>
  distinct(condition, source) |>
  count(condition) |>
  filter(n > 1) |>
  pull(condition)

condition_sets <- paper_sets |>
  filter(condition %in% multi_paper) |>
  group_by(condition, set_type) |>
  summarise(gene_symbol = list(sort(unique(gene_symbol))), .groups = "drop") |>
  mutate(set_id    = paste0(tolower(condition), "_union"),
         set_label = paste0(condition, " (all sources)"),
         source    = paste0("union: ",
                            map_chr(condition, ~ paste(sort(unique(
                              paper_sets$source[paper_sets$condition == .x])),
                              collapse = " + "))),
         level     = "condition") |>
  unnest(gene_symbol) |>
  dplyr::select(set_id, set_label, condition, source, level, set_type, gene_symbol)

gene_sets <- bind_rows(paper_sets, condition_sets) |>
  arrange(set_type, condition, level, set_id, gene_symbol)

# ---- Verify against published counts ---------------------------------------
# Hard stop rather than a warning: a silently changed set size means every
# downstream odds ratio is wrong, and the manuscript quotes these numbers.
expected <- c(scz_trubetskoy = 120, scz_chick = 8, ndd_fu_dd477 = 477,
              epi_ilae = 29, epi_zhang = 133, t2d_xue = 33, ibd_alegbe = 26,
              cad_aragam = 220)

observed <- gene_sets |> filter(level == "paper") |> count(set_id) |> deframe()

for (s in names(expected)) {
  if (is.na(observed[s]) || observed[s] != expected[s]) {
    stop(sprintf("Set '%s': expected %d genes, parsed %s. Sheet layout changed.",
                 s, expected[s], ifelse(is.na(observed[s]), "0", observed[s])))
  }
}

cat("Gene sets built:\n")
print(gene_sets |>
        group_by(set_id, set_label, condition, level, set_type) |>
        summarise(n_genes = n_distinct(gene_symbol), .groups = "drop") |>
        arrange(set_type, condition, level) |>
        as.data.frame(), right = FALSE)

# Overlap between neuro sets, reported because the tests are NOT independent.
neuro <- gene_sets |> filter(set_type == "neuro", level == "paper")
pairs <- combn(unique(neuro$set_id), 2, simplify = FALSE)
overlap <- map_dfr(pairs, function(p) {
  a <- neuro$gene_symbol[neuro$set_id == p[1]]
  b <- neuro$gene_symbol[neuro$set_id == p[2]]
  tibble(set_a = p[1], set_b = p[2],
         n_shared = length(intersect(a, b)),
         jaccard  = length(intersect(a, b)) / length(union(a, b)))
}) |> arrange(desc(n_shared))

cat("\nPairwise overlap between neuro paper-level sets:\n")
print(as.data.frame(overlap), right = FALSE)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_csv(gene_sets, out_path)
cat(sprintf("\nWrote %s (%d rows, %d sets)\n",
            out_path, nrow(gene_sets), n_distinct(gene_sets$set_id)))

cat("\nSession info:\n"); print(sessionInfo())
