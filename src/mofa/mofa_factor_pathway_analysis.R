# mofa_factor_pathway_analysis.R
# ---------------------------------------------------------------------------
# Metabolite pathway analysis (ORA + signed MSEA) on MOFA factor weights.
# Default: Factor 2, metabolite view.
#
# Inputs:  results/mofa/mofa_model.hdf5, results/untargeted/metabolite_annotations.csv,
#          doc/untargeted/combined_annotations_untargeted_developing_2026.xlsx
# Outputs: Figures and CSVs in results/mofa/pathway/
#
# Upstream:  mofa_integration.Rmd
# Downstream: None
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(fgsea)
  library(KEGGREST)
  library(cowplot)
  library(ggrepel)
  library(svglite)
})

source(file.path(root, "src/gene-metabolite/pseudobulk_functions.R"))

set.seed(42)

# ---- project root ------------------------------------------------------------
# Anchor output paths to the repo root (.git marker) instead of a cwd-relative
# "../.." — running this interactively can leave the working directory
# somewhere unexpected, which otherwise silently writes results under src/.
find_project_root <- function(marker = ".git") {
  d <- normalizePath(getwd())
  repeat {
    if (dir.exists(file.path(d, marker))) return(d)
    parent <- dirname(d)
    if (parent == d) {
      stop("Could not locate project root (no ", marker, " found above ", getwd(), ")")
    }
    d <- parent
  }
}
root        <- find_project_root()
results_dir <- file.path(root, "results")
doc_dir     <- file.path(root, "doc")

# ==============================================================================
# 0. Parameters
# ==============================================================================

FACTOR        <- 2L                 # MOFA factor to analyse
VIEW          <- "metabolites"      # MOFA view
TOP_N_SET     <- c(25L, 50L, 100L)  # cutoffs to run
TOP_N_PRIMARY <- 50L                # cutoff reported in the manuscript
MIN_SET       <- 5L                 # min pathway members (after background intersect)
MAX_SET       <- 500L
NPERM         <- 10000L

factor_name <- paste0("Factor", FACTOR)

mofa_dir <- file.path(results_dir, "mofa")
out_dir  <- file.path(mofa_dir, "pathway", tolower(factor_name))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Project root : %s\n", root))
cat(sprintf("Target       : %s, view = %s\n", factor_name, VIEW))
cat(sprintf("Output        : %s\n\n", out_dir))

# ==============================================================================
# 1. Load MOFA weights
# ==============================================================================

weights_all <- read_csv(file.path(mofa_dir, "mofa_metabolite_weights.csv"),
                        show_col_types = FALSE)

# get_weights(as.data.frame = TRUE) returns: view, feature, factor, value
stopifnot(all(c("view", "feature", "factor", "value") %in% names(weights_all)))

w_fac <- weights_all |>
  dplyr::filter(view == VIEW, factor == factor_name) |>
  dplyr::select(feature, weight = value) |>
  dplyr::filter(!is.na(weight)) |>
  dplyr::arrange(dplyr::desc(abs(weight)))

cat(sprintf("Features with a %s weight: %d\n", factor_name, nrow(w_fac)))
cat(sprintf("Weight range: %.3f to %.3f | %d positive, %d negative\n\n",
            min(w_fac$weight), max(w_fac$weight),
            sum(w_fac$weight > 0), sum(w_fac$weight < 0)))

# ==============================================================================
# 2. Sign orientation — what does "+" mean for this factor?
#
# MOFA factor signs are arbitrary (the model is invariant to jointly flipping
# a factor's values and its weights), so "positive weight" is meaningless until
# anchored to something external. Print the factor's GW / hormone associations
# so the +/- ORA results can be read directionally.
# ==============================================================================

cat("=== Sign orientation for", factor_name, "===\n")

fac_path  <- file.path(mofa_dir, "mofa_factor_values.csv")
horm_path <- file.path(mofa_dir, "mofa_hormone_correlations.csv")

if (file.exists(fac_path)) {
  fac_vals <- read_csv(fac_path, show_col_types = FALSE)
  fv <- fac_vals |> dplyr::filter(factor == factor_name)

  if (nrow(fv) > 0 && "GW" %in% names(fv)) {
    ct_gw <- suppressWarnings(cor.test(fv$value, fv$GW, method = "spearman"))
    cat(sprintf("  %s vs GW : rho = %+.3f, p = %.3g (n = %d)\n",
                factor_name, ct_gw$estimate, ct_gw$p.value, sum(!is.na(fv$GW))))
    cat(sprintf("  -> a POSITIVE weight means the metabolite is higher in samples with %s GW\n",
                ifelse(ct_gw$estimate > 0, "higher", "lower")))
  }
  if (nrow(fv) > 0 && "Sex" %in% names(fv)) {
    cat("  Factor value by sex:\n")
    print(fv |> dplyr::group_by(Sex) |>
            dplyr::summarise(median = median(value, na.rm = TRUE),
                             n = dplyr::n(), .groups = "drop") |>
            as.data.frame(), row.names = FALSE)
  }
} else {
  cat("  mofa_factor_values.csv not found — skipping GW/Sex orientation.\n")
}

if (file.exists(horm_path)) {
  horm <- read_csv(horm_path, show_col_types = FALSE) |>
    dplyr::filter(if ("factor" %in% names(read_csv(horm_path, n_max = 1,
                                                   show_col_types = FALSE)))
      factor == factor_name else TRUE)
  cat("  Post-hoc hormone correlations:\n")
  print(as.data.frame(horm), row.names = FALSE)
} else {
  cat("  mofa_hormone_correlations.csv not found — skipping hormone orientation.\n")
}
cat("\n")

# ==============================================================================
# 3. Reconcile MOFA feature names back to annotation `Name`
#
# mofa_integration.qmd renamed the metabolite view rows from Compound.ID to
# transliterated Names (alpha/beta/gamma/delta/upsilon, +/- for the plus-minus
# sign, iconv ASCII//TRANSLIT for the rest), appending " (Compound.ID)" where
# the Name was duplicated, and falling back to Compound.ID where no Name
# existed. Those strings will NOT join to the raw `Name` column, so map back:
#   (a) feature ends in " (CXXXX)"       -> use the parenthesised Compound.ID
#   (b) feature is itself a Compound.ID  -> direct
#   (c) otherwise                        -> match on transliterated Name
# ==============================================================================

annot_full <- read_csv(
  file.path(results_dir, "untargeted", "metabolite_annotations_flagged.csv"),
  show_col_types = FALSE
)

transliterate_name <- function(x) {
  x <- gsub("\u00b1", "+/-", x)
  x <- gsub("\u03b1|\u0391", "alpha-", x)
  x <- gsub("\u03b2|\u0392", "beta-",  x)
  x <- gsub("\u03b3|\u0393", "gamma-", x)
  x <- gsub("\u03b4|\u0394", "delta-", x)
  x <- gsub("\u03c5|\u03a5", "upsilon-", x)
  iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
}

annot_map <- annot_full |>
  dplyr::mutate(mofa_name = transliterate_name(Name)) |>
  dplyr::select(Compound.ID, Name, mofa_name)

# (a) strip a trailing " (Compound.ID)" disambiguator, if present
w_fac <- w_fac |>
  dplyr::mutate(
    paren_id = str_match(feature, "\\s\\(([^()]+)\\)$")[, 2],
    id_from_paren = ifelse(!is.na(paren_id) & paren_id %in% annot_map$Compound.ID,
                           paren_id, NA_character_),
    # (b) feature is itself a Compound.ID
    id_direct = ifelse(feature %in% annot_map$Compound.ID, feature, NA_character_),
    Compound.ID = dplyr::coalesce(id_direct, id_from_paren)
  )

# (c) everything still unmatched -> join on transliterated Name
name_lookup <- annot_map |>
  dplyr::distinct(mofa_name, .keep_all = TRUE) |>   # ambiguous names resolved via (a)
  dplyr::select(mofa_name, id_from_name = Compound.ID)

w_fac <- w_fac |>
  dplyr::left_join(name_lookup, by = c("feature" = "mofa_name")) |>
  dplyr::mutate(Compound.ID = dplyr::coalesce(Compound.ID, id_from_name)) |>
  dplyr::left_join(annot_map |> dplyr::select(Compound.ID, Name), by = "Compound.ID") |>
  dplyr::select(feature, Compound.ID, Name, weight)

n_unmatched <- sum(is.na(w_fac$Name))
cat(sprintf("Name reconciliation: %d / %d features matched to annotation\n",
            nrow(w_fac) - n_unmatched, nrow(w_fac)))
if (n_unmatched > 0) {
  cat("  UNMATCHED (excluded from pathway analysis):\n")
  print(w_fac |> dplyr::filter(is.na(Name)) |> dplyr::pull(feature))
}

w_fac <- w_fac |> dplyr::filter(!is.na(Name)) |> dplyr::distinct(Name, .keep_all = TRUE)

# Background = metabolites that entered the MOFA model, not all measured
background_all <- w_fac$Name
cat(sprintf("Background (MOFA metabolite view): %d metabolites\n\n",
            length(background_all)))

write_csv(w_fac, file.path(out_dir, paste0(tolower(factor_name), "_weights_annotated.csv")))

# ==============================================================================
# 4. Build pathway sets  (same construction as metabolite_pathway_analysis.qmd)
# ==============================================================================

annot_bg <- annot_full |> dplyr::filter(Name %in% background_all)

pathway_info <- read_excel(
  file.path(doc_dir, "untargeted", "combined_annotations_untargeted_developing_2026.xlsx"),
  sheet = "Pathways"
) |>
  dplyr::rename(pathway_id = PathwayId, pathway_name = pathwayName,
                pathway_source = pathwaySource) |>
  dplyr::select(pathway_id, pathway_name, pathway_source)

# ---- 4a. parse the existing `Pathways` column -------------------------------
pathways_parsed <- annot_bg |>
  dplyr::filter(!is.na(Pathways) & Pathways != "") |>
  dplyr::select(Compound.ID, Name, KEGG, Pathways) |>
  separate_rows(Pathways, sep = ";") |>
  dplyr::mutate(
    pathway_id = str_trim(Pathways),
    pathway_source = dplyr::case_when(
      str_detect(pathway_id, "^map\\d+") ~ "kegg",
      str_detect(pathway_id, "^R-HSA-")  ~ "reactome",
      str_detect(pathway_id, "^WP\\d+$") ~ "wiki",
      TRUE                               ~ NA_character_
    )
  ) |>
  dplyr::select(-Pathways) |>
  dplyr::filter(!is.na(pathway_id), pathway_id != "", !is.na(pathway_source))

cat(sprintf("Metabolites with existing pathway annotation: %d of %d\n",
            n_distinct(pathways_parsed$Name), nrow(annot_bg)))

# ---- 4b. KEGGREST re-map for metabolites with a KEGG ID but no Pathways -----
metabolites_missing <- annot_bg |>
  dplyr::filter(is.na(Pathways) | Pathways == "") |>
  dplyr::filter(!is.na(KEGG) & KEGG != "") |>
  dplyr::select(Compound.ID, Name, KEGG) |>
  separate_rows(KEGG, sep = ";") |>
  dplyr::mutate(KEGG = str_trim(KEGG),
                KEGG_valid = str_detect(KEGG, "^C\\d{5}$")) |>
  dplyr::filter(!is.na(KEGG), KEGG != "")

valid_ids <- metabolites_missing |> dplyr::filter(KEGG_valid) |>
  dplyr::pull(KEGG) |> unique()

cat(sprintf("KEGG IDs to query: %d\n", length(valid_ids)))

query_kegg_pathways <- function(kegg_ids) {
  results <- vector("list", length(kegg_ids))
  for (i in seq_along(kegg_ids)) {
    kid <- kegg_ids[i]
    results[[i]] <- tryCatch({
      links <- keggLink("pathway", paste0("cpd:", kid))
      tibble(KEGG = kid, pathway_id = str_remove(links, "path:"))
    }, error = function(e) {
      message("  Query failed for ", kid, ": ", e$message)
      tibble(KEGG = character(), pathway_id = character())
    })
    if (i %% 10 == 0) { cat("  ", i, "/", length(kegg_ids), "\n"); flush.console() }
    Sys.sleep(0.35)   # stay inside KEGG's rate limit
  }
  bind_rows(results)
}

kegg_link_df <- if (length(valid_ids) > 0) query_kegg_pathways(valid_ids) else
  tibble(KEGG = character(), pathway_id = character())

kegg_new <- if (nrow(kegg_link_df) > 0) {
  metabolites_missing |>
    dplyr::inner_join(kegg_link_df, by = "KEGG") |>
    dplyr::mutate(pathway_source = "kegg") |>
    dplyr::select(Compound.ID, Name, KEGG, pathway_source, pathway_id)
} else {
  tibble(Compound.ID = character(), Name = character(), KEGG = character(),
         pathway_source = character(), pathway_id = character())
}
cat(sprintf("New KEGG assignments: %d rows | %d metabolites\n",
            nrow(kegg_new), n_distinct(kegg_new$Name)))

# ---- 4c. combine, label, and cut to 5–500 members ---------------------------
pathways_combined <- bind_rows(
  pathways_parsed,
  kegg_new |> dplyr::select(Compound.ID, Name, pathway_source, pathway_id)
) |>
  dplyr::left_join(pathway_info |> dplyr::select(pathway_id, pathway_name),
                   by = "pathway_id") |>
  dplyr::mutate(
    pathway_label = dplyr::case_when(
      !is.na(pathway_name) & pathway_name != "" ~ pathway_name,
      pathway_source == "kegg"     ~ paste0("KEGG:", pathway_id),
      pathway_source == "reactome" ~ paste0("Reactome:", pathway_id),
      pathway_source == "wiki"     ~ paste0("WikiPathways:", pathway_id),
      pathway_source == "pfocr"    ~ paste0("PFOCR:", pathway_id),
      TRUE ~ pathway_id
    )
  )

pathway_sets <- pathways_combined |>
  dplyr::distinct(Name, pathway_id, pathway_label, pathway_source) |>
  dplyr::filter(Name %in% background_all) |>
  dplyr::mutate(pathway_source = factor(
    pathway_source,
    levels = c("kegg", "reactome", "wiki", "pfocr"),
    labels = c("KEGG", "Reactome", "WikiPathways", "PFOCR")
  )) |>
  dplyr::group_by(pathway_label, pathway_source) |>
  dplyr::summarise(members = list(unique(Name)), .groups = "drop") |>
  dplyr::mutate(
    n_members = lengths(members),
    # append source so cross-database duplicates stay distinct
    pathway_label_src = paste0(pathway_label, " (", pathway_source, ")")
  ) |>
  dplyr::filter(n_members >= MIN_SET, n_members <= MAX_SET)

pathway_source_lookup <- pathway_sets |>
  dplyr::select(pathway_label_src, pathway_source)

pathway_list <- setNames(pathway_sets$members, pathway_sets$pathway_label_src)

cat(sprintf("\nPathway sets for analysis: %d | median size %.0f | range %d-%d\n\n",
            length(pathway_list), median(lengths(pathway_list)),
            min(lengths(pathway_list)), max(lengths(pathway_list))))

# ==============================================================================
# 5. Analysis functions
# ==============================================================================

# ---- ORA (Fisher's exact, one-sided greater) --------------------------------
run_ora <- function(foreground, background, pathway_list,
                    pathway_source_lookup, analysis_label) {
  cat("\n===", analysis_label, "===\n")
  cat("Foreground:", length(foreground), "| Background:", length(background), "\n")
  if (length(foreground) < 1) { cat("Empty foreground — skipping.\n"); return(NULL) }

  bg_only <- background[!background %in% foreground]

  results <- purrr::map_dfr(names(pathway_list), function(pw) {
    members <- pathway_list[[pw]]
    a <- sum(foreground %in% members)
    b <- length(foreground) - a
    c <- sum(bg_only %in% members)
    d <- length(bg_only) - c
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    tibble(
      pathway    = pw,
      n_hits     = a,
      set_size   = a + c,
      odds_ratio = as.numeric(ft$estimate),
      pval       = ft$p.value,
      hit_names  = paste(foreground[foreground %in% members], collapse = ";")
    )
  }) |>
    dplyr::filter(n_hits > 0) |>
    dplyr::mutate(padj = p.adjust(pval, method = "BH")) |>
    dplyr::left_join(pathway_source_lookup, by = c("pathway" = "pathway_label_src")) |>
    dplyr::arrange(padj, pval) |>
    dplyr::mutate(analysis = analysis_label)

  cat("Pathways with >=1 hit:", nrow(results),
      "| p < 0.05:", sum(results$pval < 0.05, na.rm = TRUE),
      "| FDR < 0.05:", sum(results$padj < 0.05, na.rm = TRUE), "\n")
  results
}

# ---- Signed MSEA (fgsea on the raw weight vector) ---------------------------
# scoreType = "std" because the statistic is signed and both tails are
# meaningful: NES > 0 = set concentrated among positively-weighted metabolites.
run_msea_signed <- function(weight_df, pathway_list, analysis_label,
                            nperm = NPERM, seed = 42) {
  cat("\n===", analysis_label, "===\n")

  ranked <- weight_df |>
    dplyr::arrange(dplyr::desc(weight)) |>
    dplyr::select(Name, weight) |>
    deframe()

  cat("Metabolites in ranked list:", length(ranked), "\n")

  pl <- lapply(pathway_list, function(x) intersect(x, names(ranked)))
  pl <- pl[lengths(pl) >= MIN_SET]
  cat("Pathway sets after intersection:", length(pl), "\n")
  if (length(pl) == 0) { cat("No sets with sufficient coverage — skipping.\n"); return(NULL) }

  set.seed(seed)
  results <- fgsea(pathways = pl, stats = ranked, nPermSimple = nperm,
                   minSize = MIN_SET, maxSize = MAX_SET, scoreType = "std") |>
    as_tibble() |>
    dplyr::arrange(padj, pval) |>
    dplyr::left_join(pathway_source_lookup, by = c("pathway" = "pathway_label_src")) |>
    dplyr::mutate(analysis = analysis_label)

  cat("Pathways tested:", nrow(results),
      "| p < 0.05:", sum(results$pval < 0.05, na.rm = TRUE),
      "| FDR < 0.05:", sum(results$padj < 0.05, na.rm = TRUE), "\n")
  results
}

# ---- Plots -------------------------------------------------------------------
strip_src <- function(x) {
  str_remove(x, " \\(KEGG\\)$| \\(Reactome\\)$| \\(WikiPathways\\)$| \\(PFOCR\\)$")
}

plot_ora_dotplot <- function(results, plot_title, top_n = 20) {
  if (is.null(results) || nrow(results) == 0) return(NULL)
  pd <- results |>
    dplyr::filter(pval < 0.05) |>
    head(top_n) |>
    dplyr::mutate(
      stars = dplyr::case_when(padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""),
      pathway_base = strip_src(pathway),
      pathway_base = fct_reorder(pathway_base, -log10(pval))
    )
  if (nrow(pd) == 0) return(NULL)

  ggplot(pd, aes(x = -log10(pval), y = pathway_base)) +
    geom_point(aes(size = n_hits, colour = pathway_source)) +
    geom_text(aes(label = stars), hjust = -1.2, size = 4) +
    scale_size_continuous(range = c(2, 7)) +
    labs(title = plot_title, x = expression(-log[10](italic(p))), y = NULL,
         size = "Hits", colour = "Source",
         caption = "* FDR < 0.05   ** FDR < 0.01") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
}

plot_msea_bar <- function(results, plot_title, top_n = 20) {
  if (is.null(results) || nrow(results) == 0) return(NULL)
  pd <- results |>
    dplyr::filter(pval < 0.05) |>
    head(top_n) |>
    dplyr::mutate(
      stars = dplyr::case_when(padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""),
      pathway_base = fct_reorder(strip_src(pathway), NES),
      direction = ifelse(NES > 0, "Positive weights", "Negative weights")
    )
  if (nrow(pd) == 0) return(NULL)

  ggplot(pd, aes(x = NES, y = pathway_base, fill = direction)) +
    geom_col() +
    geom_text(aes(label = stars, hjust = ifelse(NES > 0, -0.3, 1.3)), size = 4) +
    scale_fill_manual(values = c("Positive weights" = "#B2182B",
                                 "Negative weights" = "#2166AC")) +
    labs(title = plot_title, x = "NES", y = NULL, fill = NULL,
         caption = "* FDR < 0.05   ** FDR < 0.01") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
}

# ==============================================================================
# 6. Signed MSEA — no cutoff, uses the full weight ranking
# ==============================================================================

msea_res <- run_msea_signed(
  w_fac, pathway_list,
  analysis_label = sprintf("%s signed MSEA (metabolite weights)", factor_name)
)

if (!is.null(msea_res)) {
  # leadingEdge is a list-column — flatten before writing
  write_csv(
    msea_res |> dplyr::mutate(leadingEdge = purrr::map_chr(leadingEdge, paste, collapse = ";")),
    file.path(out_dir, sprintf("%s_msea_signed.csv", tolower(factor_name)))
  )
  p <- plot_msea_bar(msea_res, sprintf("%s — signed MSEA on metabolite weights", factor_name))
  if (!is.null(p)) {
    save_dual_format(p, out_dir, sprintf("%s_msea_signed", tolower(factor_name)),
                     width = 10, height = 7)
  }
}

# ==============================================================================
# 7. ORA across top-N cutoffs x {all, positive, negative}
#
# The +/- sets are a partition of the same top-N set (not two independent
# top-N lists), so the three analyses are directly comparable and the union of
# the signed foregrounds is exactly the "all" foreground.
# ==============================================================================

ora_all <- list()

for (N in TOP_N_SET) {
  top_set <- w_fac |> dplyr::arrange(dplyr::desc(abs(weight))) |> head(N)

  fg_all <- top_set$Name
  fg_pos <- top_set |> dplyr::filter(weight > 0) |> dplyr::pull(Name)
  fg_neg <- top_set |> dplyr::filter(weight < 0) |> dplyr::pull(Name)

  cat(sprintf("\n---------- top %d by |weight|: %d positive / %d negative ----------\n",
              N, length(fg_pos), length(fg_neg)))
  cat(sprintf("  |weight| cutoff at rank %d: %.4f\n", N, min(abs(top_set$weight))))

  sub_dir <- file.path(out_dir, paste0("top", N))
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

  write_csv(top_set, file.path(sub_dir, sprintf("top%d_weights.csv", N)))

  res_list <- list(
    all      = run_ora(fg_all, background_all, pathway_list, pathway_source_lookup,
                       sprintf("%s top %d |weight| (all)", factor_name, N)),
    positive = run_ora(fg_pos, background_all, pathway_list, pathway_source_lookup,
                       sprintf("%s top %d |weight| (positive)", factor_name, N)),
    negative = run_ora(fg_neg, background_all, pathway_list, pathway_source_lookup,
                       sprintf("%s top %d |weight| (negative)", factor_name, N))
  )

  for (nm in names(res_list)) {
    res <- res_list[[nm]]
    if (is.null(res) || nrow(res) == 0) next

    res <- res |> dplyr::mutate(direction = nm, top_n = N)
    ora_all[[paste0(nm, "_", N)]] <- res

    write_csv(res, file.path(sub_dir, sprintf("ora_top%d_%s.csv", N, nm)))

    p <- plot_ora_dotplot(res, sprintf("%s — top %d |weight| (%s)", factor_name, N, nm))
    if (!is.null(p)) {
      save_dual_format(p, sub_dir, sprintf("ora_top%d_%s", N, nm),
                       width = 10, height = 7)
    }
  }
}

# ---- combined table ----------------------------------------------------------
ora_combined <- bind_rows(ora_all)
write_csv(ora_combined, file.path(out_dir, sprintf("%s_ora_all_cutoffs.csv",
                                                   tolower(factor_name))))

# ---- side-by-side panel at the primary cutoff -------------------------------
prim <- ora_combined |> dplyr::filter(top_n == TOP_N_PRIMARY)
if (nrow(prim) > 0) {
  panels <- purrr::compact(purrr::map(c("all", "positive", "negative"), function(d) {
    plot_ora_dotplot(prim |> dplyr::filter(direction == d),
                     sprintf("%s (top %d)", d, TOP_N_PRIMARY))
  }))
  if (length(panels) > 0) {
    p_panel <- plot_grid(plotlist = panels, ncol = length(panels))
    title_gg <- cowplot::ggdraw() +
      cowplot::draw_label(sprintf("%s metabolite ORA — all vs. positive vs. negative weights",
                                  factor_name), fontface = "bold", size = 13)
    p_final <- plot_grid(title_gg, p_panel, ncol = 1, rel_heights = c(0.06, 1))
    save_dual_format(p_final, out_dir,
                     sprintf("%s_ora_panel_top%d", tolower(factor_name), TOP_N_PRIMARY),
                     width = 6 * length(panels), height = 7)
  }
}

# ---- cutoff sensitivity: how stable is each pathway across N? ---------------
if (nrow(ora_combined) > 0) {
  stability <- ora_combined |>
    dplyr::filter(pval < 0.05) |>
    dplyr::group_by(direction, pathway) |>
    dplyr::summarise(n_cutoffs = n_distinct(top_n),
                     min_padj  = min(padj, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(direction, dplyr::desc(n_cutoffs), min_padj)
  write_csv(stability, file.path(out_dir, sprintf("%s_ora_cutoff_stability.csv",
                                                  tolower(factor_name))))
  cat("\nPathways nominally significant at all", length(TOP_N_SET), "cutoffs:\n")
  print(as.data.frame(stability |> dplyr::filter(n_cutoffs == length(TOP_N_SET))),
        row.names = FALSE)
}

cat(sprintf("\nDone. Outputs in %s\n", out_dir))
sessionInfo()

# ---- AI assistance disclosure ------------------------------------------------
# Code in this script was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.

# ---- session info ------------------------------------------------------------
cat("\n\n---- Session Info ----\n")
print(sessionInfo())
