# ============================================================
# Single-Metabolite Profile Script
# ============================================================
# Produces a 4-panel figure for one metabolite:
#   Panel 1 – Sex comparison (violin + box + jitter)
#   Panel 2 – GW trajectory (overall, no sex coloring)
#   Panel 3 – Sex × GW (plot_trajectory: per-sex lm)
#   Panel 4 – Significance summary across all analyses
#
# Run order: must follow Filtering → Exogenous Flagging → all
#            analysis scripts so CSVs exist on disk.
# ============================================================
setwd("/scratch/jdlab/sneha/developing-brain-metabolomics/src/untargeted")
# ── USER SETTINGS ────────────────────────────────────────────
metabolite_query <- "Prolinamide"   # Full Name OR Compound.ID
dataset          <- "batch2"          # "batch1" | "batch2" | "newbatch2"
results_base     <- "../../results/untargeted"
output_dir       <- "."               # Directory for saved PNG
# ─────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biobase)
  library(cowplot)
  library(ggtext)
})
source("functions.R")

# ============================================================
# 1. Dataset configuration
# ============================================================
dataset_config <- list(
  batch1 = list(
    peak_file       = "batch1_peak_area_values.csv",
    metadata_file   = "batch1_sample_metadata.csv",
    annotation_file = "batch1_metabolite_annotations.csv",
    csv_folder      = "batch1_filtering/csv",
    subgroup_prefix = "batch1",
    has_analysis    = FALSE,
    has_pathway     = FALSE
  ),
  batch2 = list(
    peak_file       = "batch2_peak_area_clean.csv",
    metadata_file   = "batch2_sample_metadata.csv",
    annotation_file = "batch2_metabolite_annotations.csv",
    csv_folder      = "batch2_clean/csv",
    subgroup_prefix = "batch2",
    has_analysis    = TRUE,
    has_pathway     = TRUE   # MSEA / ORA only run for batch2
  ),
  newbatch2 = list(
    peak_file       = "newbatch2_peak_area_clean.csv",
    metadata_file   = "newbatch2_sample_metadata.csv",
    annotation_file = "newbatch2_metabolite_annotations.csv",
    csv_folder      = "newbatch2_clean/csv",
    subgroup_prefix = "newbatch2",
    has_analysis    = TRUE,
    has_pathway     = FALSE
  )
)

if (!dataset %in% names(dataset_config))
  stop("dataset must be one of: 'batch1', 'batch2', 'newbatch2'")

cfg     <- dataset_config[[dataset]]
csv_dir <- file.path(results_base, cfg$csv_folder)
pfx     <- cfg$subgroup_prefix

# ============================================================
# 2. Load data
# ============================================================
cat("Loading data for dataset:", dataset, "\n")

peak_area <- read_csv(file.path(results_base, cfg$peak_file),
                      show_col_types = FALSE)
metadata  <- read_csv(file.path(results_base, cfg$metadata_file),
                      show_col_types = FALSE)
annot     <- read_csv(file.path(results_base, cfg$annotation_file),
                      show_col_types = FALSE)

# ============================================================
# 3. Resolve metabolite (Name or Compound.ID, partial OK)
# ============================================================
resolve_metabolite <- function(query, annot_df) {
  # Exact Name match
  if (query %in% annot_df$Name) {
    idx <- which(annot_df$Name == query)[1]
    return(list(name = annot_df$Name[idx], cid = annot_df$Compound.ID[idx]))
  }
  # Exact Compound.ID match
  if (query %in% annot_df$Compound.ID) {
    idx <- which(annot_df$Compound.ID == query)[1]
    return(list(name = annot_df$Name[idx], cid = annot_df$Compound.ID[idx]))
  }
  # Partial Name match (case-insensitive)
  idx <- grep(query, annot_df$Name, ignore.case = TRUE)
  if (length(idx) == 0) stop("Metabolite not found: ", query)
  if (length(idx) > 1) {
    cat("Multiple partial matches — using first. All matches:\n")
    cat(paste(annot_df$Name[idx], collapse = "\n"), "\n\n")
  }
  list(name = annot_df$Name[idx[1]], cid = annot_df$Compound.ID[idx[1]])
}

hit      <- resolve_metabolite(metabolite_query, annot)
met_name <- hit$name
met_cid  <- hit$cid

cat("Metabolite : ", met_name, "\n")
cat("Compound ID:", met_cid,  "\n\n")

# ============================================================
# 4. Build ExpressionSet  (mirrors main analysis scripts)
# ============================================================
compounds_in_peak <- unique(peak_area$Compound.ID)

annot_filtered <- annot %>%
  filter(Compound.ID %in% compounds_in_peak)

final_matrix <- peak_area %>%
  dplyr::select(starts_with("Sample"))

metabolite_names <- annot_filtered$Name
final_matrix     <- apply(final_matrix, 2, as.numeric)
rownames(final_matrix) <- metabolite_names

metadata <- metadata %>%
  filter(Sample %in% colnames(final_matrix)) %>%
  arrange(match(Sample, colnames(final_matrix)))

sample_meta_df          <- as.data.frame(metadata)
rownames(sample_meta_df) <- sample_meta_df$Sample

met_meta_df          <- as.data.frame(annot_filtered)
rownames(met_meta_df) <- met_meta_df$Name

eset <- ExpressionSet(
  assayData  = final_matrix,
  phenoData  = AnnotatedDataFrame(sample_meta_df),
  featureData = AnnotatedDataFrame(met_meta_df)
)

# Log2 transform (consistent with analysis scripts)
exprs(eset) <- log2(exprs(eset) + 1)

# Check metabolite is present in the filtered set
if (!met_name %in% rownames(eset)) {
  stop("'", met_name, "' was not retained after filtering. ",
       "It may have been removed during the filtering step.")
}

# ============================================================
# 5. Extract single-metabolite data
# ============================================================
intensity  <- exprs(eset)[met_name, ]
pheno      <- pData(eset)

plot_df <- data.frame(
  Sample    = pheno$Sample,
  GW        = pheno$GW,
  Sex       = pheno$Sex,
  Intensity = intensity
)

sex_colors <- c("F" = "#FFC20A", "M" = "#571C90")

# ============================================================
# 6. Panel 1 – Sex comparison
# ============================================================
p1 <- ggplot(plot_df, aes(x = Sex, y = Intensity, fill = Sex, color = Sex)) +
  geom_violin(alpha = 0.35, trim = FALSE, linewidth = 0.5) +
  geom_boxplot(width = 0.18, alpha = 0.7, outlier.shape = NA, linewidth = 0.5,
               color = "grey30") +
  geom_jitter(width = 0.08, size = 1.8, alpha = 0.7) +
  scale_fill_manual(values  = sex_colors) +
  scale_color_manual(values = sex_colors) +
  labs(
    x = "Sex", y = "log2 Peak Area",
    title = "Sex Comparison"
  ) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5)
  )

# ============================================================
# 7. Panel 2 – GW trajectory (overall, no sex split)
# ============================================================
p2 <- ggplot(plot_df, aes(x = GW, y = Intensity)) +
  geom_point(size = 2.5, alpha = 0.75, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "#1f78b4",
              linewidth = 0.8, fill = "#a6cee3", alpha = 0.25) +
  labs(
    x = "Gestational Week", y = "log2 Peak Area",
    title = "GW Trajectory"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# ============================================================
# 8. Panel 3 – Sex × GW  (plot_trajectory from functions.R)
# ============================================================
annot_confidence <- annot %>%
  filter(Name == met_name) %>%
  pull(Best.Match) %>%
  .[1]

p3 <- plot_trajectory(
  metabolite_name      = met_name,
  eset                 = eset,
  title_suffix         = "Sex × GW",
  annotation_confidence = annot_confidence
)

# ============================================================
# 9. Panel 4 – Significance summary
# ============================================================

# ── 9a. Helper: safe CSV read ───────────────────────────────
safe_read <- function(path) {
  if (file.exists(path)) {
    read_csv(path, show_col_types = FALSE)
  } else {
    NULL
  }
}

# ── 9b. Lookup helper ───────────────────────────────────────
# Returns a one-row tibble: (test, category, status, value_label, direction)
# status: "FDR" | "Nominal" | "NS" | "VIP>2" | "VIP>1" | "In pathway" | "NA"
limma_row <- function(label, category, df, fdr_col = "adj.P.Val",
                      p_col = "P.Value", fc_col = "logFC") {
  if (is.null(df)) return(tibble(test=label, category=category,
                                  status="NA", value="–", direction="–"))
  row <- df %>% filter(Name == met_name | Compound.ID == met_cid)
  if (nrow(row) == 0) return(tibble(test=label, category=category,
                                     status="NA", value="–", direction="–"))
  fdr <- row[[fdr_col]][1]
  pv  <- row[[p_col]][1]
  fc  <- if (fc_col %in% names(row)) row[[fc_col]][1] else NA_real_

  status <- if (!is.na(fdr) && fdr < 0.1) "FDR" else
            if (!is.na(pv)  && pv  < 0.05) "Nominal" else "NS"

  dir_label <- if (!is.na(fc)) {
    if (fc > 0) "↑" else "↓"
  } else "–"

  value_label <- if (!is.na(fdr)) sprintf("%.3f", fdr) else "–"

  tibble(test=label, category=category, status=status,
         value=value_label, direction=dir_label)
}

wilcox_row <- function(label, category, df,
                       fdr_col = "FDR", p_col = "p_value", fc_col = "logFC") {
  limma_row(label, category, df, fdr_col, p_col, fc_col)
}

vip_row <- function(label, category, df, vip_col = "VIP_comp1") {
  if (is.null(df)) return(tibble(test=label, category=category,
                                  status="NA", value="–", direction="–"))
  row <- df %>% filter(Name == met_name | Compound.ID == met_cid)
  if (nrow(row) == 0) return(tibble(test=label, category=category,
                                     status="NA", value="–", direction="–"))
  vip <- row[[vip_col]][1]
  status <- if (!is.na(vip) && vip > 2) "VIP>2" else
            if (!is.na(vip) && vip > 1) "VIP>1" else "NS"
  tibble(test=label, category=category, status=status,
         value=sprintf("%.2f", vip), direction="–")
}

msea_row <- function(label, category, df) {
  if (is.null(df)) return(tibble(test=label, category=category,
                                  status="NA", value="–", direction="–",
                                  pathway_hits=list(character(0))))
  # Check if metabolite Name appears in leadingEdge of any FDR < 0.1 pathway
  sig_paths <- df %>%
    filter(padj < 0.1) %>%
    filter(map_lgl(strsplit(leadingEdge, ";"),
                   ~ met_name %in% trimws(.x)))
  status <- if (nrow(sig_paths) > 0) "In pathway" else "NS"
  pw_hits <- sig_paths$pathway
  tibble(test=label, category=category, status=status,
         value=as.character(nrow(sig_paths)), direction="–",
         pathway_hits=list(pw_hits))
}

ora_row <- function(label, category, df) {
  if (is.null(df)) return(tibble(test=label, category=category,
                                  status="NA", value="–", direction="–",
                                  pathway_hits=list(character(0))))
  # Metabolite appears in hit_names of a FDR-significant pathway
  sig_paths <- df %>%
    filter(padj < 0.1) %>%
    filter(map_lgl(strsplit(hit_names, ";"),
                   ~ met_name %in% trimws(.x)))
  status <- if (nrow(sig_paths) > 0) "In pathway" else "NS"
  pw_hits <- sig_paths$pathway
  tibble(test=label, category=category, status=status,
         value=as.character(nrow(sig_paths)), direction="–",
         pathway_hits=list(pw_hits))
}

# ── 9c. Build summary rows ──────────────────────────────────

if (!cfg$has_analysis) {

  summary_rows <- tibble(
    test      = "No analysis available for Batch 1",
    category  = "Main",
    status    = "NA",
    value     = "–",
    direction = "–"
  )
  pathway_names <- character(0)

} else {

  # Main analysis CSVs
  gw_limma    <- safe_read(file.path(csv_dir, paste0(pfx, "_gw_limma_results.csv")))
  sex_limma   <- safe_read(file.path(csv_dir, paste0(pfx, "_sex_limma_results.csv")))
  sex_wilcox  <- safe_read(file.path(csv_dir, paste0(pfx, "_sex_differential_analysis.csv")))
  sex_vip     <- safe_read(file.path(csv_dir, paste0(pfx, "_sex_plsda_vip.csv")))
  gw_vip      <- safe_read(file.path(csv_dir, paste0(pfx, "_gw_plsda_vip.csv")))
  sex_adj     <- safe_read(file.path(csv_dir, paste0(pfx, "_sex_adjusted_limma_results.csv")))
  gw_adj      <- safe_read(file.path(csv_dir, paste0(pfx, "_gw_adjusted_limma_results.csv")))
  int_limma   <- safe_read(file.path(csv_dir, paste0(pfx, "_interaction_limma_results.csv")))

  # Subgroup CSVs – try "all" scope first, fall back to first available
  find_subgroup_csv <- function(subtype, suffix) {
    base_dir  <- file.path(results_base, paste0(pfx, "_", subtype))
    all_path  <- file.path(base_dir, "all", "csv",
                           paste0(pfx, "_", subtype, "_all_", suffix, ".csv"))
    if (file.exists(all_path)) return(safe_read(all_path))
    # Fall back: look for any scope CSV matching the suffix
    candidates <- list.files(base_dir, recursive = TRUE,
                             pattern = paste0("_", suffix, "\\.csv$"),
                             full.names = TRUE)
    if (length(candidates) > 0) return(safe_read(candidates[1]))
    NULL
  }

  steroid_gw     <- find_subgroup_csv("steroids",          "gw_limma")
  nt_gw          <- find_subgroup_csv("neurotransmitters",  "gw_limma")
  nt_sex_limma   <- find_subgroup_csv("neurotransmitters",  "sex_limma")
  nt_sex_wilcox  <- find_subgroup_csv("neurotransmitters",  "sex_wilcoxon")

  # Pathway CSVs (batch2 only)
  if (cfg$has_pathway) {
    msea_sex <- safe_read(file.path(csv_dir, "batch2_msea_sex.csv"))
    msea_gw  <- safe_read(file.path(csv_dir, "batch2_msea_gw.csv"))
    ora_sex  <- safe_read(file.path(csv_dir, "batch2_ora_sex_limma_nominal.csv"))
    ora_gw_n <- safe_read(file.path(csv_dir, "batch2_ora_gw_limma_nominal.csv"))
    ora_gw_f <- safe_read(file.path(csv_dir, "batch2_ora_gw_limma_fdr.csv"))
    ora_vip1 <- safe_read(file.path(csv_dir, "batch2_ora_sex_plsda_vip1.csv"))
    ora_vip2 <- safe_read(file.path(csv_dir, "batch2_ora_sex_plsda_vip2.csv"))
  } else {
    msea_sex <- msea_gw <- ora_sex <- ora_gw_n <- ora_gw_f <-
      ora_vip1 <- ora_vip2 <- NULL
  }

  # Collect rows
  summary_rows <- bind_rows(
    # ── Main ──────────────────────────────────────────────────
    limma_row("GW limma",            "Main", gw_limma),
    limma_row("Sex limma",           "Main", sex_limma),
    wilcox_row("Sex Wilcoxon",       "Main", sex_wilcox),
    vip_row("Sex PLS-DA VIP",        "Main", sex_vip),
    vip_row("GW PLS-DA VIP",         "Main", gw_vip),
    limma_row("Sex limma (adj. GW)", "Main", sex_adj),
    limma_row("GW limma (adj. Sex)", "Main", gw_adj),
    limma_row("Sex × GW interaction","Main", int_limma),
    # ── Subgroup ──────────────────────────────────────────────
    limma_row("Steroid GW limma",    "Subgroup", steroid_gw),
    limma_row("NT GW limma",         "Subgroup", nt_gw),
    limma_row("NT Sex limma",        "Subgroup", nt_sex_limma),
    wilcox_row("NT Sex Wilcoxon",    "Subgroup", nt_sex_wilcox,
               fdr_col="FDR", p_col="p_value", fc_col="logFC"),
    # ── Pathway ───────────────────────────────────────────────
    msea_row("MSEA Sex",             "Pathway", msea_sex),
    msea_row("MSEA GW",              "Pathway", msea_gw),
    ora_row("ORA Sex nominal",       "Pathway", ora_sex),
    ora_row("ORA GW nominal",        "Pathway", ora_gw_n),
    ora_row("ORA GW FDR",            "Pathway", ora_gw_f),
    ora_row("ORA VIP > 1",           "Pathway", ora_vip1),
    ora_row("ORA VIP > 2",           "Pathway", ora_vip2)
  )

  # Collect significant pathway names for annotation text
  pathway_names <- character(0)
  if ("pathway_hits" %in% names(summary_rows)) {
    pathway_names <- summary_rows %>%
      filter(status == "In pathway") %>%
      pull(pathway_hits) %>%
      unlist() %>%
      unique()
  }
}

# ── 9d. Status factor & colors ──────────────────────────────
status_levels <- c("FDR", "Nominal", "NS", "VIP>2", "VIP>1", "In pathway", "NA")
status_colors <- c(
  "FDR"        = "#2ca02c",
  "Nominal"    = "#ff7f0e",
  "NS"         = "grey82",
  "VIP>2"      = "#1f77b4",
  "VIP>1"      = "#aec7e8",
  "In pathway" = "#9467bd",
  "NA"         = "white"
)

summary_rows <- summary_rows %>%
  mutate(
    status   = factor(status, levels = status_levels),
    category = factor(category, levels = c("Main", "Subgroup", "Pathway")),
    test     = factor(test, levels = rev(unique(test)))   # bottom-to-top
  )

# ── 9e. Build tile plot ─────────────────────────────────────
p4 <- ggplot(summary_rows,
             aes(x = category, y = test, fill = status)) +
  geom_tile(color = "white", linewidth = 0.6, width = 0.85, height = 0.85) +
  geom_text(aes(label = ifelse(status %in% c("NA", "NS"), "", value)),
            size = 3, color = "white", fontface = "bold") +
  geom_text(aes(label = ifelse(direction %in% c("↑","↓"), direction, "")),
            nudge_x = 0.28, size = 3.2, color = "grey20") +
  scale_fill_manual(values  = status_colors,
                    name    = "Result",
                    drop    = FALSE) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "Analysis Summary") +
  theme(
    axis.text.x       = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 9),
    panel.grid        = element_blank(),
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    legend.position   = "right",
    legend.key.size   = unit(0.55, "cm"),
    legend.text       = element_text(size = 8)
  )

# Pathway annotation text (if any)
if (length(pathway_names) > 0) {
  pw_text <- paste("Significant pathways:", paste(pathway_names, collapse = " | "))
  p4 <- p4 +
    labs(caption = str_wrap(pw_text, width = 90)) +
    theme(plot.caption = element_text(size = 7.5, hjust = 0, color = "#9467bd"))
}

# ============================================================
# 10. Assemble & save
# ============================================================
top_row <- plot_grid(p1, p2, p3, nrow = 1, labels = c("A","B","C"),
                     label_size = 11)

full_fig <- plot_grid(
  top_row, p4,
  nrow        = 2,
  rel_heights = c(1, 1.1),
  labels      = c("", "D"),
  label_size  = 11
)

title_grob <- ggdraw() +
  draw_label(met_name, fontface = "bold", size = 13, hjust = 0.5) +
  draw_label(paste0("Dataset: ", dataset, "  |  Compound ID: ", met_cid),
             size = 8.5, color = "grey40", hjust = 0.5, y = 0.2)

final_plot <- plot_grid(title_grob, full_fig,
                        ncol = 1, rel_heights = c(0.06, 1))

out_file <- file.path(output_dir,
                      paste0("profile_",
                             gsub("[^A-Za-z0-9_]", "_", met_name),
                             "_", dataset, ".png"))

ggsave(out_file, final_plot, width = 15, height = 12, dpi = 180)
cat("\nSaved:", out_file, "\n")
