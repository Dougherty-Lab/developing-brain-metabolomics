# missingness_visualization_functions.R
# ---------------------------------------------------------------------------
# Missingness visualization helpers: histograms of per-metabolite sample
# missingness, overall and stratified by sex.
#
# Sourced by: Clean_Visualization.Rmd,
#             Sex_and_GW_Analysis.Rmd,
#             Neurotransmitter_Analysis.Rmd,
#             Steroid_Analysis.Rmd
# ---------------------------------------------------------------------------

# ============================================================================
# Missingness Visualization Functions
# ============================================================================

# Plot: # of metabolites (y) vs # of samples missing (x) — combined
plot_missingness_by_sample_pct <- function(missing_stats, title_prefix = "",
                                           png_path = NULL, svg_path = NULL,
                                           width = 10, height = 8, dpi = 300) {
  plot_data <- missing_stats
  n_total <- max(plot_data$total_samples)

  p <- ggplot(plot_data, aes(x = missing_count)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black", boundary = -0.5) +
    theme_minimal() +
    labs(
      title = paste0(title_prefix, "Metabolite Missingness Across Samples"),
      x = "Number of Samples Missing Metabolite",
      y = "Number of Metabolites"
    ) +
    scale_x_continuous(limits = c(-0.5, n_total + 0.5), breaks = seq(0, n_total, 1))

  save_dual_format(
    p, paste0(title_prefix, "missingness_by_sample_pct"),
    png_path, svg_path, width, height, dpi
  )
  return(p)
}

# Plot: # of metabolites (y) vs # of samples missing (x) — sex-split
plot_missingness_by_sample_pct_by_sex <- function(missing_stats_by_sex, title_prefix = "",
                                                  png_path = NULL, svg_path = NULL,
                                                  width = 10, height = 10, dpi = 300) {
  plot_data <- missing_stats_by_sex
  n_total <- max(plot_data$total_samples)

  p <- ggplot(plot_data, aes(x = missing_count, fill = Sex)) +
    geom_histogram(
      binwidth = 1, color = "black", position = "identity",
      alpha = 0.6, boundary = -0.5
    ) +
    theme_minimal() +
    labs(
      title = paste0(title_prefix, "Metabolite Missingness Across Samples by Sex"),
      x = "Number of Samples Missing Metabolite",
      y = "Number of Metabolites"
    ) +
    scale_x_continuous(limits = c(-0.5, n_total + 0.5), breaks = seq(0, n_total, 1)) +
    facet_wrap(~Sex, ncol = 1)

  save_dual_format(
    p, paste0(title_prefix, "missingness_by_sample_pct_by_sex"),
    png_path, svg_path, width, height, dpi
  )
  return(p)
}

# Wrapper function to generate all missingness plots
plot_all_missingness <- function(missing_stats, missing_stats_by_sex = NULL,
                                 sample_metadata = NULL,
                                 title_prefix = "",
                                 png_path = NULL, svg_path = NULL,
                                 width = 10, height = 8, dpi = 300) {
  plots <- list(
    missingness_by_sample_pct = plot_missingness_by_sample_pct(
      missing_stats, title_prefix, png_path, svg_path, width, height, dpi
    )
  )

  # Add sex-stratified plot if data provided
  if (!is.null(missing_stats_by_sex)) {
    plots$missingness_by_sample_pct_by_sex <- plot_missingness_by_sample_pct_by_sex(
      missing_stats_by_sex, title_prefix, png_path, svg_path, width, height + 2, dpi
    )
  }

  return(plots)
}

# ---- AI assistance disclosure ------------------------------------------------
# Code in this file was developed with assistance from Claude (Anthropic).
# All AI-generated code was reviewed, validated, and adapted by the author.
