# ============================================================================
# Missingness Visualization Functions
# ============================================================================

plot_missingness_histogram <- function(missing_stats, title_prefix = "", 
                                       png_path = NULL, svg_path = NULL,
                                       width = 10, height = 8, dpi = 300) {
  p <- ggplot(missing_stats, aes(x = missing_proportion)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Metabolite Missingness"),
         x = "Proportion Missing",
         y = "Number of Metabolites") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  
  save_dual_format(p, paste0(title_prefix, "missingness_histogram"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_missingness_histogram_by_sex <- function(missing_stats_by_sex, title_prefix = "", 
                                              png_path = NULL, svg_path = NULL,
                                              width = 10, height = 10, dpi = 300) {
  p <- ggplot(missing_stats_by_sex, aes(x = missing_proportion, fill = Sex)) +
    geom_histogram(binwidth = 0.05, color = "black", position = "identity", alpha = 0.6) +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Metabolite Missingness by Sex"),
         x = "Proportion Missing",
         y = "Number of Metabolites") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    facet_wrap(~Sex, ncol = 1)
  
  save_dual_format(p, paste0(title_prefix, "missingness_histogram_by_sex"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_missingness_count_distribution <- function(missing_stats, title_prefix = "", 
                                                png_path = NULL, svg_path = NULL,
                                                width = 10, height = 8, dpi = 300) {
  # Create summary of missingness counts
  missingness_count_summary <- missing_stats %>%
    group_by(missing_count) %>%
    summarise(n_metabolites = n(), .groups = "drop")
  
  p <- ggplot(missingness_count_summary, aes(x = missing_count, y = n_metabolites)) +
    geom_bar(stat = "identity", fill = "darkorange", color = "black") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Missing Sample Counts"),
         x = "Number of Missing Samples",
         y = "Number of Metabolites")
  
  save_dual_format(p, paste0(title_prefix, "missingness_count_distribution"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_missingness_count_distribution_by_sex <- function(missing_stats_by_sex, title_prefix = "", 
                                                       png_path = NULL, svg_path = NULL,
                                                       width = 12, height = 8, dpi = 300) {
  # Create summary of missingness counts by sex
  missingness_count_summary <- missing_stats_by_sex %>%
    group_by(missing_count, Sex) %>%
    summarise(n_metabolites = n(), .groups = "drop")
  
  p <- ggplot(missingness_count_summary, aes(x = missing_count, y = n_metabolites, fill = Sex)) +
    geom_bar(stat = "identity", color = "black", position = "dodge") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Missing Sample Counts by Sex"),
         x = "Number of Missing Samples",
         y = "Number of Metabolites")
  
  save_dual_format(p, paste0(title_prefix, "missingness_count_distribution_by_sex"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_missingness_vs_gw <- function(missing_stats, sample_metadata, title_prefix = "", 
                                   png_path = NULL, svg_path = NULL,
                                   width = 10, height = 8, dpi = 300) {
  # This function requires additional context - for now, create a summary plot
  # that could incorporate GW if metadata is provided with metabolite-level GW info
  
  p <- ggplot(missing_stats, aes(x = missing_proportion, fill = cut(missing_proportion, 
                                                                      breaks = c(0, 0.25, 0.5, 0.75, 1),
                                                                      labels = c("0-25%", "25-50%", "50-75%", "75-100%")))) +
    geom_histogram(binwidth = 0.05, color = "black") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Missingness Categories"),
         x = "Proportion Missing",
         y = "Number of Metabolites",
         fill = "Missingness Category") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  
  save_dual_format(p, paste0(title_prefix, "missingness_categories"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_sample_completeness <- function(missing_stats, sample_metadata, title_prefix = "", 
                                     png_path = NULL, svg_path = NULL,
                                     width = 12, height = 8, dpi = 300) {
  # Calculate completeness per sample from missingness data
  # Note: This would need to be restructured from per-metabolite to per-sample
  
  p <- ggplot(missing_stats, aes(x = detected_samples)) +
    geom_histogram(binwidth = 5, fill = "purple", color = "black") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Sample Completeness"),
         x = "Number of Metabolites Detected",
         y = "Frequency")
  
  save_dual_format(p, paste0(title_prefix, "sample_completeness"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

# Wrapper function to generate all missingness plots
plot_all_missingness <- function(missing_stats, missing_stats_by_sex = NULL, 
                                 sample_metadata = NULL,
                                 title_prefix = "", 
                                 png_path = NULL, svg_path = NULL,
                                 width = 10, height = 8, dpi = 300) {
  
  plots <- list(
    histogram = plot_missingness_histogram(missing_stats, title_prefix, png_path, svg_path, width, height, dpi),
    count_distribution = plot_missingness_count_distribution(missing_stats, title_prefix, png_path, svg_path, width, height, dpi),
    sample_completeness = plot_sample_completeness(missing_stats, sample_metadata, title_prefix, png_path, svg_path, width + 2, height, dpi)
  )
  
  # Add sex-stratified plots if data provided
  if (!is.null(missing_stats_by_sex)) {
    plots$histogram_by_sex <- plot_missingness_histogram_by_sex(missing_stats_by_sex, title_prefix, png_path, svg_path, width, height + 2, dpi)
    plots$count_distribution_by_sex <- plot_missingness_count_distribution_by_sex(missing_stats_by_sex, title_prefix, png_path, svg_path, width + 2, height, dpi)
  }
  
  return(plots)
}
