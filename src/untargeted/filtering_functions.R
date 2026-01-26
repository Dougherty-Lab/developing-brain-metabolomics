# =============================================================================
# OUTLIER REMOVAL FUNCTIONS
# =============================================================================

#' Identify and Remove Sample Outliers Based on Total Intensity
#'
#' @param data Data frame containing metabolomics peak data
#' @param metadata Data frame containing sample metadata
#' @param threshold Z-score threshold for outlier detection (default = 3)
#' @param sample_cols Columns containing sample intensity data (default = 4:ncol(data))
#' @param save_plot Logical, whether to save diagnostic plots (default = TRUE)
#' @param output_prefix Prefix for output files (default = "intensity_outliers")
#' @param png_path Path to save PNG/SVG plots
#' @param csv_path Path to save CSV files
#'
#' @return List containing filtered data and outlier information
remove_intensity_outliers <- function(data, 
                                      metadata,
                                      threshold = 3, 
                                      sample_cols = 4:ncol(data),
                                      save_plot = TRUE,
                                      output_prefix = "intensity_outliers",
                                      png_path = "../../results/untargeted/Batch 2",
                                      csv_path = "../../results/untargeted/Batch 2") {
  
  # Extract sample columns
  sample_data <- data[, sample_cols]
  
  # Calculate total intensity per sample
  sample_sums <- colSums(sample_data, na.rm = TRUE)
  
  # Calculate z-scores
  z_scores <- scale(sample_sums)
  z_scores_vec <- as.vector(z_scores)
  
  # Identify outliers
  outlier_mask <- abs(z_scores_vec) > threshold
  outlier_samples <- names(sample_data)[outlier_mask]
  
  # Create diagnostic data frame
  diagnostic_df <- data.frame(
    Sample = names(sample_data),
    TotalIntensity = sample_sums,
    ZScore = z_scores_vec,
    IsOutlier = outlier_mask
  )
  
  # Merge with metadata for additional context
  diagnostic_df <- diagnostic_df %>%
    left_join(metadata, by = "Sample")
  
  # Print summary
  cat("\n=== Intensity Outlier Removal Summary ===\n")
  cat("Z-score threshold:", threshold, "\n")
  cat("Total samples:", ncol(sample_data), "\n")
  cat("Outliers detected:", sum(outlier_mask), "\n")
  if(sum(outlier_mask) > 0) {
    cat("\nOutlier samples:\n")
    print(diagnostic_df %>% filter(IsOutlier) %>% 
            select(Sample, TotalIntensity, ZScore))
  }
  cat("Samples retained:", sum(!outlier_mask), "\n\n")
  
  # Create diagnostic plot
  if(save_plot) {
    p <- ggplot(diagnostic_df, aes(x = reorder(Sample, TotalIntensity), 
                                   y = TotalIntensity, 
                                   color = IsOutlier)) +
      geom_point(size = 3) +
      geom_hline(yintercept = mean(sample_sums) + threshold * sd(sample_sums), 
                 linetype = "dashed", color = "red") +
      geom_hline(yintercept = mean(sample_sums) - threshold * sd(sample_sums), 
                 linetype = "dashed", color = "red") +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                         labels = c("Retained", "Outlier")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = "Sample Total Intensity Distribution",
           subtitle = paste0("Z-score threshold: ±", threshold),
           x = "Sample",
           y = "Total Intensity",
           color = "Status")
    
    # Save as SVG
    ggsave(file.path(png_path, paste0(output_prefix, "_diagnostic_plot.svg")), 
           p, width = 12, height = 6)
    
    # Also save as PNG
    ggsave(file.path(png_path, paste0(output_prefix, "_diagnostic_plot.png")), 
           p, width = 12, height = 6, dpi = 300)
  }
  
  # Filter data to remove outliers
  filtered_data <- data %>%
    select(-all_of(outlier_samples))
  
  # Save diagnostic information
  write.csv(diagnostic_df, 
            file.path(csv_path, paste0(output_prefix, "_diagnostic_info.csv")), 
            row.names = FALSE)
  
  # Return results
  return(list(
    filtered_data = filtered_data,
    outlier_samples = outlier_samples,
    diagnostic_df = diagnostic_df,
    n_outliers = sum(outlier_mask),
    n_retained = sum(!outlier_mask)
  ))
}


#' Identify and Remove Sample Outliers Based on Number of Detected Peaks
#'
#' @param data Data frame containing metabolomics peak data
#' @param metadata Data frame containing sample metadata
#' @param threshold Z-score threshold for outlier detection (default = 3)
#' @param sample_cols Columns containing sample intensity data (default = 4:ncol(data))
#' @param na_threshold Value below which to consider as missing/not detected (default = 0)
#' @param save_plot Logical, whether to save diagnostic plots (default = TRUE)
#' @param output_prefix Prefix for output files (default = "detection_outliers")
#' @param png_path Path to save PNG/SVG plots
#' @param csv_path Path to save CSV files
#'
#' @return List containing filtered data and outlier information
remove_detection_outliers <- function(data, 
                                      metadata,
                                      threshold = 3, 
                                      sample_cols = 4:ncol(data),
                                      na_threshold = 0,
                                      save_plot = TRUE,
                                      output_prefix = "detection_outliers",
                                      png_path = "../../results/untargeted/Batch 2",
                                      csv_path = "../../results/untargeted/Batch2") {
  
  # Extract sample columns
  sample_data <- data[, sample_cols]
  
  # Calculate number of detected peaks per sample (non-NA and > na_threshold)
  peaks_detected <- colSums(!is.na(sample_data) & sample_data > na_threshold)
  
  # Calculate z-scores
  z_scores <- scale(peaks_detected)
  z_scores_vec <- as.vector(z_scores)
  
  # Identify outliers
  outlier_mask <- abs(z_scores_vec) > threshold
  outlier_samples <- names(sample_data)[outlier_mask]
  
  # Create diagnostic data frame
  diagnostic_df <- data.frame(
    Sample = names(sample_data),
    PeaksDetected = peaks_detected,
    ZScore = z_scores_vec,
    IsOutlier = outlier_mask,
    DetectionRate = peaks_detected / nrow(sample_data) * 100
  )
  
  # Merge with metadata
  diagnostic_df <- diagnostic_df %>%
    left_join(metadata, by = "Sample")
  
  # Print summary
  cat("\n=== Detection Outlier Removal Summary ===\n")
  cat("Z-score threshold:", threshold, "\n")
  cat("Total samples:", ncol(sample_data), "\n")
  cat("Total peaks:", nrow(sample_data), "\n")
  cat("Outliers detected:", sum(outlier_mask), "\n")
  if(sum(outlier_mask) > 0) {
    cat("\nOutlier samples:\n")
    print(diagnostic_df %>% filter(IsOutlier) %>% 
            select(Sample, PeaksDetected, DetectionRate, ZScore))
  }
  cat("Samples retained:", sum(!outlier_mask), "\n")
  cat("Mean detection rate:", round(mean(diagnostic_df$DetectionRate), 2), "%\n\n")
  
  # Create diagnostic plot
  if(save_plot) {
    p <- ggplot(diagnostic_df, aes(x = reorder(Sample, PeaksDetected), 
                                   y = PeaksDetected, 
                                   color = IsOutlier)) +
      geom_point(size = 3) +
      geom_hline(yintercept = mean(peaks_detected) + threshold * sd(peaks_detected), 
                 linetype = "dashed", color = "red") +
      geom_hline(yintercept = mean(peaks_detected) - threshold * sd(peaks_detected), 
                 linetype = "dashed", color = "red") +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                         labels = c("Retained", "Outlier")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = "Number of Detected Peaks per Sample",
           subtitle = paste0("Z-score threshold: ±", threshold),
           x = "Sample",
           y = "Number of Peaks Detected",
           color = "Status")
    
    # Save as SVG
    ggsave(file.path(png_path, paste0(output_prefix, "_diagnostic_plot.svg")), 
           p, width = 12, height = 6)
    
    # Also save as PNG
    ggsave(file.path(png_path, paste0(output_prefix, "_diagnostic_plot.png")), 
           p, width = 12, height = 6, dpi = 300)
  }
  
  # Filter data to remove outliers
  filtered_data <- data %>%
    select(-all_of(outlier_samples))
  
  # Save diagnostic information
  write.csv(diagnostic_df, 
            file.path(csv_path, paste0(output_prefix, "_diagnostic_info.csv")), 
            row.names = FALSE)
  
  # Return results
  return(list(
    filtered_data = filtered_data,
    outlier_samples = outlier_samples,
    diagnostic_df = diagnostic_df,
    n_outliers = sum(outlier_mask),
    n_retained = sum(!outlier_mask)
  ))
}


# RT DENSITY-BASED THRESHOLD DETECTION FUNCTION
# Automatically identifies RT cutoffs based on metabolite density distributions
# Uses first significant density decrease as threshold


# ============================================================================
# MAIN FUNCTION: Detect RT threshold from density
# ============================================================================

detect_rt_threshold_from_density <- function(
    compound_annotation,
    column_mode_combinations = NULL,  # If NULL, will detect all combinations
    bin_size = 10,                    # in seconds (fixed at 10s - most reliable)
    smoothing_window = 3,              # for smoothing density curve
    min_density_drop_fraction = 0.3,   # require 30% drop from peak to detect threshold
    plot_diagnostics = TRUE,
    output_dir = NULL,                # Legacy parameter - will be overridden if specific paths given
    png_path = NULL,
    svg_path = NULL,
    csv_path = NULL
) {
  
  # Handle legacy output_dir parameter
  if (!is.null(output_dir) && is.null(png_path) && is.null(svg_path) && is.null(csv_path)) {
    png_path <- output_dir
    svg_path <- output_dir
    csv_path <- output_dir
  }
  
  # Set defaults if still null
  if (is.null(png_path)) png_path <- "../../results/untargeted/RT_threshold_diagnostics"
  if (is.null(svg_path)) svg_path <- "../../results/untargeted/RT_threshold_diagnostics"
  if (is.null(csv_path)) csv_path <- "../../results/untargeted/RT_threshold_diagnostics"
  
  # Create output directories if needed
  for (path in c(png_path, svg_path, csv_path)) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Extract Column and Mode from Compound.ID if not provided
  data_with_colmode <- compound_annotation %>%
    mutate(
      Column = sub("-[^-]+-[0-9]+$", "", Compound.ID),
      Mode = sub(".*-(.*?)-.*", "\\1", Compound.ID),
      RT_min = `RT[min]`  # ensure RT column is accessible
    ) %>%
    select(Compound.ID, Column, Mode, RT_min)
  
  # Identify all column-mode combinations if not provided
  if (is.null(column_mode_combinations)) {
    column_mode_combinations <- data_with_colmode %>%
      distinct(Column, Mode) %>%
      arrange(Column, Mode)
  }
  
  # Store results for all combinations and bin sizes
  all_thresholds <- list()
  all_density_data <- list()
  
  # ========================================================================
  # PROCESS EACH COLUMN-MODE COMBINATION
  # ========================================================================
  
  for (i in seq_len(nrow(column_mode_combinations))) {
    col_name <- column_mode_combinations$Column[i]
    mode_name <- column_mode_combinations$Mode[i]
    
    cat("\n--- Processing:", col_name, "-", mode_name, "---\n")
    
    # Filter data for this column-mode
    subset_data <- data_with_colmode %>%
      filter(Column == col_name, Mode == mode_name) %>%
      pull(RT_min) %>%
      na.omit()
    
    n_compounds <- length(subset_data)
    cat("  Total compounds:", n_compounds, "\n")
    
    if (n_compounds < 10) {
      cat("  WARNING: Very few compounds. Results may be unreliable.\n")
    }
    
    # Get RT range
    rt_min <- min(subset_data, na.rm = TRUE)
    rt_max <- max(subset_data, na.rm = TRUE)
    cat("  RT range:", rt_min, "-", rt_max, "minutes\n")
    
    # ====================================================================
    # ANALYZE WITH 10-SECOND BIN SIZE (most reliable)
    # ====================================================================
    
    bin_min <- bin_size / 60  # convert to minutes
    
    cat("  \n  Bin size:", bin_size, "seconds (", bin_min, "min)\n")
    
    # Create bins - ensure they fully span the data range with small margin
    # This prevents histogram errors when data falls outside breaks
    bin_start <- floor(rt_min * 100) / 100  # Floor to 2 decimal places
    bin_end <- ceiling(rt_max * 100) / 100   # Ceiling to 2 decimal places
    
    bin_breaks <- seq(bin_start - bin_min, bin_end + bin_min, by = bin_min)
    bin_centers <- bin_breaks[-length(bin_breaks)] + bin_min / 2
    
    # Calculate density (counts per bin)
    density_raw <- hist(subset_data, breaks = bin_breaks, plot = FALSE)$counts
    
    # Smooth the density curve
    if (length(density_raw) > smoothing_window) {
      density_smooth <- zoo::rollmean(density_raw, k = smoothing_window, fill = NA)
    } else {
      density_smooth <- density_raw
    }
    
    # ===== IMPROVED THRESHOLD DETECTION =====
    # Find the FIRST SIGNIFICANT DROP in the density curve
    # Don't wait for the global peak - look for any sharp drop early on
    
    # Calculate the gradient (first derivative of density)
    # Negative values = decreasing density
    density_gradient <- diff(density_smooth)
    
    # Define significant drop as > 10% of the max absolute gradient magnitude
    max_gradient_magnitude <- max(abs(density_gradient), na.rm = TRUE)
    significant_threshold <- max_gradient_magnitude * 0.1
    
    # Find ALL significant drops (negative gradients exceeding threshold)
    significant_drops <- which(density_gradient < -significant_threshold)
    
    if (length(significant_drops) > 0) {
      # Take the FIRST significant drop
      first_drop_idx <- significant_drops[1]
      
      if (first_drop_idx > 0 && first_drop_idx < length(bin_centers)) {
        rt_threshold_candidate <- bin_centers[first_drop_idx]
      } else {
        rt_threshold_candidate <- bin_centers[1]
      }
    } else {
      # No significant drop found, use where density is lowest in first half
      first_half_idx <- seq(1, length(density_smooth) %/% 2)
      lowest_idx <- which.min(density_smooth[first_half_idx])[1]
      rt_threshold_candidate <- bin_centers[lowest_idx]
    }
    
    rt_threshold_strict <- NA
    
    # Store results
    result_key <- paste0(col_name, "-", mode_name)
    
    all_thresholds[[result_key]] <- tibble(
      Column = col_name,
      Mode = mode_name,
      BinSize_sec = bin_size,
      BinSize_min = bin_min,
      RT_threshold = rt_threshold_candidate,
      N_compounds = n_compounds,
      RT_range = paste0(round(rt_min, 2), "-", round(rt_max, 2))
    )
    
    all_density_data[[result_key]] <- tibble(
      Column = col_name,
      Mode = mode_name,
      BinSize_sec = bin_size,
      BinSize_min = bin_min,
      BinCenter = bin_centers,
      Density_raw = density_raw,
      Density_smooth = density_smooth,
      RT_threshold = rt_threshold_candidate
    )
    
    cat("    Threshold: ", round(rt_threshold_candidate, 3), "min\n")
  }
  
  # ========================================================================
  # COMBINE ALL RESULTS
  # ========================================================================
  
  thresholds_df <- bind_rows(all_thresholds)
  density_data_long <- bind_rows(all_density_data)
  
  cat("\n\n=== SUMMARY OF DETECTED THRESHOLDS ===\n")
  print(thresholds_df)
  
  # ========================================================================
  # CREATE DIAGNOSTIC VISUALIZATIONS
  # ========================================================================
  
  if (plot_diagnostics) {
    cat("\n\nGenerating diagnostic visualizations...\n")
    
    # For each column-mode combination, create comprehensive diagnostic plots
    for (i in seq_len(nrow(column_mode_combinations))) {
      col_name <- column_mode_combinations$Column[i]
      mode_name <- column_mode_combinations$Mode[i]
      
      # Filter data and density for this combination
      subset_data <- data_with_colmode %>%
        filter(Column == col_name, Mode == mode_name) %>%
        pull(RT_min) %>%
        na.omit()
      
      density_subset <- density_data_long %>%
        filter(Column == col_name, Mode == mode_name)
      
      rt_threshold <- density_subset$RT_threshold[1]
      
      # Calculate bin width for plotting
      bin_width <- unique(density_subset$BinSize_min)[1] * 0.8
      
      p <- ggplot(density_subset, aes(x = BinCenter, y = Density_smooth)) +
        geom_col(aes(y = Density_raw), fill = "lightgray", alpha = 0.6, width = bin_width) +
        geom_line(color = "steelblue", size = 1, na.rm = TRUE) +
        geom_vline(xintercept = rt_threshold, color = "red", linetype = "dashed", size = 1) +
        geom_point(aes(y = Density_smooth), color = "steelblue", size = 2) +
        labs(
          title = paste0(col_name, " - ", mode_name, " (10s bins)"),
          x = "Retention Time (min)",
          y = "Metabolite Count",
          subtitle = paste0("Threshold: ", round(rt_threshold, 3), " min")
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          axis.text = element_text(size = 10)
        )
      
      # Save figure
      filename_png <- file.path(png_path, 
                                paste0("RT_density_", col_name, "_", mode_name, ".png"))
      filename_svg <- file.path(svg_path, 
                                paste0("RT_density_", col_name, "_", mode_name, ".svg"))
      
      ggsave(filename_png, plot = p, width = 10, height = 5, dpi = 300)
      ggsave(filename_svg, plot = p, width = 10, height = 5)
      
      cat("  Saved:", basename(filename_png), "\n")
    }
    
    cat("Diagnostic plots saved to:\n")
    cat("  PNG:", png_path, "\n")
    cat("  SVG:", svg_path, "\n")
  }
  
  # Save summary threshold data as CSV
  summary_csv_file <- file.path(csv_path, "RT_threshold_summary.csv")
  write.csv(thresholds_df, summary_csv_file, row.names = FALSE)
  cat("Threshold summary saved to:", summary_csv_file, "\n")
  
  # ========================================================================
  # RETURN RESULTS
  # ========================================================================
  
  return(list(
    thresholds_summary = thresholds_df,
    density_data = density_data_long,
    recommended_thresholds = thresholds_df %>%
      select(Column, Mode, RT_cutoff = RT_threshold) %>%
      distinct()
  ))
}


# ============================================================================
# HELPER FUNCTION: Apply detected thresholds to dataset
# ============================================================================

apply_density_detected_thresholds <- function(
    dataset,
    compound_annotation,
    threshold_results
) {
  
  # Extract the recommended thresholds
  recommended_thresholds <- threshold_results$recommended_thresholds
  
  cat("\n=== APPLYING DENSITY-DETECTED THRESHOLDS ===\n")
  cat("Using", nrow(recommended_thresholds), "column-mode combinations\n\n")
  print(recommended_thresholds)
  
  # Prepare compound annotation with column-mode info
  compound_with_cutoff <- compound_annotation %>%
    mutate(
      Column = sub("-[^-]+-[0-9]+$", "", Compound.ID),
      Mode = sub(".*-(.*?)-.*", "\\1", Compound.ID),
      RT_min = `RT[min]`
    ) %>%
    left_join(recommended_thresholds, by = c("Column", "Mode")) %>%
    filter(RT_min >= RT_cutoff) %>%
    select(Compound.ID)
  
  # Filter dataset
  dataset_filtered <- dataset %>%
    filter(Compound.ID %in% compound_with_cutoff$Compound.ID)
  
  # Summary by column-mode
  filtering_summary <- compound_annotation %>%
    mutate(
      Column = sub("-[^-]+-[0-9]+$", "", Compound.ID),
      Mode = sub(".*-(.*?)-.*", "\\1", Compound.ID),
      RT_min = `RT[min]`
    ) %>%
    left_join(recommended_thresholds, by = c("Column", "Mode")) %>%
    group_by(Column, Mode, RT_cutoff) %>%
    summarise(
      total = n(),
      below_cutoff = sum(RT_min < RT_cutoff, na.rm = TRUE),
      above_cutoff = sum(RT_min >= RT_cutoff, na.rm = TRUE),
      percent_retained = round(100 * above_cutoff / total, 1),
      .groups = "drop"
    )
  
  # ========================================================================
  # PRINT COMPREHENSIVE SUMMARY
  # ========================================================================
  
  cat("\n", strrep("=", 60), "\n")
  cat("=== RT THRESHOLD FILTERING SUMMARY ===\n")
  cat(strrep("=", 60), "\n\n")
  
  cat("DETECTION METHOD: Density-based (first significant drop)\n")
  cat("BIN SIZE: 10 seconds\n\n")
  
  cat("OVERALL STATISTICS:\n")
  cat("  Total compounds before filtering:", nrow(dataset), "\n")
  cat("  Total compounds after filtering:", nrow(dataset_filtered), "\n")
  cat("  Compounds removed:", nrow(dataset) - nrow(dataset_filtered), "\n")
  cat("  Percent retained:", round(100 * nrow(dataset_filtered) / nrow(dataset), 1), "%\n\n")
  
  cat("THRESHOLDS APPLIED BY COLUMN-MODE:\n")
  cat(strrep("-", 60), "\n")
  print(filtering_summary)
  cat(strrep("-", 60), "\n\n")
  
  # Calculate detection rate if dataset has sample columns
  if (ncol(dataset) > 3) {
    sample_cols <- 4:ncol(dataset)
    total_peak_measurements <- nrow(dataset_filtered) * length(sample_cols)
    detection_rate <- (sum(!is.na(dataset_filtered[, sample_cols])) / total_peak_measurements) * 100
    
    cat("DETECTION STATISTICS (filtered dataset):\n")
    cat("  Total compounds:", nrow(dataset_filtered), "\n")
    cat("  Total samples:", length(sample_cols), "\n")
    cat("  Total peak measurements:", total_peak_measurements, "\n")
    cat("  Mean detection rate:", round(detection_rate, 2), "%\n\n")
  }
  
  cat(strrep("=", 60), "\n\n")
  
  return(list(
    filtered_dataset = dataset_filtered,
    thresholds_applied = recommended_thresholds,
    filtering_summary = filtering_summary
  ))
}


# ============================================================================
# EXAMPLE USAGE (commented out for sourcing)
# ============================================================================

# # Load your data
# source("path_to_your_data_loading_script.R")
# 
# # Detect thresholds from density with separate output paths
# threshold_results <- detect_rt_threshold_from_density(
#   compound_annotation = Batch2_Untargeted_Raw_compound_annotation,
#   min_density_drop_fraction = 0.3,
#   png_path = "../../results/untargeted/plots/png",
#   svg_path = "../../results/untargeted/plots/svg",
#   csv_path = "../../results/untargeted/data"
# )
# 
# # Or using legacy output_dir parameter (all files go to same place):
# threshold_results <- detect_rt_threshold_from_density(
#   compound_annotation = Batch2_Untargeted_Raw_compound_annotation,
#   output_dir = "../../results/untargeted/RT_threshold_diagnostics"
# )
# 
# # View recommended thresholds
# threshold_results$recommended_thresholds
# 
# # Apply to your dataset
# filtering_results <- apply_density_detected_thresholds(
#   dataset = Batch2_annotated_dataset,
#   compound_annotation = Batch2_Untargeted_Raw_compound_annotation,
#   threshold_results = threshold_results
# )
# 
# # Use filtered dataset
# Batch2_annotated_dataset_RT_filtered <- filtering_results$filtered_dataset
# METABOLITE MISSINGNESS FILTERING
# Filter out metabolites not detected in at least 75% of samples


# ============================================================================
# MAIN FILTERING FUNCTION: Remove metabolites by missingness threshold
# ============================================================================

filter_metabolites_by_missingness <- function(
    dataset,
    compound_annotation,
    detection_threshold = 0.75,  # Keep metabolites detected in ≥75% of samples
    sample_cols = NULL,
    verbose = TRUE
) {
  
  # ========================================================================
  # DETERMINE SAMPLE COLUMNS
  # ========================================================================
  
  if (is.null(sample_cols)) {
    # Assume columns 4 onward are samples
    sample_cols <- 4:ncol(dataset)
  }
  
  sample_names <- colnames(dataset)[sample_cols]
  n_samples <- length(sample_cols)
  
  if (verbose) {
    cat("\nDetermining sample columns...\n")
    cat("  Sample columns identified:", n_samples, "\n")
  }
  
  # ========================================================================
  # CALCULATE MISSINGNESS PER METABOLITE
  # ========================================================================
  
  metabolite_id_col <- colnames(dataset)[1]
  
  missingness_stats <- data.frame(
    metabolite_id = dataset[[metabolite_id_col]],
    missing_count = rowSums(is.na(dataset[, sample_cols])),
    total_samples = n_samples
  ) %>%
    mutate(
      missing_proportion = missing_count / total_samples,
      detected_samples = total_samples - missing_count,
      detection_rate = detected_samples / total_samples
    )
  
  # ========================================================================
  # IDENTIFY METABOLITES TO KEEP/REMOVE
  # ========================================================================
  
  min_detected <- detection_threshold * n_samples
  
  metabolites_to_keep <- missingness_stats %>%
    filter(detection_rate >= detection_threshold) %>%
    pull(metabolite_id)
  
  metabolites_to_remove <- missingness_stats %>%
    filter(detection_rate < detection_threshold) %>%
    pull(metabolite_id)
  
  # ========================================================================
  # FILTER DATASET
  # ========================================================================
  
  dataset_filtered <- dataset %>%
    filter(!!sym(metabolite_id_col) %in% metabolites_to_keep)
  
  # ========================================================================
  # FILTER ANNOTATION IF PROVIDED
  # ========================================================================
  
  annotation_filtered <- NULL
  if (!is.null(compound_annotation)) {
    annotation_filtered <- compound_annotation %>%
      filter(Compound.ID %in% metabolites_to_keep)
  }
  
  # ========================================================================
  # GENERATE COMPREHENSIVE SUMMARY
  # ========================================================================
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("=== MISSINGNESS-BASED FILTERING SUMMARY ===\n")
    cat(strrep("=", 70), "\n\n")
    
    cat("FILTERING CRITERIA:\n")
    cat("  Detection threshold: ≥", detection_threshold * 100, "% of samples\n")
    cat("  Minimum detected samples per metabolite:", min_detected, "\n")
    cat("  Total samples analyzed:", n_samples, "\n\n")
    
    cat("OVERALL STATISTICS:\n")
    cat("  Total metabolites before filtering:", nrow(dataset), "\n")
    cat("  Total metabolites after filtering:", nrow(dataset_filtered), "\n")
    cat("  Metabolites removed:", length(metabolites_to_remove), "\n")
    cat("  Percent retained:", round(100 * nrow(dataset_filtered) / nrow(dataset), 2), "%\n\n")
    
    # Detection rate statistics
    cat("DETECTION RATE STATISTICS (before filtering):\n")
    cat("  Mean detection rate:", round(mean(missingness_stats$detection_rate) * 100, 2), "%\n")
    cat("  Median detection rate:", round(median(missingness_stats$detection_rate) * 100, 2), "%\n")
    cat("  Min detection rate:", round(min(missingness_stats$detection_rate) * 100, 2), "%\n")
    cat("  Max detection rate:", round(max(missingness_stats$detection_rate) * 100, 2), "%\n\n")
    
    cat("DETECTION RATE STATISTICS (after filtering):\n")
    filtered_stats <- missingness_stats %>%
      filter(metabolite_id %in% metabolites_to_keep)
    cat("  Mean detection rate:", round(mean(filtered_stats$detection_rate) * 100, 2), "%\n")
    cat("  Median detection rate:", round(median(filtered_stats$detection_rate) * 100, 2), "%\n")
    cat("  Min detection rate:", round(min(filtered_stats$detection_rate) * 100, 2), "%\n")
    cat("  Max detection rate:", round(max(filtered_stats$detection_rate) * 100, 2), "%\n\n")
    
    # Missingness distribution
    cat("METABOLITES REMOVED BY DETECTION RATE:\n")
    removed_stats <- missingness_stats %>%
      filter(metabolite_id %in% metabolites_to_remove) %>%
      arrange(detection_rate)
    
    rate_bins <- c(0, 0.25, 0.50, 0.75)
    for (i in seq_len(length(rate_bins))) {
      if (i < length(rate_bins)) {
        lower <- rate_bins[i]
        upper <- rate_bins[i + 1]
        count <- sum(removed_stats$detection_rate >= lower & removed_stats$detection_rate < upper)
        cat("  ", formatC(lower * 100, width = 4, format = "f", digits = 1), "-", 
            formatC(upper * 100, width = 4, format = "f", digits = 1), "% detected: ", count, "\n")
      }
    }
    cat("\n")
    
    # Sample-wise detection rate in filtered dataset
    if (ncol(dataset_filtered) > 3) {
      sample_detection <- colSums(!is.na(dataset_filtered[, sample_cols])) / nrow(dataset_filtered)
      cat("SAMPLE-WISE DETECTION RATE (filtered dataset):\n")
      cat("  Mean across samples:", round(mean(sample_detection) * 100, 2), "%\n")
      cat("  Range:", round(min(sample_detection) * 100, 2), "-", 
          round(max(sample_detection) * 100, 2), "%\n\n")
    }
    
    cat(strrep("=", 70), "\n\n")
  }
  
  # ========================================================================
  # RETURN RESULTS
  # ========================================================================
  
  return(list(
    filtered_dataset = dataset_filtered,
    filtered_annotation = annotation_filtered,
    missingness_stats = missingness_stats,
    metabolites_removed = metabolites_to_remove,
    metabolites_kept = metabolites_to_keep,
    n_removed = length(metabolites_to_remove),
    n_kept = length(metabolites_to_keep),
    detection_threshold = detection_threshold,
    summary = list(
      n_before = nrow(dataset),
      n_after = nrow(dataset_filtered),
      n_removed = length(metabolites_to_remove),
      pct_retained = round(100 * nrow(dataset_filtered) / nrow(dataset), 2),
      mean_detection_before = round(mean(missingness_stats$detection_rate) * 100, 2),
      mean_detection_after = round(mean(filtered_stats$detection_rate) * 100, 2)
    )
  ))
}


# ============================================================================
# HELPER FUNCTION: Visualize missingness distribution
# ============================================================================

plot_missingness_distribution <- function(
    missingness_stats,
    detection_threshold = 0.75,
    output_file = NULL,
    png_path = NULL,
    svg_path = NULL
) {
  
  # Create directories if needed
  if (!is.null(png_path) && !dir.exists(png_path)) {
    dir.create(png_path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!is.null(svg_path) && !dir.exists(svg_path)) {
    dir.create(svg_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  p <- ggplot(missingness_stats, aes(x = detection_rate * 100, fill = detection_rate >= detection_threshold)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "black", size = 0.3) +
    geom_vline(xintercept = detection_threshold * 100, color = "red", linetype = "dashed", size = 1) +
    scale_fill_manual(
      values = c("TRUE" = "#1b9e77", "FALSE" = "#d73027"),
      labels = c("TRUE" = "Retained", "FALSE" = "Removed"),
      guide = "legend"
    ) +
    scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "top",
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "Metabolite Detection Rate Distribution",
      subtitle = paste0("Threshold: ≥", detection_threshold * 100, "% | ", 
                        "Retained: ", sum(missingness_stats$detection_rate >= detection_threshold), " | ",
                        "Removed: ", sum(missingness_stats$detection_rate < detection_threshold)),
      x = "Detection Rate (%)",
      y = "Number of Metabolites",
      fill = "Status"
    )
  
  # Save to separate PNG and SVG paths if provided
  if (!is.null(png_path)) {
    filename_png <- file.path(png_path, "missingness_distribution.png")
    ggsave(filename_png, plot = p, width = 10, height = 6, dpi = 300)
    cat("PNG plot saved to:", filename_png, "\n")
  }
  
  if (!is.null(svg_path)) {
    filename_svg <- file.path(svg_path, "missingness_distribution.svg")
    ggsave(filename_svg, plot = p, width = 10, height = 6)
    cat("SVG plot saved to:", filename_svg, "\n")
  }
  
  # Also support legacy output_file parameter for backward compatibility
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)
    cat("Plot saved to:", output_file, "\n")
  }
  
  return(p)
}


# ============================================================================
# EXAMPLE USAGE (commented out for sourcing)
# ============================================================================

# # Load your data and missingness calculation function
# source("path_to/missingness_functions.R")
# 
# # Filter metabolites by missingness
# missingness_results <- filter_metabolites_by_missingness(
#   dataset = Batch2_annotated_dataset_RT_filtered,
#   compound_annotation = Batch2_Untargeted_Raw_compound_annotation,
#   detection_threshold = 0.75,  # Keep metabolites in ≥75% of samples
#   verbose = TRUE
# )
# 
# # Access filtered dataset
# Batch2_annotated_dataset_missingness_filtered <- missingness_results$filtered_dataset
# 
# # Visualize the distribution with separate PNG and SVG paths
# p_missingness <- plot_missingness_distribution(
#   missingness_stats = missingness_results$missingness_stats,
#   detection_threshold = 0.75,
#   png_path = "../../results/untargeted/plots/png",
#   svg_path = "../../results/untargeted/plots/svg"
# )
# 
# # Or use legacy output_file parameter:
# # p_missingness <- plot_missingness_distribution(
# #   missingness_stats = missingness_results$missingness_stats,
# #   detection_threshold = 0.75,
# #   output_file = "../../results/untargeted/plots/missingness_distribution.png"
# # )
# 
# # View summary
# missingness_results$summary


# METABOLITE DUPLICATE DETECTION
# Identify duplicate metabolites by Name across the dataset


# ============================================================================
# MAIN FUNCTION: Detect duplicate metabolites by Name
# ============================================================================

detect_metabolite_duplicates <- function(
    dataset,
    compound_annotation = NULL,
    name_column = "Name",
    csv_path = NULL,
    verbose = TRUE
) {
  
  # ========================================================================
  # CREATE OUTPUT DIRECTORY IF NEEDED
  # ========================================================================
  
  if (!is.null(csv_path) && !dir.exists(csv_path)) {
    dir.create(csv_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ========================================================================
  # IDENTIFY DUPLICATES IN DATASET
  # ========================================================================
  
  # Get metabolite identifiers and names
  if (!is.null(compound_annotation)) {
    # Use compound annotation if provided
    metabolite_info <- compound_annotation %>%
      select(Compound.ID, all_of(name_column)) %>%
      distinct()
    
    # Join with dataset
    dataset_with_names <- dataset %>%
      select(colnames(dataset)[1]) %>%
      left_join(
        metabolite_info,
        by = setNames("Compound.ID", colnames(dataset)[1])
      )
  } else {
    # Assume dataset has a Name column
    dataset_with_names <- dataset %>%
      select(colnames(dataset)[1], all_of(name_column))
  }
  
  # Count occurrences of each metabolite name
  duplicate_summary <- dataset_with_names %>%
    group_by(!!sym(name_column)) %>%
    summarise(
      n_occurrences = n(),
      compound_ids = paste(!!sym(colnames(dataset_with_names)[1]), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(is_duplicate = n_occurrences > 1) %>%
    arrange(desc(n_occurrences))
  
  # Separate duplicates and unique
  duplicates <- duplicate_summary %>%
    filter(is_duplicate)
  
  unique_metabolites <- duplicate_summary %>%
    filter(!is_duplicate)
  
  # ========================================================================
  # SAVE DUPLICATES TO CSV IF PATH PROVIDED
  # ========================================================================
  
  csv_file <- NULL
  if (!is.null(csv_path) && nrow(duplicates) > 0) {
    csv_file <- file.path(csv_path, "metabolite_duplicates.csv")
    write.csv(duplicates, csv_file, row.names = FALSE)
  }
  
  # ========================================================================
  # GENERATE SUMMARY
  # ========================================================================
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("=== METABOLITE DUPLICATE DETECTION ===\n")
    cat(strrep("=", 70), "\n\n")
    
    cat("OVERALL STATISTICS:\n")
    cat("  Total metabolite identifiers:", nrow(dataset), "\n")
    cat("  Unique metabolite names:", nrow(duplicate_summary), "\n")
    cat("  Non-duplicate metabolites (appear once):", nrow(unique_metabolites), "\n")
    cat("  Duplicate metabolite names:", nrow(duplicates), "\n")
    cat("  Total duplicate identifiers:", nrow(dataset) - nrow(unique_metabolites), "\n\n")
    
    if (nrow(duplicates) > 0) {
      cat("DUPLICATE SUMMARY:\n")
      cat(strrep("-", 70), "\n")
      print(duplicates %>% select(!!sym(name_column), n_occurrences))
      cat(strrep("-", 70), "\n\n")
      
      if (!is.null(csv_file)) {
        cat("Detailed duplicate information saved to:", csv_file, "\n\n")
      }
    } else {
      cat("NO DUPLICATES FOUND!\n\n")
    }
    
    cat(strrep("=", 70), "\n\n")
  }
  
  # ========================================================================
  # RETURN RESULTS
  # ========================================================================
  
  return(list(
    duplicate_summary = duplicate_summary,
    duplicates = duplicates,
    unique_metabolites = unique_metabolites,
    n_total = nrow(dataset),
    n_unique_names = nrow(duplicate_summary),
    n_duplicates = nrow(duplicates),
    n_duplicate_identifiers = nrow(dataset) - nrow(unique_metabolites),
    has_duplicates = nrow(duplicates) > 0,
    csv_file_saved = csv_file
  ))
}


# ============================================================================
# HELPER FUNCTION: Quick duplicate check (minimal output)
# ============================================================================

check_duplicates <- function(dataset, compound_annotation = NULL, name_column = "Name") {
  
  result <- detect_metabolite_duplicates(
    dataset = dataset,
    compound_annotation = compound_annotation,
    name_column = name_column,
    verbose = FALSE
  )
  
  # Print minimal summary statistics only
  cat("\n=== OVERALL STATISTICS ===\n")
  cat("  Total metabolite identifiers:", result$n_total, "\n")
  cat("  Unique metabolite names:", result$n_unique_names, "\n")
  cat("  Non-duplicate metabolites (appear once):", nrow(result$unique_metabolites), "\n")
  cat("  Duplicate metabolite names:", result$n_duplicates, "\n")
  cat("  Total duplicate identifiers:", result$n_duplicate_identifiers, "\n\n")
  
  invisible(result)
}


# ============================================================================
# EXAMPLE USAGE (commented out for sourcing)
# ============================================================================

# source("detect_duplicates.R")
# 
# # Full duplicate detection with CSV output
# dup_results <- detect_metabolite_duplicates(
#   dataset = Batch2_annotated_dataset,
#   compound_annotation = Batch2_Untargeted_Raw_compound_annotation,
#   name_column = "Name",
#   csv_path = "../../results/untargeted/data",
#   verbose = TRUE
# )
# 
# # Quick check after each filtering step
# check_duplicates(Batch2_annotated_dataset_RT_filtered, name_column = "Name")
# check_duplicates(Batch2_annotated_dataset_missingness_filtered, name_column = "Name")
# 
# # Access results
# dup_results$duplicates          # Dataframe of duplicate metabolites
# dup_results$has_duplicates      # TRUE/FALSE
# dup_results$n_duplicate_identifiers  # How many IDs are duplicates
# dup_results$csv_file_saved      # Path to saved CSV file

# METABOLITE DUPLICATE REMOVAL
# Remove duplicates by selecting the metabolite with highest mean intensity

# ============================================================================
# MAIN FUNCTION: Remove duplicates by highest intensity
# ============================================================================

remove_duplicates_by_intensity <- function(
    dataset,
    compound_annotation = NULL,
    name_column = "Name",
    sample_cols = NULL,
    verbose = TRUE
) {
  
  # ========================================================================
  # DETERMINE SAMPLE COLUMNS
  # ========================================================================
  
  if (is.null(sample_cols)) {
    # Assume columns 4 onward are samples in the original dataset
    sample_cols <- 4:ncol(dataset)
    sample_names <- colnames(dataset)[sample_cols]
  } else {
    sample_names <- colnames(dataset)[sample_cols]
  }
  
  n_samples <- length(sample_cols)
  
  # ========================================================================
  # GET METABOLITE NAMES
  # ========================================================================
  
  metabolite_id_col <- colnames(dataset)[1]
  
  if (!is.null(compound_annotation)) {
    # Use compound annotation if provided
    metabolite_info <- compound_annotation %>%
      select(Compound.ID, all_of(name_column), Formula) %>%
      distinct()
    
    dataset_with_names <- dataset %>%
      select(all_of(c(metabolite_id_col, sample_names))) %>%
      left_join(
        metabolite_info,
        by = setNames("Compound.ID", metabolite_id_col)
      )
  } else {
    # Assume dataset has a Name column
    dataset_with_names <- dataset %>%
      select(all_of(c(metabolite_id_col, sample_names, name_column)))
  }
  
  # ========================================================================
  # CALCULATE MEAN INTENSITY FOR EACH METABOLITE
  # ========================================================================
  
  dataset_with_intensity <- dataset_with_names %>%
    mutate(
      mean_intensity = rowMeans(
        select(., all_of(sample_names)), 
        na.rm = TRUE
      )
    )
  
  # ========================================================================
  # IDENTIFY AND REMOVE DUPLICATES
  # ========================================================================
  
  # Group by name and select row with highest mean intensity
  dataset_deduplicated <- dataset_with_intensity %>%
    group_by(!!sym(name_column)) %>%
    arrange(desc(mean_intensity), .by_group = TRUE) %>%
    slice(1) %>%  # Keep only first row (highest intensity)
    ungroup() %>%
    select(-mean_intensity)  # Remove temporary intensity column
  
  # Get the removed metabolites for reporting
  metabolites_removed <- dataset_with_intensity %>%
    group_by(!!sym(name_column)) %>%
    filter(n() > 1) %>%  # Only groups with duplicates
    arrange(desc(mean_intensity), .by_group = TRUE) %>%
    slice(2:n()) %>%  # All except the first (highest intensity)
    ungroup() %>%
    select(!!sym(metabolite_id_col), !!sym(name_column), mean_intensity) %>%
    arrange(!!sym(name_column), desc(mean_intensity))
  
  # Count duplicate groups
  n_dup_groups <- nrow(metabolites_removed %>% distinct(!!sym(name_column)))
  
  # ========================================================================
  # REORDER COLUMNS TO MATCH ORIGINAL DATASET
  # ========================================================================
  
  # Keep the metabolite ID, annotation columns, then sample columns
  # Order should be: Compound.ID, Name, Formula, ..., then samples
  annotation_cols <- setdiff(colnames(dataset_with_intensity), c(metabolite_id_col, sample_names, "mean_intensity"))
  
  dataset_deduplicated <- dataset_deduplicated %>%
    select(all_of(c(metabolite_id_col, annotation_cols, sample_names)))
  
  # ========================================================================
  # GENERATE COMPREHENSIVE SUMMARY
  # ========================================================================
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("=== DUPLICATE REMOVAL BY INTENSITY SUMMARY ===\n")
    cat(strrep("=", 70), "\n\n")
    
    cat("OVERALL STATISTICS:\n")
    cat("  Total metabolite identifiers before:", nrow(dataset), "\n")
    cat("  Total metabolite identifiers after:", nrow(dataset_deduplicated), "\n")
    cat("  Duplicates removed:", nrow(dataset) - nrow(dataset_deduplicated), "\n")
    cat("  Percent retained:", round(100 * nrow(dataset_deduplicated) / nrow(dataset), 2), "%\n\n")
    
    # Duplicates breakdown
    cat("DUPLICATE GROUPS RESOLVED:\n")
    cat("  Metabolite names with duplicates:", n_dup_groups, "\n")
    cat("  Duplicate identifiers removed:", nrow(metabolites_removed), "\n\n")
    
    # Distribution of duplicates removed
    dup_distribution <- metabolites_removed %>%
      group_by(!!sym(name_column)) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(n) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(n)
    
    cat("DISTRIBUTION OF DUPLICATES:\n")
    for (i in seq_len(nrow(dup_distribution))) {
      n_dups <- dup_distribution$n[i]
      n_metabolites <- dup_distribution$count[i]
      cat("  ", n_metabolites, " metabolite name(s) with ", n_dups, " duplicate(s) removed\n", sep = "")
    }
    cat("\n")
    
    cat(strrep("=", 70), "\n\n")
  }
  
  # ========================================================================
  # RETURN RESULTS
  # ========================================================================
  
  return(list(
    filtered_dataset = dataset_deduplicated,
    metabolites_removed = metabolites_removed,
    n_before = nrow(dataset),
    n_after = nrow(dataset_deduplicated),
    n_removed = nrow(dataset) - nrow(dataset_deduplicated),
    pct_retained = round(100 * nrow(dataset_deduplicated) / nrow(dataset), 2),
    n_duplicate_groups = n_dup_groups,
    summary = list(
      n_before = nrow(dataset),
      n_after = nrow(dataset_deduplicated),
      n_removed = nrow(dataset) - nrow(dataset_deduplicated),
      pct_retained = round(100 * nrow(dataset_deduplicated) / nrow(dataset), 2),
      n_duplicate_groups = n_dup_groups
    )
  ))
}


# ============================================================================
# HELPER FUNCTION: Quick duplicate removal check
# ============================================================================

quick_dedup <- function(dataset, compound_annotation = NULL, name_column = "Name", sample_cols = NULL) {
  
  result <- remove_duplicates_by_intensity(
    dataset = dataset,
    compound_annotation = compound_annotation,
    name_column = name_column,
    sample_cols = sample_cols,
    verbose = FALSE
  )
  
  # Print minimal summary statistics only
  cat("\n=== OVERALL STATISTICS ===\n")
  cat("  Total metabolite identifiers before:", result$n_before, "\n")
  cat("  Total metabolite identifiers after:", result$n_after, "\n")
  cat("  Duplicates removed:", result$n_removed, "\n")
  cat("  Percent retained:", result$pct_retained, "%\n")
  cat("  Duplicate groups resolved:", result$n_duplicate_groups, "\n\n")
  
  invisible(result)
}


# ============================================================================
# EXAMPLE USAGE (commented out for sourcing)
# ============================================================================

# source("remove_duplicates_by_intensity.R")
# 
# # Remove duplicates by selecting highest intensity
# dedup_results <- remove_duplicates_by_intensity(
#   dataset = Batch2_annotated_dataset_missingness_filtered,
#   compound_annotation = Batch2_Untargeted_Raw_compound_annotation,
#   name_column = "Name",
#   verbose = TRUE
# )
# 
# # Assign the deduplicated dataset
# Batch2_annotated_dataset_final <- dedup_results$filtered_dataset
# 
# # Quick check (minimal output)
# quick_dedup(Batch2_annotated_dataset_missingness_filtered, 
#             compound_annotation = Batch2_Untargeted_Raw_compound_annotation)
# 
# # View removed metabolites
# dedup_results$metabolites_removed
# 
# # Access summary
# dedup_results$summary