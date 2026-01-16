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