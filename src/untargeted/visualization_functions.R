# Metabolomics Processing and Analysis Functions

# ============================================================================
# Calculate Missingness
# ============================================================================
# Function to calculate percentage of samples each peak is present in
calculate_missingness <- function(data) {
  # Calculate proportion of NA values for each metabolite
  missing_prop <- rowMeans(is.na(data[,-1]))
  
  # Add missingness to metabolite names
  missing_df <- data.frame(
    metabolite = data[,1],
    missing_proportion = missing_prop
  )
  
  return(missing_df)
}

# ============================================================================
# Numeric
# ============================================================================
#function to convert to numeric
convert_to_numeric <- function(x) {
  x <- as.character(x)
  x[x == "N/A"] <- NA
  as.numeric(x)}

# ============================================================================
# Impute 1/5th minimum
# ============================================================================
# Imputation function to replace NaN with 1/5th of the minimum value in each row
min_value_impute <- function(column) {
  min_value <- min(column, na.rm = TRUE)
  column[is.na(column)] <- min_value / 5
  return(column)
}
# ============================================================================
# MTAC Datafile Cleaning
# ============================================================================
# 1.  Remove asterisks from metabolite names
# 2.  Make data set completely numeric
# 3.  Remove rows that are NaN across all samples

# Data Cleaning Function for Metabolomics Data
clean_metabolomics_data <- function(peak_data, 
                                    sample_metadata, 
                                    compound_annotation, 
                                    dataset_name = "dataset",
                                    sample_col_pattern = "Area\\. ",
                                    sample_col_replacement = "Sample",
                                    first_sample_col = 4) {
  
  # Input validation
  if (!is.data.frame(peak_data)) stop("peak_data must be a data frame")
  if (!is.data.frame(sample_metadata)) stop("sample_metadata must be a data frame")
  if (!is.data.frame(compound_annotation)) stop("compound_annotation must be a data frame")
  
  cat(paste0("\n=== Cleaning ", dataset_name, " ===\n"))
  
  # Store original dimensions
  orig_rows <- nrow(peak_data)
  orig_cols <- ncol(peak_data)
  
  # Clean compound annotation
  cat("Cleaning compound annotations...\n")
  compound_annotation_clean <- compound_annotation %>%
    mutate(Name = gsub("\\*$", "", Name))
  
  # Clean peak data
  cat("Cleaning peak intensity data...\n")
  
  # Determine last sample column
  last_sample_col <- ncol(peak_data)
  
  # Clean the data
  peak_data_clean <- peak_data %>%
    # Remove asterisks from metabolite names
    mutate(Name = gsub("\\*$", "", Name)) %>%
    # Convert sample columns to numeric
    mutate(across(all_of(first_sample_col:last_sample_col), convert_to_numeric)) %>%
    # Remove rows that are completely NA across all samples
    filter(if_any(all_of(first_sample_col:last_sample_col), ~!is.na(.)))
  
  # Replace NA with NaN in sample columns for downstream functions
  peak_data_clean <- cbind(
    peak_data_clean[, 1:(first_sample_col-1)], 
    peak_data_clean[, first_sample_col:last_sample_col] %>% replace(is.na(.), NaN)
  )
  
  # Rename sample columns
  cat("Renaming sample columns...\n")
  sample_col_indices <- first_sample_col:last_sample_col
  new_col_names <- names(peak_data_clean)[sample_col_indices] %>%
    sub(sample_col_pattern, sample_col_replacement, .)
  names(peak_data_clean)[sample_col_indices] <- new_col_names
  
  # Clean sample metadata
  cat("Cleaning sample metadata...\n")
  
  # Check if metadata has 'SampleName' column and rename to 'Sample'
  if ("SampleName" %in% names(sample_metadata)) {
    sample_metadata_clean <- sample_metadata %>%
      rename(Sample = SampleName)
  } else {
    sample_metadata_clean <- sample_metadata
  }
  
  # Filter metadata to only samples present in peak data
  sample_cols_in_data <- names(peak_data_clean)[sample_col_indices]
  sample_metadata_clean <- sample_metadata_clean %>%
    filter(Sample %in% sample_cols_in_data)
  
  # Report cleaning summary
  removed_rows <- orig_rows - nrow(peak_data_clean)
  n_samples <- length(sample_col_indices)
  
  cat("\n--- Cleaning Summary ---\n")
  cat(paste0("Dataset: ", dataset_name, "\n"))
  cat(paste0("Original metabolites: ", orig_rows, "\n"))
  cat(paste0("Cleaned metabolites: ", nrow(peak_data_clean), "\n"))
  cat(paste0("Removed metabolites (all NA): ", removed_rows, "\n"))
  cat(paste0("Number of samples: ", n_samples, "\n"))
  cat(paste0("Samples in metadata: ", nrow(sample_metadata_clean), "\n"))
  
  # Check for samples in data but not in metadata
  missing_metadata <- setdiff(sample_cols_in_data, sample_metadata_clean$Sample)
  if (length(missing_metadata) > 0) {
    warning(paste0("Warning: ", length(missing_metadata), 
                   " samples in peak data are missing from metadata: ",
                   paste(missing_metadata, collapse = ", ")))
  }
  
  cat("\nCleaning complete!\n")
  
  # Create variable names based on dataset_name
  peak_var_name <- paste0(dataset_name, "_peak_data")
  metadata_var_name <- paste0(dataset_name, "_sample_metadata")
  annotation_var_name <- paste0(dataset_name, "_compound_annotation")
  
  # Assign to global environment with custom names
  assign(peak_var_name, peak_data_clean, envir = .GlobalEnv)
  assign(metadata_var_name, sample_metadata_clean, envir = .GlobalEnv)
  assign(annotation_var_name, compound_annotation_clean, envir = .GlobalEnv)
  
  cat("\nCleaned data assigned to:\n")
  cat(paste0("  - ", peak_var_name, "\n"))
  cat(paste0("  - ", metadata_var_name, "\n"))
  cat(paste0("  - ", annotation_var_name, "\n"))
  
  # Return list of cleaned data (also available via assignment)
  invisible(list(
    peak_data = peak_data_clean,
    sample_metadata = sample_metadata_clean,
    compound_annotation = compound_annotation_clean,
    dataset_name = dataset_name,
    cleaning_info = list(
      original_metabolites = orig_rows,
      cleaned_metabolites = nrow(peak_data_clean),
      removed_metabolites = removed_rows,
      n_samples = n_samples,
      first_sample_col = first_sample_col,
      last_sample_col = last_sample_col
    )
  ))
}
# ============================================================================
# Data Preparation Function
# ============================================================================

#' Convert wide-format metabolomics data to long format with metadata
#'
#' @param dataset Wide-format dataset with metabolite rows and sample columns
#' @param sample_metadata Metadata dataframe with Sample, Sex, GW columns
#' @param sample_cols Column indices or names containing sample data (default: 4:27)
#' @param sample_pattern Optional regex pattern to identify sample columns (e.g., "^Sample")
#' @param log_transform Logical, whether to log10 transform intensity values (default: TRUE)
#' @param remove_nonfinite Logical, whether to remove non-finite values (default: TRUE)
#' @param gw_buckets Named list of GW bucket definitions (default: 16-18, 19-21, 22-24)
#' @return Dataframe in long format with intensity values and metadata
prepare_long_data <- function(dataset, 
                              sample_metadata,
                              data_name = "data_long",
                              sample_cols = NULL,
                              sample_pattern = "^Sample",
                              log_transform = TRUE,
                              remove_nonfinite = TRUE,
                              gw_buckets = list(
                                "16-18" = c(16, 18),
                                "19-21" = c(19, 21),
                                "22-24" = c(22, 24)
                              )) {
  
  # Determine sample columns
  if (is.null(sample_cols)) {
    # Use pattern to find sample columns
    sample_cols <- grep(sample_pattern, names(dataset), value = FALSE)
    if (length(sample_cols) == 0) {
      stop("No sample columns found. Please specify sample_cols or adjust sample_pattern.")
    }
  }
  
  # Get sample names
  sample_names <- names(dataset)[sample_cols]
  
  cat("Preparing long-format data...\n")
  cat("  Found", length(sample_names), "samples\n")
  cat("  Dataset contains", nrow(dataset), "metabolites\n")
  
  # Reshape to long format
  data_long <- dataset %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "Sample",
      values_to = "Intensity"
    )
  
  cat("  Initial observations:", nrow(data_long), "\n")
  
  # Remove non-finite values if requested
  if (remove_nonfinite) {
    n_before <- nrow(data_long)
    data_long <- data_long %>%
      filter(is.finite(Intensity))
    n_removed <- n_before - nrow(data_long)
    if (n_removed > 0) {
      cat("  Removed", n_removed, "non-finite values\n")
    }
  }
  
  # Log-transform if requested
  if (log_transform) {
    data_long <- data_long %>%
      mutate(log10 = log10(Intensity))
    cat("  Log10 transformed intensity values\n")
  }
  
  # Filter metadata to samples present in dataset
  sample_metadata_filtered <- sample_metadata %>%
    filter(Sample %in% sample_names)
  
  if (nrow(sample_metadata_filtered) != length(sample_names)) {
    warning(sprintf(
      "Metadata available for %d/%d samples. Some samples may lack metadata.",
      nrow(sample_metadata_filtered),
      length(sample_names)
    ))
  }
  
  # Create GW buckets if not already present
  if (!"GW_bucket" %in% names(sample_metadata_filtered) && "GW" %in% names(sample_metadata_filtered)) {
    
    # Build case_when conditions from gw_buckets list
    bucket_conditions <- purrr::map2(
      names(gw_buckets),
      gw_buckets,
      ~quo(GW >= !!.y[1] & GW <= !!.y[2] ~ !!.x)
    )
    
    sample_metadata_filtered <- sample_metadata_filtered %>%
      mutate(GW_bucket = case_when(
        GW >= gw_buckets[[1]][1] & GW <= gw_buckets[[1]][2] ~ names(gw_buckets)[1],
        GW >= gw_buckets[[2]][1] & GW <= gw_buckets[[2]][2] ~ names(gw_buckets)[2],
        GW >= gw_buckets[[3]][1] & GW <= gw_buckets[[3]][2] ~ names(gw_buckets)[3],
        TRUE ~ NA_character_
      ))
    
    cat("  Created GW buckets:", paste(names(gw_buckets), collapse = ", "), "\n")
  }
  
  # Select metadata columns to join (all available relevant columns)
  metadata_cols <- intersect(
    names(sample_metadata_filtered),
    c("Sample", "Sex", "GW", "GW_bucket", "Batch", "Group", "Condition")
  )
  
  # Add metadata to long format data
  data_long <- data_long %>%
    left_join(
      sample_metadata_filtered %>% select(all_of(metadata_cols)),
      by = "Sample"
    )
  
  cat("  Added metadata columns:", paste(setdiff(metadata_cols, "Sample"), collapse = ", "), "\n")
  cat("  Final observations:", nrow(data_long), "\n")
  
  # Check for missing metadata
  if (any(is.na(data_long$Sex)) || any(is.na(data_long$GW))) {
    warning("Some observations are missing Sex or GW metadata")
  }
  
  # Assign to global environment with custom names
  assign(data_name, data_long, envir = .GlobalEnv)
  
  cat("Done!\n\n")
  
  return(data_long)
}

#' Create custom GW buckets for metadata
#'
#' @param metadata Metadata dataframe with GW column
#' @param buckets Named list of GW bucket definitions
#' @return Metadata with GW_bucket column added
create_gw_buckets <- function(metadata, 
                              buckets = list(
                                "16-18" = c(16, 18),
                                "19-21" = c(19, 21),
                                "22-24" = c(22, 24)
                              )) {
  
  if (!"GW" %in% names(metadata)) {
    stop("Metadata must contain a 'GW' column")
  }
  
  metadata <- metadata %>%
    mutate(GW_bucket = case_when(
      GW >= buckets[[1]][1] & GW <= buckets[[1]][2] ~ names(buckets)[1],
      GW >= buckets[[2]][1] & GW <= buckets[[2]][2] ~ names(buckets)[2],
      GW >= buckets[[3]][1] & GW <= buckets[[3]][2] ~ names(buckets)[3],
      TRUE ~ NA_character_
    ))
  
  return(metadata)
}# ============================================================================
# Helper Function for Saving Plots
# ============================================================================

#' Save a ggplot as both SVG and PNG to separate directories
#'
#' @param plot ggplot object to save
#' @param filename Base filename (without extension)
#' @param png_path Directory path to save PNG files
#' @param svg_path Directory path to save SVG files
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @param dpi Resolution for PNG (default: 300)

# ============================================================================
# Image Saving
# ============================================================================
save_dual_format <- function(plot, filename, png_path = NULL, svg_path = NULL,
                             width = 10, height = 8, dpi = 300) {
  # Save PNG
  if (!is.null(png_path)) {
    ggsave(
      filename = paste0(filename, ".png"),
      plot = plot,
      path = png_path,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
  
  # Save SVG
  if (!is.null(svg_path)) {
    ggsave(
      filename = paste0(filename, ".svg"),
      plot = plot,
      path = svg_path,
      width = width,
      height = height,
      bg = "transparent"
    )
  }
}

# ============================================================================
# Peak Intensity Visualization Functions
# ============================================================================

plot_intensity_histogram <- function(data_long, title_prefix = "", 
                                     png_path = NULL, svg_path = NULL,
                                     width = 10, height = 8, dpi = 300) {
  p <- ggplot(data_long, aes(x = log10)) +
    geom_histogram(binwidth = 0.4, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Peak Intensity"),
         x = "Log10 Intensity",
         y = "Number of Peaks")
  
  save_dual_format(p, paste0(title_prefix, "intensity_histogram"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_intensity_histogram_by_sex <- function(data_long, title_prefix = "", 
                                            png_path = NULL, svg_path = NULL,
                                            width = 10, height = 10, dpi = 300) {
  p <- ggplot(data_long, aes(x = log10, fill = Sex)) +
    geom_histogram(binwidth = 0.4, color = "black", position = "identity", alpha = 0.6) +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Peak Intensity by Sex"),
         x = "Log10 Intensity",
         y = "Number of Peaks") +
    facet_wrap(~Sex, ncol = 1)
  
  save_dual_format(p, paste0(title_prefix, "intensity_histogram_by_sex"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_intensity_violin <- function(data_long, title_prefix = "", 
                                  png_path = NULL, svg_path = NULL,
                                  width = 12, height = 8, dpi = 300) {
  p <- ggplot(data_long, aes(x = as.factor(Sample), y = log10)) +
    geom_violin(fill = "darkorange", color = "black") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Peak Intensity per Sample"),
         x = "Sample",
         y = "Log10 Intensity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  save_dual_format(p, paste0(title_prefix, "intensity_violin"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_intensity_violin_by_sex <- function(data_long, title_prefix = "", 
                                         png_path = NULL, svg_path = NULL,
                                         width = 12, height = 8, dpi = 300) {
  p <- ggplot(data_long, aes(x = as.factor(Sample), y = log10, fill = Sex)) +
    geom_violin(color = "black") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Peak Intensity per Sample by Sex"),
         x = "Sample",
         y = "Log10 Intensity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  save_dual_format(p, paste0(title_prefix, "intensity_violin_by_sex"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_intensity_vs_gw <- function(data_long, title_prefix = "", 
                                 png_path = NULL, svg_path = NULL,
                                 width = 10, height = 8, dpi = 300) {
  # Calculate mean intensity per sample
  mean_intensity_per_sample <- data_long %>%
    group_by(Sample, Sex, GW) %>%
    summarise(mean_log10_intensity = mean(log10, na.rm = TRUE),
              .groups = "drop")
  
  p <- ggplot(mean_intensity_per_sample, aes(x = GW, y = mean_log10_intensity, color = Sex)) +
    geom_point(size = 3) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Mean Peak Intensity vs Gestational Week"),
         x = "Gestational Week",
         y = "Mean Log10 Intensity")
  
  save_dual_format(p, paste0(title_prefix, "intensity_vs_GW"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

# Wrapper function to generate all intensity plots
plot_all_intensity <- function(data_long, title_prefix = "", 
                               png_path = NULL, svg_path = NULL,
                               width = 10, height = 8, dpi = 300) {
  plots <- list(
    histogram = plot_intensity_histogram(data_long, title_prefix, png_path, svg_path, width, height, dpi),
    histogram_by_sex = plot_intensity_histogram_by_sex(data_long, title_prefix, png_path, svg_path, width, height + 2, dpi),
    violin = plot_intensity_violin(data_long, title_prefix, png_path, svg_path, width + 2, height, dpi),
    violin_by_sex = plot_intensity_violin_by_sex(data_long, title_prefix, png_path, svg_path, width + 2, height, dpi),
    vs_gw = plot_intensity_vs_gw(data_long, title_prefix, png_path, svg_path, width, height, dpi)
  )
  return(plots)
}
# ============================================================================
# Peak Detection Calculation Functions
# ============================================================================

#' Calculate peak detection metrics from dataset
#'
#' @param dataset Wide-format dataset with metabolite rows and sample columns
#' @param sample_metadata Metadata dataframe with Sample, Sex, GW columns
#' @param sample_cols Column indices or names containing sample data
#' @param gw_buckets Named list of GW bucket definitions
#' @return Dataframe with peak detection metrics per sample
calculate_peak_detection <- function(dataset, 
                                     sample_metadata,
                                     sample_cols = NULL,
                                     sample_pattern = "^Sample",
                                     gw_buckets = list(
                                       "16-18" = c(16, 18),
                                       "19-21" = c(19, 21),
                                       "22-24" = c(22, 24)
                                     )) {
  
  # Determine sample columns
  if (is.null(sample_cols)) {
    sample_cols <- grep(sample_pattern, names(dataset), value = FALSE)
    if (length(sample_cols) == 0) {
      stop("No sample columns found. Please specify sample_cols or adjust sample_pattern.")
    }
  }
  
  # Get sample names
  sample_names <- names(dataset)[sample_cols]
  
  cat("Calculating peak detection metrics...\n")
  cat("  Samples:", length(sample_names), "\n")
  cat("  Total metabolites:", nrow(dataset), "\n")
  
  # Calculate detection metrics for each sample
  peak_detection_df <- data.frame(
    Sample = sample_names,
    detected_peaks = colSums(!is.na(dataset[, sample_cols])),
    total_metabolites = nrow(dataset)
  ) %>%
    mutate(detection_rate = detected_peaks / total_metabolites)
  
  # Filter metadata to samples present in dataset
  sample_metadata_filtered <- sample_metadata %>%
    filter(Sample %in% sample_names)
  
  # Create GW buckets if not already present
  if (!"GW_bucket" %in% names(sample_metadata_filtered) && "GW" %in% names(sample_metadata_filtered)) {
    sample_metadata_filtered <- sample_metadata_filtered %>%
      mutate(GW_bucket = case_when(
        GW >= gw_buckets[[1]][1] & GW <= gw_buckets[[1]][2] ~ names(gw_buckets)[1],
        GW >= gw_buckets[[2]][1] & GW <= gw_buckets[[2]][2] ~ names(gw_buckets)[2],
        GW >= gw_buckets[[3]][1] & GW <= gw_buckets[[3]][2] ~ names(gw_buckets)[3],
        TRUE ~ NA_character_
      ))
  }
  
  # Select metadata columns to join
  metadata_cols <- intersect(
    names(sample_metadata_filtered),
    c("Sample", "Sex", "GW", "GW_bucket", "Batch", "Group", "Condition")
  )
  
  # Add metadata
  peak_detection_df <- peak_detection_df %>%
    left_join(
      sample_metadata_filtered %>% select(all_of(metadata_cols)),
      by = "Sample"
    )
  
  cat("  Mean detected peaks:", round(mean(peak_detection_df$detected_peaks), 1), "\n")
  cat("  Mean detection rate:", scales::percent(mean(peak_detection_df$detection_rate)), "\n")
  
  if ("Sex" %in% names(peak_detection_df)) {
    cat("\nBy Sex:\n")
    sex_summary <- peak_detection_df %>%
      group_by(Sex) %>%
      summarise(
        n_samples = n(),
        mean_detected = round(mean(detected_peaks), 1),
        mean_rate = scales::percent(mean(detection_rate))
      )
    print(sex_summary)
  }
  
  cat("Done!\n\n")
  
  return(peak_detection_df)
}

#' Calculate metabolite missingness statistics
#'
#' @param dataset Wide-format dataset with metabolite rows and sample columns
#' @param sample_cols Column indices or names containing sample data (optional)
#' @return Dataframe with missingness statistics per metabolite
calculate_missingness <- function(dataset, sample_cols = NULL) {
  
  # If sample_cols not specified, find all numeric columns except first few metadata columns
  if (is.null(sample_cols)) {
    # Assume columns 4 onward are samples, or use all numeric columns
    numeric_cols <- sapply(dataset, is.numeric)
    sample_cols <- which(numeric_cols)
    if (ncol(dataset) > 3 && all(numeric_cols[4:ncol(dataset)])) {
      sample_cols <- 4:ncol(dataset)
    }
  }
  
  # Get metabolite identifiers (assuming first column)
  metabolite_id_col <- names(dataset)[1]
  
  # Calculate missingness for each metabolite
  missing_stats <- data.frame(
    metabolite_id = dataset[[metabolite_id_col]],
    missing_count = rowSums(is.na(dataset[, sample_cols])),
    total_samples = length(sample_cols)
  ) %>%
    mutate(
      missing_proportion = missing_count / total_samples,
      detected_samples = total_samples - missing_count,
      detection_rate = detected_samples / total_samples
    )
  
  return(missing_stats)
}

#' Calculate missingness by sex
#'
#' @param dataset Wide-format dataset with metabolite rows and sample columns
#' @param sample_metadata Metadata dataframe with Sample and Sex columns
#' @param sample_cols Column indices or names containing sample data
#' @return Dataframe with missingness statistics by sex
calculate_missingness_by_sex <- function(dataset, sample_metadata, sample_cols = NULL) {
  
  # Determine sample columns
  if (is.null(sample_cols)) {
    sample_cols <- 4:ncol(dataset)
  }
  
  sample_names <- names(dataset)[sample_cols]
  
  # Filter metadata to available samples
  sample_metadata_filtered <- sample_metadata %>%
    filter(Sample %in% sample_names)
  
  # Split samples by sex
  female_samples <- sample_metadata_filtered %>% filter(Sex == "F") %>% pull(Sample)
  male_samples <- sample_metadata_filtered %>% filter(Sex == "M") %>% pull(Sample)
  
  # Calculate missingness for female samples
  female_cols <- which(names(dataset) %in% female_samples)
  missing_female <- calculate_missingness(dataset, female_cols)
  missing_female$Sex <- "F"
  
  # Calculate missingness for male samples
  male_cols <- which(names(dataset) %in% male_samples)
  missing_male <- calculate_missingness(dataset, male_cols)
  missing_male$Sex <- "M"
  
  # Combine
  missing_by_sex <- bind_rows(missing_female, missing_male)
  
  return(missing_by_sex)
}

#' Print summary of peak detection metrics
#'
#' @param peak_detection_df Dataframe from calculate_peak_detection
#' @return NULL (prints to console)
summarize_peak_detection <- function(peak_detection_df) {
  
  cat("\n=== PEAK DETECTION SUMMARY ===\n")
  cat("Total samples:", nrow(peak_detection_df), "\n")
  cat("Mean detected peaks:", round(mean(peak_detection_df$detected_peaks), 1), "\n")
  cat("SD detected peaks:", round(sd(peak_detection_df$detected_peaks), 1), "\n")
  cat("Range:", min(peak_detection_df$detected_peaks), "-", 
      max(peak_detection_df$detected_peaks), "\n")
  cat("Mean detection rate:", scales::percent(mean(peak_detection_df$detection_rate)), "\n")
  
  if ("Sex" %in% names(peak_detection_df)) {
    cat("\n=== BY SEX ===\n")
    sex_summary <- peak_detection_df %>%
      group_by(Sex) %>%
      summarise(
        N = n(),
        Mean_Detected = round(mean(detected_peaks), 1),
        SD_Detected = round(sd(detected_peaks), 1),
        Mean_Rate = scales::percent(mean(detection_rate))
      )
    print(sex_summary)
  }
  
  if ("GW_bucket" %in% names(peak_detection_df)) {
    cat("\n=== BY GW BUCKET ===\n")
    gw_summary <- peak_detection_df %>%
      group_by(GW_bucket) %>%
      summarise(
        N = n(),
        Mean_Detected = round(mean(detected_peaks), 1),
        SD_Detected = round(sd(detected_peaks), 1),
        Mean_Rate = scales::percent(mean(detection_rate))
      )
    print(gw_summary)
  }
  
  cat("\n")
}

# ============================================================================
# Updated Peak Detection Visualization Functions
# ============================================================================

plot_detected_peaks_histogram <- function(peak_detection_df, title_prefix = "", 
                                          png_path = NULL, svg_path = NULL,
                                          width = 10, height = 8, dpi = 300) {
  
  if (!"detected_peaks" %in% names(peak_detection_df)) {
    stop("peak_detection_df must contain 'detected_peaks' column. Use calculate_peak_detection() first.")
  }
  
  p <- ggplot(peak_detection_df, aes(x = detected_peaks)) +
    geom_histogram(binwidth = 10, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Detected Peaks Across Samples"),
         x = "Number of Detected Peaks",
         y = "Number of Samples") +
    geom_vline(aes(xintercept = mean(detected_peaks)), 
               color = "red", linetype = "dashed", size = 1)
  
  save_dual_format(p, paste0(title_prefix, "detected_peaks_histogram"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_detected_peaks_by_sex <- function(peak_detection_df, title_prefix = "", 
                                       png_path = NULL, svg_path = NULL,
                                       width = 10, height = 10, dpi = 300) {
  
  if (!"detected_peaks" %in% names(peak_detection_df) || !"Sex" %in% names(peak_detection_df)) {
    stop("peak_detection_df must contain 'detected_peaks' and 'Sex' columns")
  }
  
  p <- ggplot(peak_detection_df, aes(x = detected_peaks, fill = Sex)) +
    geom_histogram(binwidth = 10, color = "black", position = "identity", alpha = 0.6) +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Distribution of Detected Peaks by Sex"),
         x = "Number of Detected Peaks",
         y = "Number of Samples") +
    facet_wrap(~Sex, ncol = 1)
  
  save_dual_format(p, paste0(title_prefix, "detected_peaks_by_sex"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_detected_peaks_per_sample <- function(peak_detection_df, title_prefix = "", 
                                           png_path = NULL, svg_path = NULL,
                                           width = 12, height = 8, dpi = 300) {
  
  if (!"detected_peaks" %in% names(peak_detection_df)) {
    stop("peak_detection_df must contain 'detected_peaks' column")
  }
  
  p <- ggplot(peak_detection_df, aes(x = Sample, y = detected_peaks, fill = Sex)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0(title_prefix, "Number of Detected Peaks per Sample"),
         x = "Sample",
         y = "Number of Detected Peaks")
  
  save_dual_format(p, paste0(title_prefix, "detected_peaks_per_sample"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_detected_peaks_vs_gw <- function(peak_detection_df, title_prefix = "", 
                                      png_path = NULL, svg_path = NULL,
                                      width = 10, height = 8, dpi = 300) {
  
  if (!"detected_peaks" %in% names(peak_detection_df) || !"GW" %in% names(peak_detection_df)) {
    stop("peak_detection_df must contain 'detected_peaks' and 'GW' columns")
  }
  
  p <- ggplot(peak_detection_df, aes(x = GW, y = detected_peaks, color = Sex)) +
    geom_point(size = 3) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Detected Peaks vs Gestational Week"),
         x = "Gestational Week",
         y = "Number of Detected Peaks")
  
  save_dual_format(p, paste0(title_prefix, "detected_peaks_vs_GW"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

plot_detection_rate_vs_gw <- function(peak_detection_df, title_prefix = "", 
                                      png_path = NULL, svg_path = NULL,
                                      width = 10, height = 8, dpi = 300) {
  
  if (!"detection_rate" %in% names(peak_detection_df) || !"GW" %in% names(peak_detection_df)) {
    stop("peak_detection_df must contain 'detection_rate' and 'GW' columns")
  }
  
  p <- ggplot(peak_detection_df, aes(x = GW, y = detection_rate, color = Sex)) +
    geom_point(size = 3) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Detection Rate vs Gestational Week"),
         x = "Gestational Week",
         y = "Detection Rate") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  
  save_dual_format(p, paste0(title_prefix, "detection_rate_vs_GW"), 
                   png_path, svg_path, width, height, dpi)
  return(p)
}

# Wrapper function to generate all detection plots
plot_all_detection <- function(peak_detection_df, title_prefix = "", 
                               png_path = NULL, svg_path = NULL,
                               width = 10, height = 8, dpi = 300) {
  plots <- list(
    histogram = plot_detected_peaks_histogram(peak_detection_df, title_prefix, png_path, svg_path, width, height, dpi),
    by_sex = plot_detected_peaks_by_sex(peak_detection_df, title_prefix, png_path, svg_path, width, height + 2, dpi),
    per_sample = plot_detected_peaks_per_sample(peak_detection_df, title_prefix, png_path, svg_path, width + 2, height, dpi),
    vs_gw = plot_detected_peaks_vs_gw(peak_detection_df, title_prefix, png_path, svg_path, width, height, dpi),
    rate_vs_gw = plot_detection_rate_vs_gw(peak_detection_df, title_prefix, png_path, svg_path, width, height, dpi)
  )
  return(plots)
}

# ============================================================================
# Peak Features Statistical Testing Functions
# ============================================================================

#' Perform comprehensive statistical tests on sample-level peak metrics
#'
#' @param dataset Full dataset with metabolite intensity values (columns 4+ are samples)
#' @param sample_metadata Metadata with Sample, Sex, GW columns
#' @param sample_cols Column indices or names containing sample data
#' @param save_path Optional path to save results as CSV
#' @param title_prefix Optional prefix for output filenames
#' @return List containing sample_metrics dataframe and summary_results dataframe
test_peak_features <- function(dataset, 
                               sample_metadata, 
                               sample_cols = 4:27,
                               save_path = NULL,
                               title_prefix = "") {
  
  # Filter metadata to samples present in dataset
  sample_names <- names(dataset)[sample_cols]
  sample_metadata_filtered <- sample_metadata %>%
    filter(Sample %in% sample_names)
  
  # Create GW buckets if not already present
  if (!"GW_bucket" %in% names(sample_metadata_filtered)) {
    sample_metadata_filtered <- sample_metadata_filtered %>%
      mutate(GW_bucket = case_when(
        GW >= 16 & GW <= 18 ~ "16-18",
        GW >= 19 & GW <= 21 ~ "19-21",
        GW >= 22 & GW <= 24 ~ "22-24",
        TRUE ~ NA_character_
      ))
  }
  
  # Calculate sample-level metrics
  sample_metrics <- data.frame(
    Sample = sample_names,
    detected_peaks = colSums(!is.na(dataset[, sample_cols])),
    detection_rate = colSums(!is.na(dataset[, sample_cols])) / nrow(dataset),
    mean_intensity = colMeans(dataset[, sample_cols], na.rm = TRUE),
    missing_count = colSums(is.na(dataset[, sample_cols])),
    missing_proportion = colSums(is.na(dataset[, sample_cols])) / nrow(dataset)
  ) %>%
    left_join(sample_metadata_filtered %>% select(Sample, Sex, GW, GW_bucket), 
              by = "Sample")
  
  # Define metrics to test
  metrics <- c("detected_peaks", "detection_rate", "mean_intensity", 
               "missing_proportion", "missing_count")
  metric_labels <- c("Detected Peaks", "Detection Rate", "Mean Intensity", 
                     "Missing Proportion", "Missing Count")
  
  # Initialize results list
  all_models <- list()
  summary_results <- data.frame()
  
  # Run tests for each metric
  for (i in seq_along(metrics)) {
    metric <- metrics[i]
    label <- metric_labels[i]
    
    cat("\n=== ", toupper(label), " ANALYSIS ===\n", sep = "")
    
    # Continuous GW model
    formula_cont <- as.formula(paste(metric, "~ Sex * GW"))
    model_cont <- lm(formula_cont, data = sample_metrics)
    
    cat("\nLinear Model (continuous GW):\n")
    print(summary(model_cont))
    cat("\nANOVA:\n")
    print(anova(model_cont))
    
    # GW bins model
    formula_bins <- as.formula(paste(metric, "~ Sex * GW_bucket"))
    model_bins <- lm(formula_bins, data = sample_metrics)
    
    cat("\nLinear Model (GW bins):\n")
    print(summary(model_bins))
    cat("\nANOVA:\n")
    print(anova(model_bins))
    
    # Store models
    all_models[[paste0(metric, "_continuous")]] <- model_cont
    all_models[[paste0(metric, "_bins")]] <- model_bins
    
    # Extract p-values for summary table
    # Continuous model
    coef_cont <- summary(model_cont)$coefficients
    p_sex_cont <- if ("SexM" %in% rownames(coef_cont)) coef_cont["SexM", "Pr(>|t|)"] else NA
    p_gw_cont <- if ("GW" %in% rownames(coef_cont)) coef_cont["GW", "Pr(>|t|)"] else NA
    p_int_cont <- if ("SexM:GW" %in% rownames(coef_cont)) coef_cont["SexM:GW", "Pr(>|t|)"] else NA
    
    # Bins model
    anova_bins <- anova(model_bins)
    p_sex_bins <- if ("Sex" %in% rownames(anova_bins)) anova_bins["Sex", "Pr(>F)"] else NA
    p_gw_bins <- if ("GW_bucket" %in% rownames(anova_bins)) anova_bins["GW_bucket", "Pr(>F)"] else NA
    p_int_bins <- if ("Sex:GW_bucket" %in% rownames(anova_bins)) anova_bins["Sex:GW_bucket", "Pr(>F)"] else NA
    
    # Add to summary results
    summary_results <- bind_rows(
      summary_results,
      data.frame(
        Variable = rep(label, 6),
        Model = rep(c("Continuous GW", "GW Bins"), each = 3),
        Term = rep(c("Sex", "GW/GW_bucket", "Sex:GW/Sex:GW_bucket"), 2),
        P_value = c(p_sex_cont, p_gw_cont, p_int_cont,
                    p_sex_bins, p_gw_bins, p_int_bins)
      )
    )
  }
  
  # Add significance annotations
  summary_results <- summary_results %>%
    mutate(
      Significant = case_when(
        is.na(P_value) ~ "NA",
        P_value < 0.001 ~ "***",
        P_value < 0.01 ~ "**",
        P_value < 0.05 ~ "*",
        P_value < 0.1 ~ ".",
        TRUE ~ "NS"
      ),
      P_value_formatted = ifelse(is.na(P_value), "NA", 
                                 sprintf("%.4f", P_value))
    )
  
  cat("\n=== SUMMARY TABLE OF ALL STATISTICAL TESTS ===\n")
  print(summary_results)
  
  # Save results if path provided
  if (!is.null(save_path)) {
    write.csv(summary_results, 
              file.path(save_path, paste0(title_prefix, "statistical_tests_summary.csv")), 
              row.names = FALSE)
    
    write.csv(sample_metrics, 
              file.path(save_path, paste0(title_prefix, "sample_metrics_for_testing.csv")), 
              row.names = FALSE)
    
    cat("\nResults saved to:", save_path, "\n")
  }
  
  # Return results
  return(list(
    sample_metrics = sample_metrics,
    summary_results = summary_results,
    models = all_models
  ))
}

#' Run statistical tests on a single peak metric
#'
#' @param sample_metrics Dataframe with sample-level metrics and metadata
#' @param metric_name Column name of the metric to test
#' @param metric_label Human-readable label for output
#' @return List containing continuous and bins models plus p-values
test_single_metric <- function(sample_metrics, 
                               metric_name, 
                               metric_label = metric_name) {
  
  cat("\n=== ", toupper(metric_label), " ANALYSIS ===\n", sep = "")
  
  # Continuous GW model
  formula_cont <- as.formula(paste(metric_name, "~ Sex * GW"))
  model_cont <- lm(formula_cont, data = sample_metrics)
  
  cat("\nLinear Model (continuous GW):\n")
  print(summary(model_cont))
  cat("\nANOVA:\n")
  print(anova(model_cont))
  
  # GW bins model
  formula_bins <- as.formula(paste(metric_name, "~ Sex * GW_bucket"))
  model_bins <- lm(formula_bins, data = sample_metrics)
  
  cat("\nLinear Model (GW bins):\n")
  print(summary(model_bins))
  cat("\nANOVA:\n")
  print(anova(model_bins))
  
  # Extract p-values
  coef_cont <- summary(model_cont)$coefficients
  anova_bins <- anova(model_bins)
  
  results <- list(
    model_continuous = model_cont,
    model_bins = model_bins,
    p_values_continuous = list(
      Sex = if ("SexM" %in% rownames(coef_cont)) coef_cont["SexM", "Pr(>|t|)"] else NA,
      GW = if ("GW" %in% rownames(coef_cont)) coef_cont["GW", "Pr(>|t|)"] else NA,
      Interaction = if ("SexM:GW" %in% rownames(coef_cont)) coef_cont["SexM:GW", "Pr(>|t|)"] else NA
    ),
    p_values_bins = list(
      Sex = if ("Sex" %in% rownames(anova_bins)) anova_bins["Sex", "Pr(>F)"] else NA,
      GW_bucket = if ("GW_bucket" %in% rownames(anova_bins)) anova_bins["GW_bucket", "Pr(>F)"] else NA,
      Interaction = if ("Sex:GW_bucket" %in% rownames(anova_bins)) anova_bins["Sex:GW_bucket", "Pr(>F)"] else NA
    )
  )
  
  return(results)
}

#' Format p-value with significance stars
#'
#' @param p_value Numeric p-value
#' @return Formatted string with p-value and significance annotation
format_p_value <- function(p_value) {
  if (is.na(p_value)) return("NA")
  
  sig <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ ".",
    TRUE ~ "NS"
  )
  
  sprintf("%.4f %s", p_value, sig)
}

# ============================================================================
# Metabolite Abundance Calculation Functions
# ============================================================================

#' Calculate mean intensity for each metabolite across all samples
#'
#' @param dataset Wide-format dataset with metabolite rows and sample columns
#' @param sample_cols Column indices or names containing sample data
#' @param compound_annotation Optional annotation dataframe to merge
#' @return Dataframe with metabolite abundance and annotations
calculate_metabolite_abundance <- function(dataset, 
                                           sample_cols = NULL,
                                           sample_pattern = "^Sample",
                                           compound_annotation = NULL) {
  
  # Determine sample columns
  if (is.null(sample_cols)) {
    sample_cols <- grep(sample_pattern, names(dataset), value = TRUE)
    if (length(sample_cols) == 0) {
      stop("No sample columns found. Please specify sample_cols or adjust sample_pattern.")
    }
  } else if (is.numeric(sample_cols)) {
    sample_cols <- names(dataset)[sample_cols]
  }
  
  cat("Calculating metabolite abundance...\n")
  cat("  Samples:", length(sample_cols), "\n")
  cat("  Metabolites:", nrow(dataset), "\n")
  
  # Calculate mean intensity across samples
  metabolite_abundance <- dataset %>%
    mutate(mean_intensity = rowMeans(select(., all_of(sample_cols)), na.rm = TRUE)) %>%
    select(Compound.ID, Name, Formula, mean_intensity) %>%
    arrange(desc(mean_intensity)) %>%
    mutate(display_name = ifelse(is.na(Name) | Name == "", 
                                 Compound.ID, 
                                 paste0(Name, " (", Compound.ID, ")")))
  
  # Add annotations if provided
  if (!is.null(compound_annotation)) {
    metabolite_abundance <- metabolite_abundance %>%
      left_join(
        compound_annotation %>% select(Compound.ID, Name, Class, Sub.Class, Super.Class, Pathways),
        by = c("Compound.ID", "Name")
      )
    cat("  Added compound annotations\n")
  }
  
  cat("Done!\n\n")
  
  return(metabolite_abundance)
}

#' Calculate mean intensity for each metabolite stratified by a grouping variable
#'
#' @param data_long Long-format data with Intensity column
#' @param group_by Character vector of column names to group by (e.g., "Sex", "GW_bucket")
#' @param dataset Wide-format dataset to get Formula from
#' @param compound_annotation Optional annotation dataframe to merge
#' @return Dataframe with metabolite abundance by group
calculate_metabolite_abundance_by_group <- function(data_long,
                                                    group_by = "Sex",
                                                    dataset = NULL,
                                                    compound_annotation = NULL) {
  
  # Check required columns
  required_cols <- c("Compound.ID", "Name", "Intensity", group_by)
  missing_cols <- setdiff(required_cols, names(data_long))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  cat("Calculating metabolite abundance by", paste(group_by, collapse = ", "), "...\n")
  
  # Calculate mean intensity by group
  metabolite_abundance_grouped <- data_long %>%
    group_by(across(all_of(c("Compound.ID", "Name", group_by)))) %>%
    summarise(mean_intensity = mean(Intensity, na.rm = TRUE),
              .groups = "drop")
  
  # Add Formula if dataset provided
  if (!is.null(dataset) && "Formula" %in% names(dataset)) {
    metabolite_abundance_grouped <- metabolite_abundance_grouped %>%
      left_join(dataset %>% select(Compound.ID, Formula), by = "Compound.ID")
  }
  
  # Create display name
  metabolite_abundance_grouped <- metabolite_abundance_grouped %>%
    mutate(display_name = ifelse(is.na(Name) | Name == "", 
                                 Compound.ID, 
                                 paste0(Name, " (", Compound.ID, ")")))
  
  # Arrange by group and intensity
  metabolite_abundance_grouped <- metabolite_abundance_grouped %>%
    arrange(across(all_of(group_by)), desc(mean_intensity))
  
  # Add annotations if provided
  if (!is.null(compound_annotation)) {
    metabolite_abundance_grouped <- metabolite_abundance_grouped %>%
      left_join(
        compound_annotation %>% select(Compound.ID, Name, Class, Sub.Class, Super.Class, Pathways),
        by = c("Compound.ID", "Name")
      )
  }
  
  cat("Done!\n\n")
  
  return(metabolite_abundance_grouped)
}

# ============================================================================
# Top Metabolite Visualization Functions
# ============================================================================

#' Plot top N most abundant metabolites as bar chart
#'
#' @param metabolite_abundance Dataframe from calculate_metabolite_abundance
#' @param n Number of top metabolites to plot
#' @param title_prefix Prefix for plot title
#' @param png_path Path to save PNG
#' @param svg_path Path to save SVG
#' @param width Plot width
#' @param height Plot height
#' @param dpi PNG resolution
#' @return ggplot object
plot_top_metabolites_bar <- function(metabolite_abundance,
                                     n = 20,
                                     title_prefix = "",
                                     png_path = NULL,
                                     svg_path = NULL,
                                     width = 10,
                                     height = 8,
                                     dpi = 300) {
  
  top_n_metabolites <- metabolite_abundance %>%
    slice_head(n = n)
  
  p <- ggplot(top_n_metabolites, aes(x = reorder(display_name, mean_intensity), y = mean_intensity)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Top ", n, " Most Abundant Metabolites"),
         x = "Metabolite",
         y = "Mean Peak Intensity") +
    theme(axis.text.y = element_text(size = 8))
  
  save_dual_format(p, paste0(title_prefix, "top", n, "_metabolites"), 
                   png_path, svg_path, width, height, dpi)
  
  return(p)
}

#' Plot top N metabolites stratified by group as faceted bar chart
#'
#' @param metabolite_abundance_grouped Dataframe from calculate_metabolite_abundance_by_group
#' @param group_by Column name to facet by
#' @param n Number of top metabolites per group
#' @param title_prefix Prefix for plot title
#' @param png_path Path to save PNG
#' @param svg_path Path to save SVG
#' @param width Plot width
#' @param height Plot height
#' @param dpi PNG resolution
#' @return ggplot object
plot_top_metabolites_by_group <- function(metabolite_abundance_grouped,
                                          group_by = "Sex",
                                          n = 20,
                                          title_prefix = "",
                                          png_path = NULL,
                                          svg_path = NULL,
                                          width = 14,
                                          height = 10,
                                          dpi = 300) {
  
  top_n_by_group <- metabolite_abundance_grouped %>%
    group_by(across(all_of(group_by))) %>%
    slice_head(n = n) %>%
    ungroup()
  
  # Determine number of columns for faceting
  n_groups <- length(unique(top_n_by_group[[group_by]]))
  ncol <- min(n_groups, 3)
  
  # Adjust text size based on number of groups
  text_size <- ifelse(n_groups <= 2, 7, 6)
  
  p <- ggplot(top_n_by_group, 
              aes(x = reorder(display_name, mean_intensity), 
                  y = mean_intensity, 
                  fill = .data[[group_by]])) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free", ncol = ncol) +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Top ", n, " Most Abundant Metabolites by ", group_by),
         x = "Metabolite",
         y = "Mean Peak Intensity") +
    theme(axis.text.y = element_text(size = text_size))
  
  save_dual_format(p, paste0(title_prefix, "top", n, "_metabolites_by_", tolower(group_by)), 
                   png_path, svg_path, width, height, dpi)
  
  return(p)
}

# ============================================================================
# Heatmap Functions
# ============================================================================

#' Prepare heatmap data matrix from wide-format dataset
#'
#' @param dataset Wide-format dataset
#' @param metabolite_ids Vector of Compound.IDs to include
#' @param sample_cols Column indices or names containing sample data
#' @param log_transform Whether to log10 transform (adds 1 to avoid log(0))
#' @param impute Whether to impute missing values
#' @param impute_method Method for imputation: "min_fraction" or "zero"
#' @param min_fraction Fraction of minimum value to use for imputation (default: 0.2)
#' @return Matrix ready for heatmap plotting
prepare_heatmap_matrix <- function(dataset,
                                   metabolite_ids,
                                   sample_cols = NULL,
                                   sample_pattern = "^Sample",
                                   log_transform = TRUE,
                                   impute = FALSE,
                                   impute_method = "min_fraction",
                                   min_fraction = 0.2) {
  
  # Determine sample columns
  if (is.null(sample_cols)) {
    sample_cols <- grep(sample_pattern, names(dataset), value = TRUE)
  } else if (is.numeric(sample_cols)) {
    sample_cols <- names(dataset)[sample_cols]
  }
  
  # Filter to specified metabolites and prepare matrix
  heatmap_data <- dataset %>%
    filter(Compound.ID %in% metabolite_ids) %>%
    select(Compound.ID, Name, all_of(sample_cols)) %>%
    mutate(display_name = ifelse(is.na(Name) | Name == "", 
                                 Compound.ID, 
                                 paste0(Name, " (", Compound.ID, ")"))) %>%
    select(display_name, all_of(sample_cols)) %>%
    column_to_rownames("display_name") %>%
    as.matrix()
  
  # Log transform if requested
  if (log_transform) {
    heatmap_data <- log10(heatmap_data + 1)
  }
  
  # Handle NaN and Inf values
  heatmap_data[is.nan(heatmap_data)] <- NA
  heatmap_data[is.infinite(heatmap_data)] <- NA
  
  # Impute missing values if requested
  if (impute) {
    if (impute_method == "min_fraction") {
      for (i in 1:nrow(heatmap_data)) {
        row_data <- heatmap_data[i, ]
        if (any(is.na(row_data))) {
          min_value <- min(row_data, na.rm = TRUE)
          heatmap_data[i, is.na(row_data)] <- min_value * min_fraction
        }
      }
    } else if (impute_method == "zero") {
      heatmap_data[is.na(heatmap_data)] <- 0
    }
  }
  
  return(heatmap_data)
}

#' Plot heatmap of metabolite abundances
#'
#' @param heatmap_matrix Matrix from prepare_heatmap_matrix
#' @param sample_metadata Metadata for sample annotation
#' @param annotation_cols Vector of metadata columns to show
#' @param cluster_rows Whether to cluster rows
#' @param cluster_cols Whether to cluster columns
#' @param scale Scaling option: "none", "row", or "column"
#' @param title Plot title
#' @param filename Full path to save file (include extension)
#' @param width Plot width
#' @param height Plot height
#' @param fontsize_row Font size for row labels
#' @return NULL (saves plot to file)
plot_metabolite_heatmap <- function(heatmap_matrix,
                                    sample_metadata = NULL,
                                    annotation_cols = c("Sex", "GW_bucket"),
                                    cluster_rows = TRUE,
                                    cluster_cols = TRUE,
                                    scale = "row",
                                    title = "Metabolite Heatmap",
                                    filename = NULL,
                                    width = 12,
                                    height = 16,
                                    fontsize_row = 6) {
  
  # Prepare sample annotations if metadata provided
  sample_annotation <- NULL
  if (!is.null(sample_metadata)) {
    # Get samples present in matrix
    matrix_samples <- colnames(heatmap_matrix)
    
    # Filter metadata and prepare annotation
    available_annotation_cols <- intersect(annotation_cols, names(sample_metadata))
    
    if (length(available_annotation_cols) > 0) {
      sample_annotation <- sample_metadata %>%
        filter(Sample %in% matrix_samples) %>%
        select(Sample, all_of(available_annotation_cols)) %>%
        column_to_rownames("Sample")
    }
  }
  
  # Set NA color based on whether we have missing data
  na_col <- if (any(is.na(heatmap_matrix))) "gray90" else NULL
  
  # Create heatmap
  pheatmap(heatmap_matrix,
           scale = scale,
           annotation_col = sample_annotation,
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           main = title,
           fontsize_row = fontsize_row,
           filename = filename,
           width = width,
           height = height,
           na_col = na_col)
}

#' Generate all three versions of heatmap (unscaled, scaled, clustered)
#'
#' @param dataset Wide-format dataset
#' @param metabolite_abundance Abundance dataframe to get top metabolites
#' @param n Number of top metabolites to plot
#' @param sample_metadata Metadata for annotations
#' @param title_prefix Prefix for titles and filenames
#' @param png_path Path to save PNG files
#' @param svg_path Path to save SVG files (note: pheatmap only saves PNG)
#' @return List of three heatmap matrices
plot_top_metabolites_heatmaps <- function(dataset,
                                          metabolite_abundance,
                                          n = 100,
                                          sample_metadata = NULL,
                                          sample_cols = NULL,
                                          title_prefix = "",
                                          png_path = NULL,
                                          svg_path = NULL) {
  
  # Get top N metabolites
  top_metabolites <- metabolite_abundance %>%
    slice_head(n = n) %>%
    pull(Compound.ID)
  
  cat("Generating heatmaps for top", n, "metabolites...\n")
  
  # Version 1: Unclustered, unscaled (with NAs preserved)
  cat("  Creating unclustered, unscaled heatmap...\n")
  matrix_unclustered_unscaled <- prepare_heatmap_matrix(
    dataset = dataset,
    metabolite_ids = top_metabolites,
    sample_cols = sample_cols,
    log_transform = TRUE,
    impute = FALSE
  )
  
  if (!is.null(png_path)) {
    plot_metabolite_heatmap(
      heatmap_matrix = matrix_unclustered_unscaled,
      sample_metadata = sample_metadata,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      scale = "none",
      title = paste0(title_prefix, "Top ", n, " Metabolites - Unclustered, Unscaled"),
      filename = file.path(png_path, paste0(title_prefix, "top", n, "_heatmap_unclustered_unscaled.png"))
    )
  }
  
  # Version 2: Unclustered, scaled by row
  cat("  Creating unclustered, scaled heatmap...\n")
  matrix_unclustered_scaled <- matrix_unclustered_unscaled  # Same matrix
  
  if (!is.null(png_path)) {
    plot_metabolite_heatmap(
      heatmap_matrix = matrix_unclustered_scaled,
      sample_metadata = sample_metadata,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      scale = "row",
      title = paste0(title_prefix, "Top ", n, " Metabolites - Unclustered, Scaled"),
      filename = file.path(png_path, paste0(title_prefix, "top", n, "_heatmap_unclustered_scaled.png"))
    )
  }
  
  # Version 3: Clustered, scaled, with imputation
  cat("  Creating clustered, scaled heatmap with imputation...\n")
  matrix_clustered <- prepare_heatmap_matrix(
    dataset = dataset,
    metabolite_ids = top_metabolites,
    sample_cols = sample_cols,
    log_transform = TRUE,
    impute = TRUE,
    impute_method = "min_fraction",
    min_fraction = 0.2
  )
  
  if (!is.null(png_path)) {
    plot_metabolite_heatmap(
      heatmap_matrix = matrix_clustered,
      sample_metadata = sample_metadata,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "row",
      title = paste0(title_prefix, "Top ", n, " Metabolites - Clustered"),
      filename = file.path(png_path, paste0(title_prefix, "top", n, "_heatmap_clustered.png"))
    )
  }
  
  cat("Done!\n\n")
  
  return(list(
    unclustered_unscaled = matrix_unclustered_unscaled,
    unclustered_scaled = matrix_unclustered_scaled,
    clustered = matrix_clustered
  ))
}

# ============================================================================
# Annotation Category Abundance Functions
# ============================================================================

#' Calculate abundance aggregated by annotation category (Class, Sub.Class, Super.Class, etc.)
#'
#' @param metabolite_abundance Dataframe from calculate_metabolite_abundance with annotations
#' @param category Annotation column to aggregate by (e.g., "Class", "Sub.Class", "Super.Class")
#' @return Dataframe with category-level abundance
calculate_category_abundance <- function(metabolite_abundance, 
                                         category = "Class") {
  
  if (!category %in% names(metabolite_abundance)) {
    stop(paste("Category", category, "not found in metabolite_abundance. Did you add compound_annotation?"))
  }
  
  cat("Calculating abundance by", category, "...\n")
  
  # Aggregate by category
  category_abundance <- metabolite_abundance %>%
    filter(!is.na(.data[[category]]) & .data[[category]] != "") %>%
    group_by(.data[[category]]) %>%
    summarise(
      mean_intensity = mean(mean_intensity, na.rm = TRUE),
      n_metabolites = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_intensity))
  
  # Rename first column to standard name for easier downstream use
  names(category_abundance)[1] <- "category_name"
  category_abundance$category_type <- category
  
  cat("  Found", nrow(category_abundance), "unique", category, "categories\n")
  cat("Done!\n\n")
  
  return(category_abundance)
}

#' Calculate category abundance stratified by a grouping variable
#'
#' @param metabolite_abundance_grouped Dataframe from calculate_metabolite_abundance_by_group with annotations
#' @param category Annotation column to aggregate by
#' @param group_by Grouping variable (e.g., "Sex", "GW_bucket")
#' @return Dataframe with category-level abundance by group
calculate_category_abundance_by_group <- function(metabolite_abundance_grouped,
                                                  category = "Class",
                                                  group_by = "Sex") {
  
  if (!category %in% names(metabolite_abundance_grouped)) {
    stop(paste("Category", category, "not found in data"))
  }
  
  if (!group_by %in% names(metabolite_abundance_grouped)) {
    stop(paste("Group variable", group_by, "not found in data"))
  }
  
  cat("Calculating", category, "abundance by", group_by, "...\n")
  
  # Aggregate by category and group
  category_abundance_grouped <- metabolite_abundance_grouped %>%
    filter(!is.na(.data[[category]]) & .data[[category]] != "") %>%
    group_by(.data[[category]], .data[[group_by]]) %>%
    summarise(
      mean_intensity = mean(mean_intensity, na.rm = TRUE),
      n_metabolites = n(),
      .groups = "drop"
    ) %>%
    arrange(.data[[group_by]], desc(mean_intensity))
  
  # Rename category column
  names(category_abundance_grouped)[1] <- "category_name"
  category_abundance_grouped$category_type <- category
  
  cat("Done!\n\n")
  
  return(category_abundance_grouped)
}

# ============================================================================
# Category Visualization Functions
# ============================================================================

#' Plot top N categories as bar chart
#'
#' @param category_abundance Dataframe from calculate_category_abundance
#' @param n Number of top categories to plot
#' @param title_prefix Prefix for plot title
#' @param png_path Path to save PNG
#' @param svg_path Path to save SVG
#' @param width Plot width
#' @param height Plot height
#' @param dpi PNG resolution
#' @return ggplot object
plot_top_categories_bar <- function(category_abundance,
                                    n = 20,
                                    title_prefix = "",
                                    png_path = NULL,
                                    svg_path = NULL,
                                    width = 10,
                                    height = 8,
                                    dpi = 300) {
  
  # Get category type for labeling
  category_type <- unique(category_abundance$category_type)[1]
  
  top_n_categories <- category_abundance %>%
    slice_head(n = n)
  
  p <- ggplot(top_n_categories, 
              aes(x = reorder(category_name, mean_intensity), y = mean_intensity)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Top ", n, " Most Abundant Metabolite ", category_type, "es"),
         x = category_type,
         y = "Mean Peak Intensity") +
    theme(axis.text.y = element_text(size = 8))
  
  save_dual_format(p, 
                   paste0(title_prefix, "top", n, "_", tolower(category_type)), 
                   png_path, svg_path, width, height, dpi)
  
  return(p)
}

#' Plot top N categories stratified by group
#'
#' @param category_abundance_grouped Dataframe from calculate_category_abundance_by_group
#' @param group_by Column name to facet by
#' @param n Number of top categories per group
#' @param title_prefix Prefix for plot title
#' @param png_path Path to save PNG
#' @param svg_path Path to save SVG
#' @param width Plot width
#' @param height Plot height
#' @param dpi PNG resolution
#' @return ggplot object
plot_top_categories_by_group <- function(category_abundance_grouped,
                                         group_by = "Sex",
                                         n = 20,
                                         title_prefix = "",
                                         png_path = NULL,
                                         svg_path = NULL,
                                         width = 14,
                                         height = 10,
                                         dpi = 300) {
  
  # Get category type for labeling
  category_type <- unique(category_abundance_grouped$category_type)[1]
  
  top_n_by_group <- category_abundance_grouped %>%
    group_by(.data[[group_by]]) %>%
    slice_head(n = n) %>%
    ungroup()
  
  # Determine faceting parameters
  n_groups <- length(unique(top_n_by_group[[group_by]]))
  ncol <- min(n_groups, 3)
  text_size <- ifelse(n_groups <= 2, 7, 6)
  
  p <- ggplot(top_n_by_group, 
              aes(x = reorder(category_name, mean_intensity), 
                  y = mean_intensity, 
                  fill = .data[[group_by]])) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free", ncol = ncol) +
    theme_minimal() +
    labs(title = paste0(title_prefix, "Top ", n, " Most Abundant Metabolite ", 
                        category_type, "es by ", group_by),
         x = category_type,
         y = "Mean Peak Intensity") +
    theme(axis.text.y = element_text(size = text_size))
  
  save_dual_format(p, 
                   paste0(title_prefix, "top", n, "_", tolower(category_type), "_by_", tolower(group_by)), 
                   png_path, svg_path, width, height, dpi)
  
  return(p)
}

# ============================================================================
# Category Heatmap Functions
# ============================================================================

#' Prepare category-aggregated heatmap matrix
#'
#' @param dataset Wide-format dataset with metabolite rows and sample columns
#' @param metabolite_abundance_annotated Abundance data with annotations
#' @param category_ids Vector of category names to include
#' @param category Category column name (e.g., "Class", "Sub.Class")
#' @param sample_cols Column indices or names containing sample data
#' @param log_transform Whether to log10 transform
#' @param impute Whether to impute missing values
#' @param min_fraction Fraction of minimum value for imputation
#' @return Matrix ready for heatmap
prepare_category_heatmap_matrix <- function(dataset,
                                            metabolite_abundance_annotated,
                                            category_ids,
                                            category = "Class",
                                            sample_cols = NULL,
                                            sample_pattern = "^Sample",
                                            log_transform = TRUE,
                                            impute = TRUE,
                                            min_fraction = 0.2) {
  
  # Determine sample columns
  if (is.null(sample_cols)) {
    sample_cols <- grep(sample_pattern, names(dataset), value = TRUE)
  } else if (is.numeric(sample_cols)) {
    sample_cols <- names(dataset)[sample_cols]
  }
  
  # Get metabolites in the specified categories
  metabolites_in_categories <- metabolite_abundance_annotated %>%
    filter(.data[[category]] %in% category_ids) %>%
    select(Compound.ID, !!sym(category))
  
  # Aggregate intensities by category and sample
  category_heatmap_data <- metabolite_abundance_annotated %>%
    filter(.data[[category]] %in% category_ids) %>%
    left_join(dataset %>% select(Compound.ID, all_of(sample_cols)), by = "Compound.ID") %>%
    select(!!sym(category), all_of(sample_cols)) %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Intensity") %>%
    group_by(.data[[category]], Sample) %>%
    summarise(mean_intensity = mean(Intensity, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Sample, values_from = mean_intensity) %>%
    column_to_rownames(category) %>%
    as.matrix()
  
  # Handle NaN values
  category_heatmap_data[is.nan(category_heatmap_data)] <- NA
  
  # Log transform if requested
  if (log_transform) {
    category_heatmap_data <- log10(category_heatmap_data + 1)
  }
  
  # Handle Inf values
  category_heatmap_data[is.infinite(category_heatmap_data)] <- NA
  
  # Impute missing values if requested
  if (impute) {
    for (i in 1:nrow(category_heatmap_data)) {
      row_data <- category_heatmap_data[i, ]
      if (any(is.na(row_data))) {
        min_value <- min(row_data, na.rm = TRUE)
        if (is.finite(min_value)) {
          category_heatmap_data[i, is.na(row_data)] <- min_value * min_fraction
        }
      }
    }
  }
  
  return(category_heatmap_data)
}

#' Plot heatmap of top N categories
#'
#' @param dataset Wide-format dataset
#' @param metabolite_abundance_annotated Abundance data with annotations
#' @param category_abundance Category abundance dataframe
#' @param category Category type ("Class", "Sub.Class", etc.)
#' @param n Number of top categories to plot
#' @param sample_metadata Metadata for annotations
#' @param title_prefix Prefix for title and filename
#' @param png_path Path to save PNG
#' @return Matrix used for heatmap
plot_top_categories_heatmap <- function(dataset,
                                        metabolite_abundance_annotated,
                                        category_abundance,
                                        category = "Class",
                                        n = 100,
                                        sample_metadata = NULL,
                                        sample_cols = NULL,
                                        title_prefix = "",
                                        png_path = NULL) {
  
  # Get top N categories
  top_categories <- category_abundance %>%
    slice_head(n = n) %>%
    pull(category_name)
  
  cat("Generating heatmap for top", n, category, "categories...\n")
  
  # Prepare matrix
  category_matrix <- prepare_category_heatmap_matrix(
    dataset = dataset,
    metabolite_abundance_annotated = metabolite_abundance_annotated,
    category_ids = top_categories,
    category = category,
    sample_cols = sample_cols,
    log_transform = TRUE,
    impute = TRUE,
    min_fraction = 0.2
  )
  
  # Plot heatmap
  if (!is.null(png_path)) {
    plot_metabolite_heatmap(
      heatmap_matrix = category_matrix,
      sample_metadata = sample_metadata,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "row",
      title = paste0(title_prefix, "Top ", n, " Most Abundant Metabolite ", category, "es"),
      filename = file.path(png_path, paste0(title_prefix, "top", n, "_", tolower(category), "_heatmap.png")),
      fontsize_row = 6
    )
  }
  
  cat("Done!\n\n")
  
  return(category_matrix)
}

# ============================================================================
# Wrapper Functions for Complete Analysis
# ============================================================================

#' Generate all visualizations for a specific annotation category
#'
#' @param metabolite_abundance_annotated Overall abundance with annotations
#' @param metabolite_abundance_by_sex Sex-stratified abundance with annotations
#' @param metabolite_abundance_by_gw GW-stratified abundance with annotations
#' @param dataset Wide-format dataset
#' @param category Annotation category to analyze
#' @param n_bar Number of categories for bar charts
#' @param n_heatmap Number of categories for heatmap
#' @param sample_metadata Sample metadata
#' @param title_prefix Title prefix
#' @param png_path PNG save path
#' @param svg_path SVG save path
#' @return List of all plots and data
analyze_category_abundance <- function(metabolite_abundance_annotated,
                                       metabolite_abundance_by_sex,
                                       metabolite_abundance_by_gw,
                                       dataset,
                                       category = "Class",
                                       n_bar = 20,
                                       n_heatmap = 100,
                                       sample_metadata = NULL,
                                       sample_cols = NULL,
                                       title_prefix = "",
                                       png_path = NULL,
                                       svg_path = NULL) {
  
  cat("\n=== ANALYZING", toupper(category), "ABUNDANCE ===\n\n")
  
  # Calculate abundances
  cat("Step 1: Calculating abundance metrics...\n")
  category_abundance <- calculate_category_abundance(
    metabolite_abundance = metabolite_abundance_annotated,
    category = category
  )
  
  category_abundance_by_sex <- calculate_category_abundance_by_group(
    metabolite_abundance_grouped = metabolite_abundance_by_sex,
    category = category,
    group_by = "Sex"
  )
  
  category_abundance_by_gw <- calculate_category_abundance_by_group(
    metabolite_abundance_grouped = metabolite_abundance_by_gw,
    category = category,
    group_by = "GW_bucket"
  )
  
  # Generate visualizations
  cat("Step 2: Generating bar charts...\n")
  plot_overall <- plot_top_categories_bar(
    category_abundance = category_abundance,
    n = n_bar,
    title_prefix = title_prefix,
    png_path = png_path,
    svg_path = svg_path
  )
  
  plot_by_sex <- plot_top_categories_by_group(
    category_abundance_grouped = category_abundance_by_sex,
    group_by = "Sex",
    n = n_bar,
    title_prefix = title_prefix,
    png_path = png_path,
    svg_path = svg_path
  )
  
  plot_by_gw <- plot_top_categories_by_group(
    category_abundance_grouped = category_abundance_by_gw,
    group_by = "GW_bucket",
    n = n_bar,
    title_prefix = title_prefix,
    png_path = png_path,
    svg_path = svg_path,
    width = 18
  )
  
  cat("Step 3: Generating heatmap...\n")
  heatmap_matrix <- plot_top_categories_heatmap(
    dataset = dataset,
    metabolite_abundance_annotated = metabolite_abundance_annotated,
    category_abundance = category_abundance,
    category = category,
    n = n_heatmap,
    sample_metadata = sample_metadata,
    sample_cols = sample_cols,
    title_prefix = title_prefix,
    png_path = png_path
  )
  
  cat("=== COMPLETED", toupper(category), "ANALYSIS ===\n\n")
  
  return(list(
    abundance = category_abundance,
    abundance_by_sex = category_abundance_by_sex,
    abundance_by_gw = category_abundance_by_gw,
    plots = list(
      overall = plot_overall,
      by_sex = plot_by_sex,
      by_gw = plot_by_gw
    ),
    heatmap_matrix = heatmap_matrix
  ))
}

