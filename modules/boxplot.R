# Boxplot Generation for Streamlit App
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Source common functions - using the correct file name
source("viz_file_handling.R")

# Custom color palette - EXACTLY as in the standalone script
CUSTOM_COLORS <- c("#FF3200", "#73D2CD", "#890620", "#DCCCA3", "#762A83", "#E9A800")

# Function to read count matrix for boxplots
read_count_matrix_boxplot <- function(file_path, type) {
  if (tools::file_ext(file_path) == "csv") {
    count_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    count_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  return(count_data)
}

# Function to calculate text sizes based on dimensions - IMPROVED VERSION
calculate_text_sizes <- function(width, height) {
  # Calculate the area of the image in megapixels
  area_mp <- (width * height) / 1000000
  
  # Base scaling factors (more responsive to image size)
  # For a "standard" 3000x2600 image (~7.8 MP), we want reasonable base sizes
  base_area <- (1200 * 800) / 1000000  # ~7.8 MP
  
  # Calculate scaling factor based on area ratio with logarithmic scaling
  # This prevents text from becoming too large on huge images
  area_ratio <- area_mp / base_area
  scale_factor <- 1 + log(area_ratio) * 0.3  # More gradual scaling
  
  # Base text sizes
  base_title_size <- 9
  base_axis_text_size <- 8
  base_axis_title_size <- 9
  base_legend_text_size <- 8
  
  # Calculate scaled text sizes with reasonable limits
  title_size <- max(6, min(36, base_title_size * scale_factor))
  axis_text_size <- max(4, min(24, base_axis_text_size * scale_factor))
  axis_title_size <- max(5, min(28, base_axis_title_size * scale_factor))
  legend_text_size <- max(4, min(20, base_legend_text_size * scale_factor))
  
  cat("Area scaling factor:", scale_factor, "\n")
  
  return(list(
    title_size = title_size,
    axis_text_size = axis_text_size,
    axis_title_size = axis_title_size,
    legend_text_size = legend_text_size
  ))
}

# Function to create boxplots
create_boxplots_streamlit <- function(args) {
  # Parse arguments
  count_file <- args[1]
  metadata_file <- args[2]
  output_dir <- args[3]
  plot_type <- args[4]  # "gene" or "transcript"
  genes_input <- args[5]
  width <- as.numeric(args[6])
  height <- as.numeric(args[7])
  dpi <- as.numeric(args[8])  # Now hardcoded to 300
  
  cat("Starting boxplot generation...\n")
  cat("Dimensions:", width, "x", height, "pixels\n")
  
  # Calculate text sizes based on dimensions
  text_sizes <- calculate_text_sizes(width, height)
  cat("Calculated text sizes:\n")
  cat("  Title:", text_sizes$title_size, "\n")
  cat("  Axis text:", text_sizes$axis_text_size, "\n")
  cat("  Axis title:", text_sizes$axis_title_size, "\n")
  cat("  Legend text:", text_sizes$legend_text_size, "\n")
  
  # Read genes from input or file
  if (file.exists(genes_input)) {
    genes_to_plot <- readLines(genes_input, warn = FALSE)
    genes_to_plot <- genes_to_plot[genes_to_plot != ""]  # Remove empty lines
    cat("Read", length(genes_to_plot), "genes from file\n")
  } else {
    genes_to_plot <- unlist(strsplit(genes_input, ","))
    cat("Read", length(genes_to_plot), "genes from input\n")
  }
  
  # Clean gene names
  genes_to_plot <- trimws(genes_to_plot)
  
  # Read count matrix
  cat("Reading count matrix...\n")
  count_data <- read_count_matrix_boxplot(count_file, plot_type)
  cat("Count data dimensions:", dim(count_data), "\n")
  
  # Read metadata
  cat("Reading metadata...\n")
  metadata <- read_metadata_common(metadata_file)
  col_info <- detect_metadata_columns(metadata)
  cat("Detected sample column:", col_info$sample_col, "\n")
  cat("Detected condition column:", col_info$condition_col, "\n")
  
  # Check if required columns exist
  if (is.null(col_info$sample_col)) {
    stop("Could not detect sample column in metadata. Please ensure metadata has a sample identifier column.")
  }
  
  if (is.null(col_info$condition_col)) {
    stop("Could not detect condition column in metadata. Please ensure metadata has a condition column.")
  }
  
  # Prepare count data based on type
  if (plot_type == "gene") {
    if (!"gene_symbol" %in% colnames(count_data)) {
      stop("Column 'gene_symbol' not found in gene count data.")
    }
    id_col <- "gene_symbol"
  } else {
    if (!"transcript_name" %in% colnames(count_data)) {
      stop("Column 'transcript_name' not found in transcript count data.")
    }
    id_col <- "transcript_name"
  }
  
  # Get sample columns (all columns except ID columns)
  id_columns <- c()
  if ("gene_id" %in% colnames(count_data)) id_columns <- c(id_columns, "gene_id")
  if ("transcript_id" %in% colnames(count_data)) id_columns <- c(id_columns, "transcript_id")
  id_columns <- c(id_columns, id_col)
  
  sample_columns <- setdiff(colnames(count_data), id_columns)
  cat("Sample columns found:", length(sample_columns), "\n")
  
  # Reshape data to long format
  cat("Reshaping data to long format...\n")
  long_data <- count_data %>%
    select(all_of(c(id_col, sample_columns))) %>%
    pivot_longer(
      cols = -all_of(id_col),
      names_to = "Sample",
      values_to = "Expression"
    ) %>%
    rename(Gene_IDs = all_of(id_col))
  
  # Merge with metadata
  cat("Merging with metadata...\n")
  long_data <- long_data %>%
    left_join(metadata, by = c("Sample" = col_info$sample_col)) %>%
    filter(!is.na(!!sym(col_info$condition_col)))
  
  # Rename condition column for consistency
  long_data <- long_data %>%
    rename(condition = !!sym(col_info$condition_col))
  
  # Filter for specified genes
  long_data <- long_data %>% filter(Gene_IDs %in% genes_to_plot)
  cat("Filtered data rows:", nrow(long_data), "\n")
  
  if (nrow(long_data) == 0) {
    stop("None of the specified genes were found in the count data.")
  }
  
  # Get unique genes and conditions
  unique_genes <- unique(long_data$Gene_IDs)
  unique_conditions <- unique(long_data$condition)
  cat("Unique genes:", length(unique_genes), "\n")
  cat("Unique conditions:", length(unique_conditions), "\n")
  
  # EXACT SAME COLOR LOGIC AS STANDALONE SCRIPT
  # Create a color palette for conditions - use custom colors if 6 or fewer conditions
  n_conditions <- length(unique_conditions)
  if (n_conditions <= 6) {
    color_palette <- CUSTOM_COLORS[1:n_conditions]
    names(color_palette) <- unique_conditions
  } else {
    # Use default palette for more than 6 conditions
    color_palette <- scales::hue_pal()(n_conditions)
    names(color_palette) <- unique_conditions
    cat("Note: More than 6 conditions detected. Using default color palette.\n")
  }
  
  # Create plots for each gene
  plot_files <- c()
  for (gene in unique_genes) {
    gene_data <- long_data %>% filter(Gene_IDs == gene)
    
    if (nrow(gene_data) > 0) {
      # Create plot with dynamic text sizing
      gene_plot <- ggplot(gene_data, 
                          aes(x = condition, y = Expression, fill = condition)) + 
        geom_boxplot(alpha = 1) + 
        scale_fill_manual(values = color_palette) +
        xlab("Condition") + 
        ylab("Expression") + 
        ggtitle(paste(gene, "Expression Across Conditions")) +
        theme_bw() +
        theme(
          panel.border = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5),
          # Gridlines
          panel.grid.major = element_line(colour = "gray60", linewidth = 0.2),
          # Dynamic text sizing
          plot.title = element_text(face = "bold", size = text_sizes$title_size, hjust = 0.5),
          axis.text = element_text(size = text_sizes$axis_text_size),
          axis.title = element_text(size = text_sizes$axis_title_size),
          legend.text = element_text(size = text_sizes$legend_text_size),
          legend.title = element_text(size = text_sizes$legend_text_size)
        )
      
      # Save plot - using hardcoded DPI of 300
      filename <- file.path(output_dir, paste0(gene, "_", plot_type, "_boxplot.png"))
      ggsave(filename, plot = gene_plot, width = width/300, height = height/300, dpi = dpi, bg = "white")
      plot_files <- c(plot_files, filename)
      cat("Created plot for", gene, "-", plot_type, "\n")
    }
  }
  
  cat("Boxplot generation completed. Created", length(plot_files), "plots.\n")
  return(plot_files)
}

# Run if called from command line
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 8) {
    tryCatch({
      create_boxplots_streamlit(args)
    }, error = function(e) {
      cat("Error in boxplot generation:", e$message, "\n")
      quit(status = 1)
    })
  } else {
    cat("Usage: Rscript boxplot.R <count_file> <metadata_file> <output_dir> <plot_type> <genes> <width> <height> <dpi>\n")
    quit(status = 1)
  }
}