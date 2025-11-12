# Heatmap Generation for Streamlit App
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(tools)
  library(paletteer)
})

# Source common functions - using the correct file name
source("viz_file_handling.R")

# Function to read count matrix
read_count_matrix_streamlit <- function(file_path, type) {
  cat("Reading count matrix from:", file_path, "\n")
  
  if (tools::file_ext(file_path) == "csv") {
    count_data <- tryCatch({
      read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      stop(paste("Failed to read CSV file:", e$message))
    })
  } else {
    count_data <- tryCatch({
      read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      stop(paste("Failed to read TSV file:", e$message))
    })
  }
  
  # Clean column names
  colnames(count_data) <- gsub("^\ufeff", "", colnames(count_data))
  colnames(count_data) <- make.names(colnames(count_data))
  
  cat("Count matrix columns:", ncol(count_data), "\n")
  cat("Count matrix rows:", nrow(count_data), "\n")
  
  if (type == "gene") {
    if (ncol(count_data) < 3) {
      stop("Gene count matrix must have at least 3 columns: gene_id, gene_symbol, and sample columns")
    }
    # Assume first column is gene_id, second is gene_symbol
    gene_symbols <- count_data[, 2]
    gene_symbols <- make.names(gene_symbols)
    count_matrix <- as.matrix(count_data[, 3:ncol(count_data)])
    rownames(count_matrix) <- gene_symbols
  } else if (type == "transcript") {
    if (ncol(count_data) < 3) {
      stop("Transcript count matrix must have at least 3 columns: transcript_id, transcript_name, and sample columns")
    }
    # Assume first column is transcript_id, second is transcript_name
    transcript_names <- count_data[, 2]
    transcript_names <- make.names(transcript_names)
    count_matrix <- as.matrix(count_data[, 3:ncol(count_data)])
    rownames(count_matrix) <- transcript_names
  }
  
  # Clean sample names
  colnames(count_matrix) <- make.names(colnames(count_matrix))
  
  cat("Final matrix dimensions:", dim(count_matrix), "\n")
  return(count_matrix)
}

# Function to get condition colors (Snippet 01)
get_condition_colors <- function(unique_vals) {
  # Fixed condition colors
  condition_colors <- c("#F46D43", "#66C2A5", "#FEE08B", "#039BE5", "#B2FF2E", "#64243E")
  
  n_conditions <- length(unique_vals)
  
  if (n_conditions <= 6) {
    # Use the fixed colors for up to 6 conditions
    colors <- condition_colors[1:n_conditions]
  } else {
    # For more than 6 conditions, create a gradient from the first 4 colors
    colors <- colorRampPalette(condition_colors[1:4])(n_conditions)
  }
  
  names(colors) <- unique_vals
  return(colors)
}

# Function to get metadata column colors (Snippet 01)
get_metadata_colors <- function(unique_vals, col_index) {
  # Define the color gradients for metadata columns
  metadata_gradients <- list(
    c("#41476B", "#C67B6F", "#FBDFA2"),  # Column 1
    c("#EB4A40", "#448C8A", "#6C2167"),  # Column 2
    c("#440154", "#22A884", "#FDE725"),  # Column 3
    c("#5E1521", "#02734A", "#03A62C"),  # Column 4
    c("#390E40", "#A63228", "#EBA42B")   # Column 5
  )
  
  n_categories <- length(unique_vals)
  
  # Cycle through the gradients if we have more than 5 columns
  gradient_index <- ((col_index - 1) %% length(metadata_gradients)) + 1
  base_gradient <- metadata_gradients[[gradient_index]]
  
  if (n_categories <= 3) {
    # For small number of categories, use the gradient colors directly
    colors <- base_gradient[1:n_categories]
  } else {
    # For more categories, create a gradient
    colors <- colorRampPalette(base_gradient)(n_categories)
  }
  
  names(colors) <- unique_vals
  return(colors)
}

# Function to get gene annotation colors (Snippet 01)
get_gene_annotation_colors <- function(unique_vals, col_index) {
  # Define the color gradients for gene annotation columns
  gene_annotation_gradients <- list(
    c("#FCDE9C", "#E34F6F", "#7C1D6F"),  # Column 1
    c("#EC702E", "#5E1521", "#350E16"),  # Column 2
    c("#CDC597", "#72972F", "#113719"),  # Column 3
    c("#7C1D6F", "#E34F6F", "#FCDE9C")   # Column 4
  )
  
  n_categories <- length(unique_vals)
  
  # Cycle through the gradients if we have more than 4 columns
  gradient_index <- ((col_index - 1) %% length(gene_annotation_gradients)) + 1
  base_gradient <- gene_annotation_gradients[[gradient_index]]
  
  if (n_categories <= 3) {
    # For small number of categories, use the gradient colors directly
    colors <- base_gradient[1:n_categories]
  } else {
    # For more categories, create a gradient
    colors <- colorRampPalette(base_gradient)(n_categories)
  }
  
  names(colors) <- unique_vals
  return(colors)
}

# Function to generate annotation colors (Snippet 01)
generate_annotation_colors <- function(annotation_df, is_row_annotation = FALSE, condition_column = NULL) {
  ann_colors <- list()
  
  # Track column indices for metadata and gene annotations
  metadata_col_index <- 0
  gene_annotation_col_index <- 0
  
  for (col_name in colnames(annotation_df)) {
    if (is.numeric(annotation_df[[col_name]])) {
      # For numeric variables, use color gradient
      ann_colors[[col_name]] <- colorRamp2(
        seq(min(annotation_df[[col_name]], na.rm = TRUE), 
            max(annotation_df[[col_name]], na.rm = TRUE), length = 5),
        viridis::viridis(5)
      )
    } else {
      # For categorical variables
      unique_vals <- unique(annotation_df[[col_name]])
      unique_vals <- unique_vals[!is.na(unique_vals)]  # Remove NA values
      n_categories <- length(unique_vals)
      
      if (n_categories == 0) {
        # Skip if no valid categories
        next
      }
      
      if (!is_row_annotation && !is.null(condition_column) && col_name == condition_column) {
        # Use condition colors for the condition column
        ann_colors[[col_name]] <- get_condition_colors(unique_vals)
      } else if (is_row_annotation) {
        # Use gene annotation colors for row annotations
        gene_annotation_col_index <- gene_annotation_col_index + 1
        ann_colors[[col_name]] <- get_gene_annotation_colors(unique_vals, gene_annotation_col_index)
      } else {
        # Use metadata colors for other column annotations
        metadata_col_index <- metadata_col_index + 1
        ann_colors[[col_name]] <- get_metadata_colors(unique_vals, metadata_col_index)
      }
    }
  }
  
  return(ann_colors)
}

# Define color palettes (Snippet 02)
get_color_palettes <- function() {
  palettes <- list(
    "Green_Teal_Purple" = function(n) {
      color_stops <- c("#D3DDBD", "#579A8D", "#2C1F3F")
      colorRampPalette(color_stops)(n)
    },
    "blue_teal" = function(n) {
      paletteer::paletteer_c("ggthemes::Blue-Teal", n)
    },
    "orange_blue_white" = function(n) {
      # Reversed version of Orange-Blue-White Diverging
      paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", n)
    },
    "red_blue_diverging" = function(n) {
      # Reversed version of Red-Blue Diverging
      paletteer::paletteer_c("ggthemes::Red-Blue Diverging", n)
    },
    "warmcool" = function(n) {
      # Reversed version of warmcool
      paletteer::paletteer_c("pals::warmcool", n)
    }
  )
  
  # Reverse the color schemes for orange_blue_white, red_blue_diverging, and warmcool
  palettes[["orange_blue_white"]] <- function(n) {
    rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", n))
  }
  
  palettes[["red_blue_diverging"]] <- function(n) {
    rev(paletteer::paletteer_c("ggthemes::Red-Blue Diverging", n))
  }
  
  palettes[["warmcool"]] <- function(n) {
    rev(paletteer::paletteer_c("pals::warmcool", n))
  }
  
  return(palettes)
}

# Function to get color palette names for display (Snippet 02)
get_palette_names <- function() {
  c(
    "Green_Teal_Purple",
    "blue_teal", 
    "orange_blue_white",
    "red_blue_diverging",
    "warmcool"
  )
}

# Function to calculate optimal heatmap dimensions
calculate_heatmap_dimensions <- function(n_rows, n_cols, annotation_cols, row_annotation_cols = 0) {
  # Base dimensions
  base_width <- 3000
  base_height <- 2000
  
  # Adjust width based on number of columns
  width_factor <- max(1, n_cols / 20)
  adjusted_width <- base_width * width_factor
  
  # Adjust height based on number of rows
  height_factor <- max(1, n_rows / 30)
  adjusted_height <- base_height * height_factor
  
  # Ensure minimum dimensions
  min_width <- 800
  min_height <- 600
  max_width <- 8000
  max_height <- 6000
  
  final_width <- max(min_width, min(adjusted_width, max_width))
  final_height <- max(min_height, min(adjusted_height, max_height))
  
  # Adjust for annotation columns
  if (length(annotation_cols) > 3) {
    final_width <- final_width * 1.2
  }
  
  # Adjust for row annotations
  if (row_annotation_cols > 0) {
    final_width <- final_width * 1.3
  }
  
  return(list(width = final_width, height = final_height))
}

# Main heatmap creation function for Streamlit
create_heatmap_streamlit <- function(args) {
  tryCatch({
    # Parse arguments
    count_file <- args[1]
    metadata_file <- args[2]
    output_dir <- args[3]
    plot_type <- args[4]  # "gene" or "transcript"
    genes_input <- args[5]
    width <- as.numeric(args[6])
    height <- as.numeric(args[7])
    dpi <- as.numeric(args[8])
    color_palette_name <- if (length(args) > 8) args[9] else "Green_Teal_Purple"
    gene_annotation_file <- if (length(args) > 9 && args[10] != "None") args[10] else NULL
    
    cat("Starting heatmap generation...\n")
    cat("Color palette:", color_palette_name, "\n")
    
    # Read genes from input or file
    if (file.exists(genes_input)) {
      genes_to_plot <- readLines(genes_input, warn = FALSE)
      genes_to_plot <- genes_to_plot[genes_to_plot != ""]  # Remove empty lines
      cat("Read", length(genes_to_plot), "genes from file\n")
    } else {
      stop("Genes list file not found:", genes_input)
    }
    
    if (length(genes_to_plot) == 0) {
      stop("No genes found in the genes list")
    }
    
    # Clean gene names
    genes_to_plot <- make.names(genes_to_plot)
    
    # Read count matrix
    cat("Reading count matrix...\n")
    count_matrix <- read_count_matrix_streamlit(count_file, plot_type)
    
    # Read metadata
    cat("Reading metadata...\n")
    metadata <- read_metadata_common(metadata_file)
    col_info <- detect_metadata_columns(metadata)
    cat("Detected sample column:", col_info$sample_col, "\n")
    cat("Detected condition column:", col_info$condition_col, "\n")
    
    # Use first column as row names if sample column exists
    if (!is.null(col_info$sample_col) && col_info$sample_col %in% colnames(metadata)) {
      rownames(metadata) <- metadata[[col_info$sample_col]]
    } else {
      # Use first column as row names
      rownames(metadata) <- metadata[, 1]
    }
    
    # Ensure sample names match
    common_samples <- intersect(colnames(count_matrix), rownames(metadata))
    cat("Common samples found:", length(common_samples), "\n")
    
    if (length(common_samples) == 0) {
      stop("No common samples found between count matrix and metadata.")
    }
    
    count_matrix <- count_matrix[, common_samples, drop = FALSE]
    metadata <- metadata[common_samples, , drop = FALSE]
    
    # Read gene annotation if provided
    gene_annotation <- NULL
    if (!is.null(gene_annotation_file) && file.exists(gene_annotation_file)) {
      cat("Reading gene annotation file...\n")
      gene_annotation <- read_gene_annotation_common(gene_annotation_file, genes_to_plot)
    }
    
    # Subset count matrix for requested features
    available_features <- genes_to_plot[genes_to_plot %in% rownames(count_matrix)]
    cat("Available features in count matrix:", length(available_features), "\n")
    
    if (length(available_features) == 0) {
      stop("None of the requested features were found in the count matrix.")
    }
    
    subset_matrix <- count_matrix[available_features, , drop = FALSE]
    
    # Scale the matrix
    cat("Scaling matrix...\n")
    scaled_matrix <- t(scale(t(subset_matrix)))
    scaled_matrix[is.na(scaled_matrix)] <- 0
    scaled_matrix[is.nan(scaled_matrix)] <- 0
    scaled_matrix[is.infinite(scaled_matrix)] <- 0
    
    # Calculate optimal dimensions
    row_annotation_count <- if (!is.null(gene_annotation)) ncol(gene_annotation) else 0
    dims <- calculate_heatmap_dimensions(nrow(subset_matrix), ncol(subset_matrix), colnames(metadata), row_annotation_count)
    
    # Override with user dimensions if provided
    if (width > 0 && height > 0) {
      dims$width <- width
      dims$height <- height
    }
    
    cat("Using dimensions:", dims$width, "x", dims$height, "pixels\n")
    
    # Generate annotation colors
    ann_colors <- generate_annotation_colors(metadata, is_row_annotation = FALSE, condition_column = col_info$condition_col)
    
    # Generate row annotation colors if row annotations are provided
    row_ann_colors <- NULL
    row_annotation <- NULL
    
    if (!is.null(gene_annotation)) {
      # Ensure row annotations match the order of features in the matrix
      gene_annotation <- gene_annotation[available_features, , drop = FALSE]
      row_ann_colors <- generate_annotation_colors(gene_annotation, is_row_annotation = TRUE)
      
      # Create row annotation if we have valid colors
      if (length(row_ann_colors) > 0) {
        row_annotation <- rowAnnotation(
          df = gene_annotation,
          col = row_ann_colors,
          annotation_name_gp = gpar(fontsize = 10),
          simple_anno_size = unit(0.4, "cm"),
          # Add borders to row annotation cells (Snippet 03)
          gp = gpar(col = "black", lwd = 0.8)
        )
      } else {
        cat("Warning: No valid row annotation colors generated. Skipping row annotations.\n")
      }
    }
    
    # Create output filename
    output_file <- file.path(output_dir, paste0("heatmap_", plot_type, ".png"))
    cat("Output file:", output_file, "\n")
    
    # Create heatmap
    cat("Creating heatmap...\n")
    png(output_file, width = dims$width, height = dims$height, res = dpi, bg = "white")
    
    # Create color gradient based on selected palette (Snippet 04)
    n_color_steps <- 100
    palettes <- get_color_palettes()
    
    if (color_palette_name %in% names(palettes)) {
      custom_colors <- palettes[[color_palette_name]](n_color_steps)
    } else {
      # Default to Green_Teal_Purple if palette not found
      cat("Warning: Color palette not found. Using default Green-Teal-Purple palette.\n")
      custom_colors <- palettes[["Green_Teal_Purple"]](n_color_steps)
    }
    
    # Determine the data range for color mapping
    data_range <- range(scaled_matrix, na.rm = TRUE)
    cat("Data range for color mapping:", data_range[1], "to", data_range[2], "\n")
    
    # Create breaks that match the number of colors
    color_breaks <- seq(data_range[1], data_range[2], length.out = n_color_steps)
    
    # Calculate optimal font sizes
    base_font_size <- 10
    row_font_size <- max(6, base_font_size - max(0, (nrow(subset_matrix) - 20) / 10))
    legend_font_size <- max(8, base_font_size - max(0, (length(ann_colors) - 3) / 2))
    
    # Create the top annotation with borders (Snippet 03)
    top_annotation <- HeatmapAnnotation(
      df = metadata, 
      col = ann_colors,
      annotation_name_gp = gpar(fontsize = legend_font_size),
      simple_anno_size = unit(0.4, "cm"),
      # Add borders to annotation cells (Snippet 03)
      gp = gpar(col = "black", lwd = 0.8)
    )
    
    # Create heatmap with borders and custom colors
    ht <- Heatmap(
      scaled_matrix,
      name = "Z-score",
      top_annotation = top_annotation,
      left_annotation = row_annotation,
      show_row_names = TRUE,
      show_column_names = FALSE,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      row_names_gp = gpar(fontsize = row_font_size),
      column_names_gp = gpar(fontsize = 10),
      col = colorRamp2(color_breaks, custom_colors),  
      row_title = NULL,           
      column_title = NULL,        
      rect_gp = gpar(col = "black", lwd = 0.6),
      heatmap_legend_param = list(
        title = "Z-score",
        title_gp = gpar(fontsize = legend_font_size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend_font_size),
        legend_width = unit(0.8, "cm"),
        direction = "vertical",
        title_position = "topcenter"
      ),
      use_raster = FALSE
    )
    
    # Draw heatmap
    draw(ht, 
         heatmap_legend_side = "right",
         annotation_legend_side = "right",
         padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    
    dev.off()
    
    cat("Heatmap successfully saved to:", output_file, "\n")
    return(output_file)
    
  }, error = function(e) {
    cat("ERROR in create_heatmap_streamlit:", e$message, "\n")
    stop(e)
  })
}

# Run if called from command line
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 8) {
    tryCatch({
      result <- create_heatmap_streamlit(args)
      cat("SUCCESS:", result, "\n")
    }, error = function(e) {
      cat("FATAL ERROR:", e$message, "\n")
      quit(status = 1)
    })
  } else {
    cat("Usage: Rscript heatmap.R <count_file> <metadata_file> <output_dir> <plot_type> <genes> <width> <height> <dpi> [color_palette] [gene_annotation_file]\n")
    quit(status = 1)
  }
}