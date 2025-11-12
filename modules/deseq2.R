# Bulk RNA-Seq Integrated Analysis Script
# Combines DESeq2 analysis with visualization functions

# Load required libraries
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(tools)
library(paletteer)
library(jsonlite)
library(grid)

# Function to standardize ALL column names across all data frames
standardize_all_column_names <- function(df) {
  if (!is.data.frame(df)) return(df)
  
  original_names <- colnames(df)
  new_colnames <- tolower(original_names)
  
  # More comprehensive cleaning
  new_colnames <- gsub("[\\s.-]+", "_", new_colnames)  # spaces, dots, hyphens to underscore
  new_colnames <- gsub("[^a-zA-Z0-9_]", "", new_colnames)  # remove special chars
  new_colnames <- gsub("_+", "_", new_colnames)  # collapse multiple underscores
  new_colnames <- gsub("^_|_$", "", new_colnames)  # trim leading/trailing underscores
  
  # Handle empty column names
  empty_mask <- new_colnames == ""
  if (any(empty_mask)) {
    new_colnames[empty_mask] <- paste0("column_", which(empty_mask))
  }
  
  colnames(df) <- new_colnames
  
  cat("🔧 COLUMN NAME STANDARDIZATION:\n")
  cat("   Original: ", paste(original_names, collapse = ", "), "\n")
  cat("   Standardized: ", paste(new_colnames, collapse = ", "), "\n")
  
  return(df)
}

# Function to validate count matrix structure and identify required columns (HEAVILY DEBUGGED VERSION)
validate_count_matrix_structure <- function(count_matrix, analysis_type) {
  cat("🔍 ===== VALIDATING COUNT MATRIX STRUCTURE =====\n")
  cat("Analysis type:", analysis_type, "\n")
  cat("Original column names:", paste(colnames(count_matrix), collapse = ", "), "\n")
  cat("Count matrix dimensions:", nrow(count_matrix), "x", ncol(count_matrix), "\n")
  
  # Standardize column names first
  count_matrix <- standardize_all_column_names(count_matrix)
  
  cat("After standardization - Column names:", paste(colnames(count_matrix), collapse = ", "), "\n")
  cat("After standardization - First few values of first column:", head(count_matrix[[1]]), "\n")
  
  if (analysis_type == "gene") {
    # EXTENSIVE DEBUGGING FOR GENE MATRIX
    cat("🔬 GENE ANALYSIS - Checking for identifier columns...\n")
    
    # Check what columns we actually have
    available_cols <- colnames(count_matrix)
    cat("Available columns:", paste(available_cols, collapse = ", "), "\n")
    
    # Check for various possible gene ID column names
    possible_id_patterns <- c("gene_id", "geneid", "id", "gene", "ensembl", "entrez", "symbol", "name")
    found_id_cols <- character()
    
    for (pattern in possible_id_patterns) {
      matching_cols <- grep(pattern, available_cols, ignore.case = TRUE, value = TRUE)
      if (length(matching_cols) > 0) {
        found_id_cols <- c(found_id_cols, matching_cols)
        cat("Found potential ID column with pattern '", pattern, "': ", paste(matching_cols, collapse = ", "), "\n", sep = "")
      }
    }
    
    # If no pattern matches, just use the first non-sample column as ID
    if (length(found_id_cols) == 0) {
      cat("⚠️ No standard ID columns found. Using first non-sample column as ID.\n")
      # Assume sample columns are numeric, find first non-numeric column
      non_sample_cols <- character()
      for (col in available_cols) {
        if (!all(grepl("^\\d+$", count_matrix[[col]]) | is.na(count_matrix[[col]]))) {
          non_sample_cols <- c(non_sample_cols, col)
        }
      }
      
      if (length(non_sample_cols) > 0) {
        found_id_cols <- non_sample_cols[1]
        cat("Using first non-sample column as ID:", found_id_cols[1], "\n")
      } else {
        cat("❌ CRITICAL: No suitable ID columns found and no non-sample columns identified\n")
        cat("Available columns:", paste(available_cols, collapse = ", "), "\n")
        stop("Gene count matrix missing required columns. Available columns: ", paste(available_cols, collapse = ", "))
      }
    }
    
    # Use the first found ID column
    id_col <- found_id_cols[1]
    secondary_id <- if (length(found_id_cols) > 1) found_id_cols[2] else NULL
    
    cat("✅ PRIMARY ID COLUMN:", id_col, "\n")
    if (!is.null(secondary_id)) cat("✅ SECONDARY ID COLUMN:", secondary_id, "\n")
    
    # Check for biotype columns
    possible_biotype_patterns <- c("biotype", "type", "feature")
    biotype_col <- NULL
    for (pattern in possible_biotype_patterns) {
      matching_cols <- grep(pattern, available_cols, ignore.case = TRUE, value = TRUE)
      if (length(matching_cols) > 0) {
        biotype_col <- matching_cols[1]
        cat("✅ BIOTYPE COLUMN:", biotype_col, "\n")
        break
      }
    }
    
    # Identify sample columns (numeric columns or columns not matching info patterns)
    info_patterns <- c("gene", "id", "symbol", "name", "biotype", "type", "chr", "start", "end", "strand", "chromosome")
    sample_cols <- character()
    
    for (col in available_cols) {
      is_info_col <- any(sapply(info_patterns, function(pattern) grepl(pattern, col, ignore.case = TRUE)))
      # Also check if column contains mostly numbers (sample counts)
      col_values <- count_matrix[[col]]
      numeric_ratio <- sum(!is.na(as.numeric(col_values))) / length(col_values)
      
      if (!is_info_col || numeric_ratio > 0.8) {
        sample_cols <- c(sample_cols, col)
      }
    }
    
    # Remove the ID column from sample columns if it got included
    sample_cols <- setdiff(sample_cols, c(id_col, secondary_id, biotype_col))
    
    cat("✅ SAMPLE COLUMNS:", paste(sample_cols, collapse = ", "), "\n")
    
    if (length(sample_cols) == 0) {
      cat("❌ ERROR: No sample columns identified\n")
      stop("No sample columns found in gene count matrix")
    }
    
    # Check genomic info
    has_chr <- any(grepl("chr", available_cols, ignore.case = TRUE))
    has_start <- any(grepl("start", available_cols, ignore.case = TRUE))
    has_end <- any(grepl("end", available_cols, ignore.case = TRUE))
    
    return(list(
      id_col = id_col,
      secondary_id = secondary_id,
      biotype_col = biotype_col,
      has_biotype = !is.null(biotype_col),
      has_genomic_info = has_chr && has_start && has_end,
      sample_cols = sample_cols,
      count_matrix = count_matrix
    ))
    
  } else { 
    # TRANSCRIPT ANALYSIS - Similar extensive debugging
    cat("🔬 TRANSCRIPT ANALYSIS - Checking for identifier columns...\n")
    
    available_cols <- colnames(count_matrix)
    cat("Available columns:", paste(available_cols, collapse = ", "), "\n")
    
    # Check for various possible transcript ID column names
    possible_id_patterns <- c("transcript_id", "transcriptid", "tx_id", "txid", "transcript", "id")
    found_id_cols <- character()
    
    for (pattern in possible_id_patterns) {
      matching_cols <- grep(pattern, available_cols, ignore.case = TRUE, value = TRUE)
      if (length(matching_cols) > 0) {
        found_id_cols <- c(found_id_cols, matching_cols)
        cat("Found potential ID column with pattern '", pattern, "': ", paste(matching_cols, collapse = ", "), "\n", sep = "")
      }
    }
    
    if (length(found_id_cols) == 0) {
      cat("⚠️ No standard transcript ID columns found. Using first non-sample column as ID.\n")
      non_sample_cols <- character()
      for (col in available_cols) {
        if (!all(grepl("^\\d+$", count_matrix[[col]]) | is.na(count_matrix[[col]]))) {
          non_sample_cols <- c(non_sample_cols, col)
        }
      }
      
      if (length(non_sample_cols) > 0) {
        found_id_cols <- non_sample_cols[1]
        cat("Using first non-sample column as ID:", found_id_cols[1], "\n")
      } else {
        cat("❌ CRITICAL: No suitable ID columns found\n")
        stop("Transcript count matrix missing required columns. Available columns: ", paste(available_cols, collapse = ", "))
      }
    }
    
    id_col <- found_id_cols[1]
    secondary_id <- if (length(found_id_cols) > 1) found_id_cols[2] else NULL
    
    cat("✅ PRIMARY ID COLUMN:", id_col, "\n")
    if (!is.null(secondary_id)) cat("✅ SECONDARY ID COLUMN:", secondary_id, "\n")
    
    # Rest of transcript analysis logic...
    # [Keep the existing transcript analysis logic but add similar debugging]
    
    # For now, let me provide a simplified version:
    possible_biotype_patterns <- c("biotype", "type", "feature")
    biotype_col <- NULL
    for (pattern in possible_biotype_patterns) {
      matching_cols <- grep(pattern, available_cols, ignore.case = TRUE, value = TRUE)
      if (length(matching_cols) > 0) {
        biotype_col <- matching_cols[1]
        break
      }
    }
    
    info_patterns <- c("transcript", "tx", "gene", "id", "biotype", "type", "chr", "start", "end", "strand")
    sample_cols <- setdiff(available_cols, c(id_col, secondary_id, biotype_col))
    
    # Filter to likely sample columns (those that might contain counts)
    final_sample_cols <- character()
    for (col in sample_cols) {
      col_values <- count_matrix[[col]]
      numeric_ratio <- sum(!is.na(as.numeric(col_values))) / length(col_values)
      if (numeric_ratio > 0.5) {
        final_sample_cols <- c(final_sample_cols, col)
      }
    }
    
    if (length(final_sample_cols) == 0) {
      cat("❌ ERROR: No sample columns identified in transcript matrix\n")
      stop("No sample columns found in transcript count matrix")
    }
    
    cat("✅ SAMPLE COLUMNS:", paste(final_sample_cols, collapse = ", "), "\n")
    
    return(list(
      id_col = id_col,
      secondary_id = secondary_id,
      biotype_col = biotype_col,
      has_biotype = !is.null(biotype_col),
      has_genomic_info = FALSE, # Simplified for now
      sample_cols = final_sample_cols,
      count_matrix = count_matrix
    ))
  }
}

# Function to handle duplicate identifiers
process_duplicate_identifiers <- function(count_matrix, id_col, secondary_id) {
  # Handle NA or empty identifiers
  na_mask <- is.na(count_matrix[[id_col]]) | count_matrix[[id_col]] == ""
  if (any(na_mask) && !is.null(secondary_id) && secondary_id %in% colnames(count_matrix)) {
    count_matrix[[id_col]][na_mask] <- count_matrix[[secondary_id]][na_mask]
  }
  
  # Handle duplicates by appending sequential numbers
  dup_mask <- duplicated(count_matrix[[id_col]]) | duplicated(count_matrix[[id_col]], fromLast = TRUE)
  if (any(dup_mask)) {
    dup_ids <- count_matrix[[id_col]][dup_mask]
    for (id in unique(dup_ids)) {
      mask <- count_matrix[[id_col]] == id
      if (sum(mask) > 1) {
        count_matrix[[id_col]][mask] <- paste0(count_matrix[[id_col]][mask], "_", seq_len(sum(mask)))
      }
    }
  }
  
  return(count_matrix)
}

# Set up logging
setup_logging <- function(output_dir) {
  log_file <- file.path(output_dir, "analysis_log.txt")
  sink(log_file, append = TRUE, split = TRUE)
  cat("=== DESeq2 Analysis Log ===\n")
  cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=================================\n\n")
  return(log_file)
}

# Function to close logging
close_logging <- function() {
  cat("\n=================================\n")
  cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  sink()
}

# Function to create directories
create_directories <- function(base_path, contrast_name, analysis_type) {
  main_dir <- file.path(base_path, analysis_type, contrast_name)
  biotype_dir <- file.path(main_dir, paste0(analysis_type, "_biotypes"))
  
  dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(biotype_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(list(main_dir = main_dir, biotype_dir = biotype_dir))
}

# Helper function to validate and set default parameters
validate_deseq2_params <- function(params) {
  if (is.null(params)) stop("DESeq2 parameters are NULL")
  
  if (is.null(params$lfc_value)) params$lfc_value <- 1
  if (is.null(params$test)) params$test <- "Wald"
  if (is.null(params$altHypothesis)) params$altHypothesis <- "greaterAbs"
  if (is.null(params$analysis_choice)) params$analysis_choice <- "1"
  if (is.null(params$contrasts)) params$contrasts <- list()
  if (is.null(params$fitType)) params$fitType <- "parametric"  # Default to parametric
  
  if (is.null(params$metadata_file)) stop("Metadata file is required for DESeq2 analysis")
  if (!file.exists(params$metadata_file)) stop("Metadata file not found: ", params$metadata_file)
  
  return(params)
}

# Function to identify and clean sample columns
identify_and_clean_sample_columns <- function(count_matrix, sample_cols, analysis_type) {
  # Clean sample names (remove _dedup suffix if present)
  clean_sample_names <- gsub("_dedup$", "", sample_cols)
  
  # Update column names in the count matrix
  for (i in seq_along(sample_cols)) {
    if (sample_cols[i] %in% colnames(count_matrix)) {
      colnames(count_matrix)[colnames(count_matrix) == sample_cols[i]] <- clean_sample_names[i]
    }
  }
  
  return(list(count_matrix = count_matrix, sample_columns = clean_sample_names))
}

# Function to validate sample matching
validate_sample_matching <- function(sample_columns, coldata, analysis_type) {
  metadata_samples <- rownames(coldata)
  
  missing_in_metadata <- setdiff(sample_columns, metadata_samples)
  missing_in_counts <- setdiff(metadata_samples, sample_columns)
  
  if (length(missing_in_metadata) > 0) {
    cat("Warning: Samples missing in metadata:", paste(missing_in_metadata, collapse = ", "), "\n")
  }
  if (length(missing_in_counts) > 0) {
    cat("Warning: Samples missing in count matrix:", paste(missing_in_counts, collapse = ", "), "\n")
  }
  
  common_samples <- intersect(sample_columns, metadata_samples)
  if (length(common_samples) == 0) {
    stop("Error: No common samples found between count matrix and metadata.")
  }
  
  cat("Found", length(common_samples), "common samples between", analysis_type, "count matrix and metadata\n")
  return(common_samples)
}

# Function to create summary file (modified for flexible biotypes)
create_summary_file <- function(deg_results, analysis_type, output_dir, contrast_name, has_biotype) {
  cat("Creating summary file...\n")
  
  if (analysis_type == "gene") {
    summary_file <- file.path(output_dir, paste0("gene_summary_", contrast_name, ".csv"))
  } else {
    summary_file <- file.path(output_dir, paste0("transcript_summary_", contrast_name, ".csv"))
  }
  
  if (has_biotype) {
    biotype_col <- if (analysis_type == "gene") "gene_biotype" else "transcript_biotype"
    
    summary_df <- deg_results %>%
      group_by(!!sym(biotype_col)) %>%
      summarise(
        total_genes = n(),
        significant_genes = sum(padj < 0.05, na.rm = TRUE),
        up_regulated = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
        down_regulated = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE)
      ) %>%
      arrange(desc(total_genes))
  } else {
    # Summary without biotype grouping
    summary_df <- data.frame(
      total_features = nrow(deg_results),
      significant_features = sum(deg_results$padj < 0.05, na.rm = TRUE),
      up_regulated = sum(deg_results$padj < 0.05 & deg_results$log2FoldChange > 0, na.rm = TRUE),
      down_regulated = sum(deg_results$padj < 0.05 & deg_results$log2FoldChange < 0, na.rm = TRUE)
    )
  }
  
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat("Summary file saved:", summary_file, "\n")
}

# Function to select top genes for heatmap (modified for flexible identifiers)
select_top_genes <- function(deg_results, normalized_counts, analysis_type, id_col, top_n = 10) {
  sig_genes <- deg_results[!is.na(deg_results$padj) & deg_results$padj < 0.05, ]
  if (nrow(sig_genes) == 0) {
    cat("No significant genes found for heatmap\n")
    return(NULL)
  }
  
  up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
  down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]
  
  up_genes <- up_genes[order(up_genes$log2FoldChange, decreasing = TRUE), ]
  down_genes <- down_genes[order(down_genes$log2FoldChange, decreasing = FALSE), ]
  
  top_up <- head(up_genes, min(top_n, nrow(up_genes)))
  top_down <- head(down_genes, min(top_n, nrow(down_genes)))
  top_genes <- rbind(top_up, top_down)
  
  if (nrow(top_genes) < 3) {
    cat("Not enough significant genes for heatmap (need at least 3)\n")
    return(NULL)
  }
  
  # Get the feature IDs for the top genes
  feature_ids <- top_genes[[id_col]]
  
  # Subset normalized counts
  if (id_col %in% colnames(normalized_counts)) {
    count_subset <- normalized_counts[normalized_counts[[id_col]] %in% feature_ids, ]
    count_subset <- count_subset[match(feature_ids, count_subset[[id_col]]), ]
    feature_names <- count_subset[[id_col]]
  } else {
    cat("Identifier column not found in normalized counts\n")
    return(NULL)
  }
  
  return(list(top_genes = top_genes, count_subset = count_subset, 
              feature_names = feature_names, feature_order = feature_ids))
}

# Function to create enhanced heatmap (modified for flexible structure)
create_enhanced_heatmap <- function(normalized_counts, deg_results, output_path, 
                                    analysis_type, coldata, contrast_name, id_col, top_n = 10) {
  tryCatch({
    conditions <- strsplit(contrast_name, "_vs_")[[1]]
    if (length(conditions) != 2) {
      cat("Warning: Contrast name format should be 'Condition1_vs_Condition2'. Using default sample order.\n")
      condition_order <- NULL
    } else {
      condition1 <- conditions[1]
      condition2 <- conditions[2]
      condition_order <- c(condition1, condition2)
    }
    
    top_data <- select_top_genes(deg_results, normalized_counts, analysis_type, id_col, top_n)
    if (is.null(top_data)) return(NULL)
    
    top_genes <- top_data$top_genes
    count_subset <- top_data$count_subset
    
    # Identify sample columns
    sample_cols <- setdiff(colnames(count_subset), colnames(deg_results))
    
    count_matrix <- as.matrix(count_subset[, sample_cols])
    rownames(count_matrix) <- make.names(count_subset[[id_col]], unique = TRUE)
    
    if (nrow(count_matrix) < 3) {
      cat("Not enough genes with normalized counts for heatmap, skipping...\n")
      return(NULL)
    }
    
    common_samples <- intersect(colnames(count_matrix), rownames(coldata))
    if (length(common_samples) == 0) {
      cat("No common samples between count matrix and coldata, skipping heatmap\n")
      return(NULL)
    }
    
    count_matrix <- count_matrix[, common_samples, drop = FALSE]
    annotation_data <- coldata[common_samples, , drop = FALSE]
    
    na_condition <- is.na(annotation_data$condition)
    if (any(na_condition)) {
      cat("Warning: Found", sum(na_condition), "samples with NA condition. Removing them.\n")
      annotation_data <- annotation_data[!na_condition, , drop = FALSE]
      count_matrix <- count_matrix[, rownames(annotation_data), drop = FALSE]
    }
    
    if (!is.null(condition_order) && all(condition_order %in% annotation_data$condition)) {
      condition1_samples <- rownames(annotation_data)[annotation_data$condition == condition1]
      condition2_samples <- rownames(annotation_data)[annotation_data$condition == condition2]
      
      if (length(condition1_samples) > 0 && length(condition2_samples) > 0) {
        sample_order <- c(condition1_samples, condition2_samples)
        count_matrix <- count_matrix[, sample_order, drop = FALSE]
        annotation_data <- annotation_data[sample_order, , drop = FALSE]
      }
    }
    
    if (ncol(count_matrix) < 2) {
      cat("Not enough samples after removing NA conditions for heatmap, skipping...\n")
      return(NULL)
    }
    
    valid_conditions <- unique(annotation_data$condition)
    condition_colors <- brewer.pal(min(8, length(valid_conditions)), "Set1")[1:length(valid_conditions)]
    names(condition_colors) <- valid_conditions
    annotation_colors <- list(condition = condition_colors)
    
    output_dir <- dirname(output_path)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    count_matrix_scaled <- t(scale(t(count_matrix)))
    count_matrix_scaled[is.na(count_matrix_scaled)] <- 0
    count_matrix_scaled[is.infinite(count_matrix_scaled)] <- 0
    
    pheatmap::pheatmap(count_matrix_scaled,
                       main = paste("Top", nrow(top_genes), "DEGs by log2FC -", contrast_name),
                       scale = "none",
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation",
                       clustering_method = "complete",
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       fontsize_row = 8,
                       fontsize_col = 8,
                       annotation_col = annotation_data,
                       rect_gp = gpar(col = "black", lwd = 0.6),
                       annotation_colors = annotation_colors,
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       filename = output_path,
                       width = 10,
                       height = 8)
    
    cat("Heatmap successfully saved:", output_path, "\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("Enhanced heatmap failed:", e$message, "\n")
    return(FALSE)
  })
}

# Function to generate enhanced volcano plot
generate_enhanced_volcano <- function(deg_results, output_path, analysis_type, contrast_name, lfc_threshold = 1) {
  tryCatch({
    p_volcano <- EnhancedVolcano(deg_results,
                                 lab = NA,
                                 x = 'log2FoldChange',
                                 y = 'pvalue',
                                 title = contrast_name,
                                 pCutoff = 0.05,
                                 FCcutoff = lfc_threshold,
                                 pointSize = 2.0,
                                 labSize = 3.0,
                                 colAlpha = 0.7,
                                 legendPosition = 'right',
                                 col = c('grey50', '#1DCD9F', 'royalblue', '#FF0B55')) + 
      theme(panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey90", linewidth = 0.25))
    
    ggsave(output_path, p_volcano, width = 12, height = 8, dpi = 300)
    
  }, error = function(e) {
    cat("Enhanced volcano plot failed:", e$message, "\n")
  })
}

# Function to perform DESeq2 analysis (MAJOR REVISION)
run_deseq2_analysis <- function(count_matrix, coldata, contrast, 
                                analysis_type, output_dirs, params, contrast_name) {
  
  cat("Running DESeq2 analysis for", analysis_type, "...\n")
  
  # Validate count matrix structure and get column information
  matrix_info <- validate_count_matrix_structure(count_matrix, analysis_type)
  count_matrix <- matrix_info$count_matrix
  id_col <- matrix_info$id_col
  has_biotype <- matrix_info$has_biotype
  biotype_col <- matrix_info$biotype_col
  sample_cols <- matrix_info$sample_cols
  
  # Process duplicate identifiers
  count_matrix <- process_duplicate_identifiers(count_matrix, id_col, matrix_info$secondary_id)
  
  # Clean sample columns
  cleaned_data <- identify_and_clean_sample_columns(count_matrix, sample_cols, analysis_type)
  count_matrix <- cleaned_data$count_matrix
  sample_columns <- cleaned_data$sample_columns
  
  # Validate sample matching
  common_samples <- validate_sample_matching(sample_columns, coldata, analysis_type)
  
  # Prepare count data for DESeq2
  count_data <- count_matrix[, common_samples, drop = FALSE]
  rownames(count_data) <- count_matrix[[id_col]]
  
  coldata_common <- coldata[common_samples, , drop = FALSE]
  coldata_common$condition <- as.factor(coldata_common$condition)
  
  # Handle different test types with configurable fitType
  fit_type <- params$fitType %||% "parametric"
  cat("Using fitType:", fit_type, "\n")
  
  if (params$test == "LRT") {
    # LRT test - handle multiple model comparisons
    if (is.null(params$lrt_comparisons) || length(params$lrt_comparisons) == 0) {
      stop("No LRT comparisons defined")
    }
    
    cat("Running LRT analysis with", length(params$lrt_comparisons), "model comparisons\n")
    
    # Store results for all comparisons
    all_results <- list()
    
    for (i in seq_along(params$lrt_comparisons)) {
      comparison <- params$lrt_comparisons[[i]]
      comparison_name <- comparison$name %||% paste0("LRT_comparison_", i)
      
      cat("Running LRT comparison:", comparison_name, "\n")
      cat("  Full model:", comparison$full, "\n")
      cat("  Reduced model:", comparison$reduced, "\n")
      
      # Create DESeq dataset with full design
      full_formula <- as.formula(comparison$full)
      reduced_formula <- as.formula(comparison$reduced)
      
      dds_lrt <- DESeqDataSetFromMatrix(countData = count_data,
                                        colData = coldata_common,
                                        design = full_formula)
      
      # Run LRT test with configurable fitType
      dds_lrt <- DESeq(dds_lrt, test = "LRT", reduced = reduced_formula, fitType = fit_type)
      res_lrt <- results(dds_lrt, alpha = 0.05)
      
      # Store results
      all_results[[comparison_name]] <- list(
        results = res_lrt,
        full_model = comparison$full,
        reduced_model = comparison$reduced
      )
      
      # For the first comparison, use it for downstream analysis
      if (i == 1) {
        dds <- dds_lrt
        res <- res_lrt
      }
    }
    
    # Save all LRT results
    saveRDS(all_results, file.path(output_dirs$main_dir, "all_lrt_results.rds"))
    
  } else {
    # Wald test - original code with configurable fitType
    dds <- DESeqDataSetFromMatrix(countData = count_data,
                                  colData = coldata_common,
                                  design = ~ condition)
    
    dds <- DESeq(dds, test = params$test, fitType = fit_type)
    res <- results(dds, contrast = contrast, alpha = 0.05, altHypothesis = params$altHypothesis)
  }
  
  saveRDS(dds, file.path(output_dirs$main_dir, "dds_object.rds"))
  
  # Create results data frame with flexible structure
  deg_results <- data.frame(
    feature_id = count_matrix[[id_col]],
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj,
    row.names = NULL
  )
  
  # Add biotype information if available
  if (has_biotype && !is.null(biotype_col)) {
    deg_results[[biotype_col]] <- count_matrix[[biotype_col]]
  }
  
  # Add secondary identifier if available
  if (!is.null(matrix_info$secondary_id) && matrix_info$secondary_id %in% colnames(count_matrix)) {
    sec_id_name <- if (analysis_type == "gene") "gene_symbol" else "transcript_name"
    deg_results[[sec_id_name]] <- count_matrix[[matrix_info$secondary_id]]
  }
  
  deg_results <- deg_results[complete.cases(deg_results), ]
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Create normalized count matrix with flexible structure
  norm_count_df <- data.frame(
    feature_id = count_matrix[[id_col]],
    normalized_counts
  )
  
  # Add secondary identifier if available
  if (!is.null(matrix_info$secondary_id) && matrix_info$secondary_id %in% colnames(count_matrix)) {
    norm_count_df[[matrix_info$secondary_id]] <- count_matrix[[matrix_info$secondary_id]]
  }
  
  # Add biotype if available
  if (has_biotype && !is.null(biotype_col)) {
    norm_count_df[[biotype_col]] <- count_matrix[[biotype_col]]
  }
  
  # Save normalized counts
  if (analysis_type == "gene") {
    write.csv(norm_count_df, file.path(params$output_dir, "normalized_gene_count_matrix.csv"), row.names = FALSE)
  } else {
    write.csv(norm_count_df, file.path(params$output_dir, "normalized_transcript_count_matrix.csv"), row.names = FALSE)
  }
  
  generate_plots(dds, deg_results, norm_count_df, output_dirs, analysis_type, params, contrast_name, coldata_common, id_col, has_biotype, biotype_col)
  
  if (has_biotype && !is.null(biotype_col)) {
    process_biotypes(deg_results, norm_count_df, output_dirs, analysis_type, params, contrast_name, coldata_common, id_col, biotype_col)
  }
  
  create_summary_file(deg_results, analysis_type, output_dirs$main_dir, contrast_name, has_biotype)
  
  write.csv(deg_results, file.path(output_dirs$main_dir, "overall_deg_results.csv"), row.names = FALSE)
  cat("Analysis completed for", analysis_type, "in", output_dirs$main_dir, "\n")
}

# Generate plots (modified for flexible structure)
generate_plots <- function(dds, deg_results, normalized_counts, output_dirs, analysis_type, params, contrast_name, coldata, id_col, has_biotype, biotype_col) {
  cat("Generating plots...\n")
  
  plot_colors <- brewer.pal(8, "Set1")
  
  tryCatch({
    vsd <- vst(dds, blind = TRUE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    
    p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
      geom_point(size = 3) +
      scale_color_manual(values = plot_colors) +
      theme_minimal() +
      ggtitle("PCA - Raw Counts") +
      theme(panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey90", linewidth = 0.25))
    
    ggsave(file.path(output_dirs$main_dir, "pca_raw.png"), p, width = 8, height = 6, dpi = 300)
  }, error = function(e) cat("PCA for raw counts failed:", e$message, "\n"))
  
  tryCatch({
    vsd_norm <- vst(dds, blind = FALSE)
    pca_data_norm <- plotPCA(vsd_norm, intgroup = "condition", returnData = TRUE)
    
    p_norm <- ggplot(pca_data_norm, aes(PC1, PC2, color = condition)) +
      geom_point(size = 3) +
      scale_color_manual(values = plot_colors) +
      theme_minimal() +
      ggtitle("PCA - Normalized Counts") +
      theme(panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey90", linewidth = 0.25))
    
    ggsave(file.path(output_dirs$main_dir, "pca_normalized.png"), p_norm, width = 8, height = 6, dpi = 300)
  }, error = function(e) cat("PCA for normalized counts failed:", e$message, "\n"))
  
  generate_enhanced_volcano(deg_results, file.path(output_dirs$main_dir, "volcano_plot.png"),
                            analysis_type, contrast_name, params$lfc_value)
  
  create_enhanced_heatmap(normalized_counts, deg_results, file.path(output_dirs$main_dir, "enhanced_heatmap.png"),
                          analysis_type, coldata, contrast_name, id_col, top_n = 10)
}

# Process by biotypes (modified to handle optional biotypes)
process_biotypes <- function(deg_results, normalized_counts, output_dirs, analysis_type, params, contrast_name, coldata, id_col, biotype_col) {
  cat("Processing by biotypes...\n")
  
  biotypes <- unique(deg_results[[biotype_col]])
  
  for (biotype in biotypes) {
    biotype_subdir <- file.path(output_dirs$biotype_dir, biotype)
    dir.create(biotype_subdir, recursive = TRUE, showWarnings = FALSE)
    
    biotype_results <- deg_results[deg_results[[biotype_col]] == biotype, ]
    write.csv(biotype_results, file.path(biotype_subdir, paste0("deg_results_", biotype, ".csv")), row.names = FALSE)
    
    if (nrow(biotype_results) > 0) {
      generate_enhanced_volcano(biotype_results, file.path(biotype_subdir, paste0("volcano_plot_", biotype, ".png")),
                                analysis_type, paste0(contrast_name, " - ", biotype), params$lfc_value)
      
      if (nrow(biotype_results) >= 3) {
        create_enhanced_heatmap(normalized_counts, biotype_results, file.path(biotype_subdir, paste0("heatmap_", biotype, ".png")),
                                analysis_type, coldata, paste0(contrast_name, " - ", biotype), id_col, top_n = 10)
      } else {
        cat("Skipping heatmap for", biotype, "- less than 3 genes\n")
      }
    }
  }
}

# Function to clean condition names
clean_condition_names <- function(conditions) {
  if (is.atomic(conditions) && length(conditions) == 1) conditions <- as.character(conditions)
  gsub("[^a-zA-Z0-9_]", "_", conditions)
}

# Function to run DESeq2 with parameters from Python (MAJOR REVISION)
run_deseq2_from_params <- function(params) {
  cat("DESeq2 Analysis with Python Parameters\n")
  cat("======================================\n")
  
  deseq2_params <- validate_deseq2_params(params$deseq2_params)
  output_dir <- params$output_dir
  
  log_file <- setup_logging(output_dir)
  
  cat("Reading metadata file...\n")
  coldata <- read.csv(deseq2_params$metadata_file, row.names = 1)
  
  if (!is.data.frame(coldata)) stop("Error: Metadata file must be a valid data frame")
  if (!"condition" %in% colnames(coldata)) {
    cat("Error: 'condition' column not found in metadata\n")
    close_logging()
    return()
  }
  
  na_conditions <- is.na(coldata$condition)
  if (any(na_conditions)) {
    cat("Warning: Found", sum(na_conditions), "samples with NA condition. Removing them.\n")
    coldata <- coldata[!na_conditions, , drop = FALSE]
  }
  
  coldata$condition <- clean_condition_names(coldata$condition)
  coldata$condition <- as.factor(coldata$condition)
  
  lfc_value <- if (!is.null(deseq2_params$lfc_value)) as.numeric(deseq2_params$lfc_value) else 1
  test <- if (!is.null(deseq2_params$test)) as.character(deseq2_params$test) else "Wald"
  altHypothesis <- if (!is.null(deseq2_params$altHypothesis)) as.character(deseq2_params$altHypothesis) else "greaterAbs"
  fit_type <- if (!is.null(deseq2_params$fitType)) as.character(deseq2_params$fitType) else "parametric"
  
  params_list <- list(lfc_value = lfc_value, test = test, fitType = fit_type, 
                      altHypothesis = altHypothesis, output_dir = output_dir)
  
  contrasts <- list()
  
  if (!is.null(deseq2_params$contrasts)) {
    if (is.data.frame(deseq2_params$contrasts)) {
      cat("Processing contrasts as data frame...\n")
      for (i in 1:nrow(deseq2_params$contrasts)) {
        contrast_row <- deseq2_params$contrasts[i, ]
        group1_clean <- clean_condition_names(contrast_row$group1)
        group2_clean <- clean_condition_names(contrast_row$group2)
        
        contrast_name <- if (!is.null(contrast_row$name) && contrast_row$name != "") contrast_row$name else paste(group1_clean, "vs", group2_clean, sep = "_")
        contrasts[[contrast_name]] <- c("condition", group1_clean, group2_clean)
        cat("Comparison added:", contrast_name, "\n")
      }
    } else if (is.list(deseq2_params$contrasts)) {
      cat("Processing contrasts as list...\n")
      for (i in seq_along(deseq2_params$contrasts)) {
        contrast <- deseq2_params$contrasts[[i]]
        
        if (is.list(contrast) && !is.null(contrast$group1) && !is.null(contrast$group2)) {
          group1_clean <- clean_condition_names(contrast$group1)
          group2_clean <- clean_condition_names(contrast$group2)
          
          contrast_name <- if (!is.null(contrast$name)) contrast$name else paste(group1_clean, "vs", group2_clean, sep = "_")
          contrasts[[contrast_name]] <- c("condition", group1_clean, group2_clean)
          cat("Comparison added:", contrast_name, "\n")
        } else if (is.character(contrast) && length(contrast) >= 2) {
          group1_clean <- clean_condition_names(contrast[1])
          group2_clean <- clean_condition_names(contrast[2])
          contrast_name <- if (length(contrast) >= 3) contrast[3] else paste(group1_clean, "vs", group2_clean, sep = "_")
          contrasts[[contrast_name]] <- c("condition", group1_clean, group2_clean)
          cat("Comparison added:", contrast_name, "\n")
        } else cat("Warning: Skipping invalid contrast format at index", i, "\n")
      }
    } else if (is.character(deseq2_params$contrasts) && length(deseq2_params$contrasts) >= 2) {
      group1_clean <- clean_condition_names(deseq2_params$contrasts[1])
      group2_clean <- clean_condition_names(deseq2_params$contrasts[2])
      contrast_name <- if (length(deseq2_params$contrasts) >= 3) deseq2_params$contrasts[3] else paste(group1_clean, "vs", group2_clean, sep = "_")
      contrasts[[contrast_name]] <- c("condition", group1_clean, group2_clean)
      cat("Comparison added:", contrast_name, "\n")
    }
  } else cat("Warning: No contrasts defined in parameters\n")
  
  analysis_choice <- deseq2_params$analysis_choice %||% "1"
  gene_matrix <- NULL
  transcript_matrix <- NULL
  
  if (analysis_choice %in% c("1", "3")) {
    if (!is.null(deseq2_params$gene_matrix_file)) {
      cat("🔍 Checking gene count matrix path:\n")
      cat("   Provided path:", deseq2_params$gene_matrix_file, "\n")
      cat("   File exists:", file.exists(deseq2_params$gene_matrix_file), "\n")
      cat("   Current working directory:", getwd(), "\n")
      
      if (file.exists(deseq2_params$gene_matrix_file)) {
        cat("📊 Reading gene count matrix from:", deseq2_params$gene_matrix_file, "\n")
        gene_matrix <- read.csv(deseq2_params$gene_matrix_file)
        
        cat("✅ Gene count matrix loaded successfully -", nrow(gene_matrix), "genes\n")
        
      } else {
        cat("❌ ERROR: Gene count matrix file not found\n")
        if (analysis_choice == "1") {
          close_logging()
          stop("Gene count matrix is required but not found.")
        } else {
          cat("⚠️  Skipping gene analysis\n")
        }
      }
      
    } else {
      cat("❌ ERROR: Gene matrix file path not provided\n")
    }
  }
  
  if (analysis_choice %in% c("2", "3")) {
    if (!is.null(deseq2_params$transcript_matrix_file)) {
      cat("🔍 Checking transcript count matrix path:\n")
      cat("   Provided path:", deseq2_params$transcript_matrix_file, "\n")
      cat("   File exists:", file.exists(deseq2_params$transcript_matrix_file), "\n")
      
      if (file.exists(deseq2_params$transcript_matrix_file)) {
        cat("📊 Reading transcript count matrix from:", deseq2_params$transcript_matrix_file, "\n")
        transcript_matrix <- read.csv(deseq2_params$transcript_matrix_file)
        
        cat("✅ Transcript count matrix loaded successfully -", nrow(transcript_matrix), "transcripts\n")
        
      } else {
        cat("❌ ERROR: Transcript count matrix file not found\n")
        if (analysis_choice == "2") {
          close_logging()
          stop("Transcript count matrix is required but not found.")
        } else {
          cat("⚠️  Skipping transcript analysis\n")
        }
      }
      
    } else {
      cat("❌ ERROR: Transcript matrix file path not provided\n")
    }
  }
  
  if (is.null(gene_matrix) && is.null(transcript_matrix)) {
    cat("❌ CRITICAL ERROR: No count matrices found for analysis\n")
    close_logging()
    stop("No count matrices found. DESeq2 analysis cannot proceed.")
  }
  
  if (length(contrasts) == 0) {
    cat("Error: No valid contrasts were created.\n")
    close_logging()
    return()
  } else cat("Successfully created", length(contrasts), "contrast(s)\n")
  
  analysis_ran <- FALSE
  
  if (analysis_choice %in% c("1", "3") && !is.null(gene_matrix)) {
    analysis_ran <- TRUE
    for (contrast_name in names(contrasts)) {
      cat("🔬 Running gene analysis for comparison:", contrast_name, "\n")
      output_dirs <- create_directories(output_dir, contrast_name, "gene")
      run_deseq2_analysis(gene_matrix, coldata, contrasts[[contrast_name]], "gene", output_dirs, params_list, contrast_name)
    }
  }
  
  if (analysis_choice %in% c("2", "3") && !is.null(transcript_matrix)) {
    analysis_ran <- TRUE
    for (contrast_name in names(contrasts)) {
      cat("🔬 Running transcript analysis for comparison:", contrast_name, "\n")
      output_dirs <- create_directories(output_dir, contrast_name, "transcript")
      run_deseq2_analysis(transcript_matrix, coldata, contrasts[[contrast_name]], "transcript", output_dirs, params_list, contrast_name)
    }
  }
  
  if (analysis_ran) {
    cat("✅ DESeq2 analysis completed successfully!\n📁 Results saved in:", output_dir, "\n")
  } else {
    cat("❌ DESeq2 analysis failed - no count matrices were processed\n")
    close_logging()
    stop("DESeq2 analysis failed - no count matrices processed")
  }
  
  close_logging()
}

# Helper function for default values
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Run the main function if script is executed directly
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (file.exists(args[1]) && tools::file_ext(args[1]) == "json") {
    run_deseq2_from_params(jsonlite::fromJSON(args[1]))
  } else if (args[1] == "deseq2") {
    # For standalone DESeq2 execution
    source("deseq2_functions.R")  # Assuming helper functions are in separate file
    run_deseq2_main()
  } else cat("Invalid argument. Use: 'deseq2' or provide a JSON parameter file.\n")
} else {
  cat("DESeq2 analysis script loaded successfully.\n")
  cat("Use: Rscript deseq2.R <parameter_file.json> to run analysis.\n")
}