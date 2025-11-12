# Common functions for visualization scripts

# Function to read metadata with encoding handling
read_metadata_common <- function(file_path) {
  library(tools)
  
  if (tools::file_ext(file_path) == "csv") {
    metadata <- tryCatch({
      read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8", 
               check.names = FALSE, strip.white = TRUE)
    }, error = function(e) {
      tryCatch({
        read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", 
                 check.names = FALSE, strip.white = TRUE)
      }, error = function(e) {
        tryCatch({
          read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1", 
                   check.names = FALSE, strip.white = TRUE)
        }, error = function(e) {
          read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1", 
                   check.names = FALSE, strip.white = TRUE, fill = TRUE)
        })
      })
    })
  } else {
    metadata <- tryCatch({
      read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8", 
                 check.names = FALSE, strip.white = TRUE)
    }, error = function(e) {
      tryCatch({
        read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", 
                   check.names = FALSE, strip.white = TRUE)
      }, error = function(e) {
        tryCatch({
          read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1", 
                     check.names = FALSE, strip.white = TRUE)
        }, error = function(e) {
          read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "latin1", 
                     check.names = FALSE, strip.white = TRUE, fill = TRUE)
        })
      })
    })
  }
  
  # Clean column names
  colnames(metadata) <- gsub("^\ufeff", "", colnames(metadata))
  colnames(metadata) <- iconv(colnames(metadata), to = "ASCII//TRANSLIT")
  colnames(metadata) <- gsub("[^a-zA-Z0-9._]", "_", colnames(metadata))
  
  return(metadata)
}

# Function to detect sample and condition columns
detect_metadata_columns <- function(metadata) {
  # Look for common sample column names
  sample_keywords <- c("sample", "id", "name", "sample_id", "sample_name", "samples")
  condition_keywords <- c("condition", "group", "treatment", "type", "class", "status")
  
  sample_col <- NULL
  condition_col <- NULL
  
  # Find sample column
  for (col in colnames(metadata)) {
    if (tolower(col) %in% sample_keywords) {
      sample_col <- col
      break
    }
  }
  
  # If no sample column found, use first column
  if (is.null(sample_col) && ncol(metadata) > 0) {
    sample_col <- colnames(metadata)[1]
  }
  
  # Find condition column
  for (col in colnames(metadata)) {
    if (tolower(col) %in% condition_keywords && col != sample_col) {
      condition_col <- col
      break
    }
  }
  
  # If no condition column found, use second column if available
  if (is.null(condition_col) && ncol(metadata) > 1) {
    condition_col <- colnames(metadata)[2]
  }
  
  return(list(sample_col = sample_col, condition_col = condition_col))
}

# Function to read gene annotation file
read_gene_annotation_common <- function(file_path, feature_names) {
  library(tools)
  
  cat("Reading gene annotation file:", file_path, "\n")
  
  if (tools::file_ext(file_path) == "csv") {
    gene_annot <- tryCatch({
      read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8", 
               check.names = FALSE, strip.white = TRUE, fill = TRUE)
    }, error = function(e) {
      tryCatch({
        read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", 
                 check.names = FALSE, strip.white = TRUE, fill = TRUE)
      }, error = function(e) {
        read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1", 
                 check.names = FALSE, strip.white = TRUE, fill = TRUE, na.strings = c("", "NA", "N/A"))
      })
    })
  } else {
    gene_annot <- tryCatch({
      read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8", 
                 check.names = FALSE, strip.white = TRUE, fill = TRUE)
    }, error = function(e) {
      tryCatch({
        read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", 
                   check.names = FALSE, strip.white = TRUE, fill = TRUE)
      }, error = function(e) {
        read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "latin1", 
                   check.names = FALSE, strip.white = TRUE, fill = TRUE, na.strings = c("", "NA", "N/A"))
      })
    })
  }
  
  # Clean column names
  colnames(gene_annot) <- gsub("^\ufeff", "", colnames(gene_annot))
  colnames(gene_annot) <- iconv(colnames(gene_annot), to = "ASCII//TRANSLIT")
  colnames(gene_annot) <- gsub("[^a-zA-Z0-9._]", "_", colnames(gene_annot))
  
  if (ncol(gene_annot) < 2) {
    cat("Warning: Gene annotation file must have at least 2 columns.\n")
    return(NULL)
  }
  
  # Clean gene names in first column
  gene_annot[, 1] <- make.names(iconv(gene_annot[, 1], to = "ASCII//TRANSLIT"))
  gene_annot[, 1] <- gsub("[^a-zA-Z0-9._]", "_", gene_annot[, 1])
  
  rownames(gene_annot) <- gene_annot[, 1]
  gene_annot <- gene_annot[, -1, drop = FALSE]
  
  # Clean all character columns
  for (col in colnames(gene_annot)) {
    if (is.character(gene_annot[[col]])) {
      gene_annot[[col]] <- iconv(gene_annot[[col]], to = "ASCII//TRANSLIT")
      gene_annot[[col]] <- gsub("[^a-zA-Z0-9._ ]", "_", gene_annot[[col]])
      gene_annot[[col]][gene_annot[[col]] == ""] <- NA
    }
  }
  
  gene_annot <- gene_annot[rowSums(is.na(gene_annot)) != ncol(gene_annot), , drop = FALSE]
  
  # Subset to only include genes that are in feature_names
  available_genes <- rownames(gene_annot)[rownames(gene_annot) %in% feature_names]
  
  if (length(available_genes) == 0) {
    cat("Warning: No matching genes found in gene annotation file.\n")
    return(NULL)
  }
  
  cat("Found annotations for", length(available_genes), "genes\n")
  return(gene_annot[available_genes, , drop = FALSE])
}