# Install all packages at once
# CRAN packages
cran_packages <- c(
  "ggplot2", "dplyr", "tidyr", "RColorBrewer", "reshape2",
  "pheatmap", "circlize", "tools", "paletteer", "jsonlite", "grid"
)

# Bioconductor packages
bioc_packages <- c("DESeq2", "EnhancedVolcano", "ComplexHeatmap")

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

# Install BiocManager if needed and then Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(bioc_packages)

# Verify installations
print("CRAN packages installed:")
print(cran_packages)

print("Bioconductor packages installed:")
print(bioc_packages)