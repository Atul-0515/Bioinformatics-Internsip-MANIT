#!/usr/bin/env Rscript

# Clean R Package Installer - Only installs missing packages
cat("Checking and installing missing R packages...\n")

# Helper function to install only missing packages
install_missing <- function(packages, source = "CRAN") {
    missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

    if (length(missing) == 0) {
        cat("All", source, "packages already installed.\n")
        return()
    }

    cat("Installing", length(missing), source, "packages:", paste(missing, collapse = ", "), "\n")

    if (source == "CRAN") {
        install.packages(missing, dependencies = TRUE, repos = "https://cloud.r-project.org/")
    } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "https://cloud.r-project.org/")
        }
        BiocManager::install(missing, update = FALSE, ask = FALSE)
    }
}

# CRAN packages
cran_packages <- c(
    "plotly", "dplyr", "tidyr", "readr", "ggplot2", "gridExtra", "RColorBrewer",
    "vegan", "jsonlite", "reshape2", "stringr", "lubridate", "corrplot",
    "pheatmap", "cluster", "factoextra", "cowplot", "ggrepel", "viridis",
    "scales", "readxl", "writexl", "openxlsx", "parallel", "doParallel",
    "igraph", "randomForest", "caret"
)

# Bioconductor packages
bioc_packages <- c(
    "VennDiagram", "ComplexHeatmap", "circlize", "ALDEx2", "phyloseq",
    "Biostrings", "GenomeInfoDb", "microbiome", "DESeq2", "edgeR",
    "ggtree", "treeio", "metagenomeSeq"
)

# Install missing packages
tryCatch(
    {
        install_missing(cran_packages, "CRAN")
        install_missing(bioc_packages, "Bioconductor")
        cat("✓ Package installation complete!\n")
    },
    error = function(e) {
        cat("✗ Error:", e$message, "\n")
    }
)
