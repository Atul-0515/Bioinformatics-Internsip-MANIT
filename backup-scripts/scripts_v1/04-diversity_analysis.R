#!/usr/bin/env Rscript

# Diversity Analysis Script for 16S rRNA Taxonomic Data
# Performs alpha and beta diversity analyses on Kraken2/Bracken output
# Author: Bioinformatics Pipeline
# Date: 2025

# Load required libraries
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(pheatmap)
})

# Function to install missing packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, repos = "https://cran.r-project.org/")
  }
}

# Check and install required packages
required_packages <- c(
  "vegan", "ggplot2", "dplyr", "tidyr", "RColorBrewer", "pheatmap"
)
install_if_missing(required_packages)

# Create output directories
create_directories <- function() {
  dirs <- c(
    "results/diversity",
    "results/diversity/plots",
    "results/diversity/tables"
  )

  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

# Load abundance data and metadata
load_data <- function() {
  cat("Loading taxonomic abundance data...\n")

  # Check if taxonomy results exist
  taxonomy_dir <- "results/taxonomy/merged_tables"
  if (!dir.exists(taxonomy_dir)) {
    stop("Taxonomy results not found. Please run 03-taxonomic_profiling.py first.")
  }

  # Load sample categories
  sample_categories_file <- "results/taxonomy/sample_categories.csv"
  if (!file.exists(sample_categories_file)) {
    stop("Sample categories file not found. Please run 03-taxonomic_profiling.py first.")
  }

  metadata <- read.csv(sample_categories_file, stringsAsFactors = FALSE)
  cat("Loaded metadata for", nrow(metadata), "samples\n")

  # Load abundance tables for different taxonomic levels
  levels <- c("species", "genus", "family", "order", "class", "phylum")
  abundance_data <- list()

  for (level in levels) {
    abundance_file <- file.path(taxonomy_dir, paste0(level, "_abundance_counts.csv"))
    relative_file <- file.path(taxonomy_dir, paste0(level, "_relative_abundance.csv"))

    if (file.exists(abundance_file) && file.exists(relative_file)) {
      abundance_data[[level]] <- list(
        counts = read.csv(abundance_file, row.names = 1, check.names = FALSE),
        relative = read.csv(relative_file, row.names = 1, check.names = FALSE)
      )
      cat(
        "Loaded", level, "level data:", nrow(abundance_data[[level]]$counts), "taxa,",
        ncol(abundance_data[[level]]$counts), "samples\n"
      )
    }
  }

  if (length(abundance_data) == 0) {
    stop("No abundance data found. Please check taxonomy results.")
  }

  return(list(abundance_data = abundance_data, metadata = metadata))
}

# Calculate alpha diversity metrics
calculate_alpha_diversity <- function(abundance_data, metadata) {
  cat("Calculating alpha diversity metrics...\n")

  alpha_results <- list()

  for (level in names(abundance_data)) {
    cat("Processing", level, "level...\n")

    # Get count data (transposed for vegan functions)
    count_matrix <- t(abundance_data[[level]]$counts)

    # Calculate diversity indices
    alpha_metrics <- data.frame(
      SampleID = rownames(count_matrix),
      Observed_Taxa = apply(count_matrix > 0, 1, sum),
      Shannon = diversity(count_matrix, index = "shannon"),
      Simpson = diversity(count_matrix, index = "simpson"),
      InvSimpson = diversity(count_matrix, index = "invsimpson"),
      Chao1 = apply(count_matrix, 1, function(x) {
        if (sum(x) == 0) {
          return(0)
        }
        estimateR(x)[2] # Chao1 estimate
      }),
      stringsAsFactors = FALSE
    )

    # Add category information
    alpha_metrics <- merge(alpha_metrics, metadata, by = "SampleID", all.x = TRUE)

    alpha_results[[level]] <- alpha_metrics

    # Save results
    write.csv(alpha_metrics,
      file.path("results/diversity/tables", paste0(level, "_alpha_diversity.csv")),
      row.names = FALSE
    )
  }

  cat("Alpha diversity calculation complete.\n")
  return(alpha_results)
}

# Plot alpha diversity
plot_alpha_diversity <- function(alpha_results) {
  cat("Creating alpha diversity plots...\n")

  # Create plots for each taxonomic level
  for (level in names(alpha_results)) {
    alpha_data <- alpha_results[[level]]

    if (!"Category" %in% colnames(alpha_data)) {
      cat("Warning: No category information found for", level, "level\n")
      next
    }

    # Prepare data for plotting
    plot_data <- alpha_data %>%
      select("SampleID", "Category", "Observed_Taxa", "Shannon", "Simpson", "InvSimpson") %>%
      pivot_longer(
        cols = c("Observed_Taxa", "Shannon", "Simpson", "InvSimpson"),
        names_to = "Metric", values_to = "Value"
      )


    # Create boxplots
    p1 <- ggplot(plot_data, aes(x = .data$Category, y = .data$Value, fill = .data$Category)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      facet_wrap(~ .data$Metric, scales = "free_y") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      ) +
      labs(
        title = paste("Alpha Diversity -", stringr::str_to_title(level), "Level"),
        x = "Sample Category", y = "Diversity Value"
      ) +
      scale_fill_brewer(type = "qual", palette = "Set2")

    # Create violin plots
    p2 <- ggplot(plot_data, aes(x = .data$Category, y = .data$Value, fill = .data$Category)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.8) +
      facet_wrap(~ .data$Metric, scales = "free_y") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      ) +
      labs(
        title = paste("Alpha Diversity Distribution -", stringr::str_to_title(level), "Level"),
        x = "Sample Category", y = "Diversity Value"
      ) +
      scale_fill_brewer(type = "qual", palette = "Set1")

    ggsave(file.path("results/diversity/plots", paste0(level, "_alpha_diversity_boxplot.png")),
      p1,
      width = 12, height = 8, dpi = 300
    )

    ggsave(file.path("results/diversity/plots", paste0(level, "_alpha_diversity_violin.png")),
      p2,
      width = 12, height = 8, dpi = 300
    )
  }

  cat("Alpha diversity plots saved.\n")
}

# Calculate beta diversity
calculate_beta_diversity <- function(abundance_data, metadata) {
  cat("Calculating beta diversity metrics...\n")

  beta_results <- list()

  for (level in names(abundance_data)) {
    cat("Processing", level, "level for beta diversity...\n")

    # Get relative abundance data (transposed for vegan)
    rel_matrix <- t(abundance_data[[level]]$relative)

    # Remove samples with zero total abundance
    sample_sums <- rowSums(rel_matrix)
    if (any(sample_sums == 0)) {
      cat("Removing", sum(sample_sums == 0), "samples with zero abundance\n")
      rel_matrix <- rel_matrix[sample_sums > 0, ]
    }

    if (nrow(rel_matrix) < 3) {
      cat("Warning: Not enough samples for beta diversity analysis at", level, "level\n")
      next
    }

    # Calculate distance matrices
    distances <- list(
      bray_curtis = vegdist(rel_matrix, method = "bray"),
      jaccard = vegdist(rel_matrix, method = "jaccard", binary = TRUE),
      euclidean = dist(rel_matrix, method = "euclidean")
    )

    # Perform PCoA (Principal Coordinates Analysis)
    pcoa_results <- list()
    for (dist_name in names(distances)) {
      pcoa <- cmdscale(distances[[dist_name]], k = min(nrow(rel_matrix) - 1, 10), eig = TRUE)

      # Calculate percentage of variance explained
      eigenvals <- pcoa$eig[pcoa$eig > 0]
      var_explained <- eigenvals / sum(eigenvals) * 100

      pcoa_df <- data.frame(
        SampleID = rownames(rel_matrix),
        PC1 = pcoa$points[, 1],
        PC2 = pcoa$points[, 2],
        stringsAsFactors = FALSE
      )

      # Add metadata
      pcoa_df <- merge(pcoa_df, metadata, by = "SampleID", all.x = TRUE)

      pcoa_results[[dist_name]] <- list(
        coordinates = pcoa_df,
        var_explained = var_explained[1:2]
      )
    }

    # Perform PERMANOVA if we have categories
    permanova_results <- NULL
    if ("Category" %in% colnames(metadata) && length(unique(metadata$Category)) > 1) {
      # Match metadata to distance matrix samples
      meta_subset <- metadata[metadata$SampleID %in% rownames(rel_matrix), ]
      meta_subset <- meta_subset[match(rownames(rel_matrix), meta_subset$SampleID), ]

      if (nrow(meta_subset) == nrow(rel_matrix)) {
        permanova_results <- list()
        for (dist_name in names(distances)) {
          perm_result <- adonis2(distances[[dist_name]] ~ Category, data = meta_subset, permutations = 999)
          permanova_results[[dist_name]] <- perm_result
        }
      }
    }

    beta_results[[level]] <- list(
      distances = distances,
      pcoa = pcoa_results,
      permanova = permanova_results
    )

    # Save distance matrices
    for (dist_name in names(distances)) {
      write.csv(
        as.matrix(distances[[dist_name]]),
        file.path(
          "results/diversity/tables",
          paste0(level, "_", dist_name, "_distance_matrix.csv")
        )
      )
    }
  }

  cat("Beta diversity calculation complete.\n")
  return(beta_results)
}

# Plot beta diversity
plot_beta_diversity <- function(beta_results) {
  cat("Creating beta diversity plots...\n")

  for (level in names(beta_results)) {
    beta_data <- beta_results[[level]]

    for (dist_name in names(beta_data$pcoa)) {
      pcoa_data <- beta_data$pcoa[[dist_name]]
      coords <- pcoa_data$coordinates
      var_exp <- pcoa_data$var_explained

      if (!"Category" %in% colnames(coords)) {
        coords$Category <- "Unknown"
      }

      # PCoA plot
      p <- ggplot(coords, aes(x = .data$PC1, y = .data$PC2, color = .data$Category)) +
        geom_point(size = 3, alpha = 0.7) +
        theme_minimal() +
        labs(
          title = paste(
            "PCoA -", stringr::str_to_title(level), "Level (",
            stringr::str_to_title(gsub("_", " ", dist_name)), "Distance)"
          ),
          x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
          y = paste0("PC2 (", round(var_exp[2], 1), "%)")
        ) +
        scale_color_brewer(type = "qual", palette = "Set1") +
        theme(legend.position = "bottom")

      # Add ellipses if we have categories
      if (length(unique(coords$Category)) > 1 && length(unique(coords$Category)) < nrow(coords)) {
        p <- p + stat_ellipse(type = "norm", level = 0.68, alpha = 0.3)
      }

      ggsave(
        file.path("results/diversity/plots", paste0(level, "_", dist_name, "_pcoa.png")),
        p,
        width = 10, height = 8, dpi = 300
      )

      # Distance heatmap
      dist_matrix <- as.matrix(beta_data$distances[[dist_name]])

      # Create annotation for heatmap
      annotation_df <- coords %>%
        select("SampleID", "Category")
      rownames(annotation_df) <- annotation_df$SampleID
      annotation_df <- annotation_df[, "Category", drop = FALSE]

      # Ensure row/column order matches
      common_samples <- intersect(rownames(dist_matrix), rownames(annotation_df))
      if (length(common_samples) > 1) {
        dist_matrix <- dist_matrix[common_samples, common_samples]
        annotation_df <- annotation_df[common_samples, , drop = FALSE]

        # Create color palette for categories
        n_categories <- length(unique(annotation_df$Category))
        colors <- RColorBrewer::brewer.pal(min(n_categories, 8), "Set2")
        names(colors) <- unique(annotation_df$Category)
        annotation_colors <- list(Category = colors)

        png(
          file.path("results/diversity/plots", paste0(level, "_", dist_name, "_heatmap.png")),
          width = 12, height = 10, units = "in", res = 300
        )

        pheatmap(dist_matrix,
          annotation_row = annotation_df,
          annotation_col = annotation_df,
          annotation_colors = annotation_colors,
          color = colorRampPalette(c("blue", "white", "red"))(100),
          main = paste(
            tools::toTitleCase(level), "Level -",
            tools::toTitleCase(gsub("_", " ", dist_name)), "Distance"
          ),
          fontsize = 8
        )

        dev.off()
      }
    }
  }

  cat("Beta diversity plots saved.\n")
}

# Generate statistical summaries
generate_statistical_summary <- function(alpha_results, beta_results) {
  cat("Generating statistical summaries...\n")

  for (level in names(alpha_results)) {
    alpha_data <- alpha_results[[level]]

    if ("Category" %in% colnames(alpha_data)) {
      # Fix variable binding warnings
      summary_stats <- alpha_data %>%
        group_by(.data$Category) %>%
        summarise(
          n_samples = n(),
          mean_observed_taxa = mean(.data$Observed_Taxa, na.rm = TRUE),
          sd_observed_taxa = sd(.data$Observed_Taxa, na.rm = TRUE),
          mean_shannon = mean(.data$Shannon, na.rm = TRUE),
          sd_shannon = sd(.data$Shannon, na.rm = TRUE),
          mean_simpson = mean(.data$Simpson, na.rm = TRUE),
          sd_simpson = sd(.data$Simpson, na.rm = TRUE),
          .groups = "drop"
        )

      # Save summary
      write.csv(summary_stats,
        file.path("results/diversity/tables", paste0(level, "_alpha_summary.csv")),
        row.names = FALSE
      )

      # Statistical tests (simplified)
      if (length(unique(alpha_data$Category)) > 1) {
        kw_shannon <- kruskal.test(Shannon ~ Category, data = alpha_data)

        # Save test results
        test_results <- data.frame(
          Test = "Kruskal-Wallis",
          Metric = "Shannon",
          Chi_squared = kw_shannon$statistic,
          P_value = kw_shannon$p.value
        )

        write.csv(test_results,
          file.path("results/diversity/tables", paste0(level, "_statistical_tests.csv")),
          row.names = FALSE
        )
      }
    }
  }

  # Beta diversity PERMANOVA results (simplified)
  if (length(beta_results) > 0) {
    permanova_summary <- data.frame()

    for (level in names(beta_results)) {
      if (!is.null(beta_results[[level]]$permanova)) {
        for (dist_name in names(beta_results[[level]]$permanova)) {
          perm_result <- beta_results[[level]]$permanova[[dist_name]]

          summary_row <- data.frame(
            Level = level,
            Distance = dist_name,
            F_statistic = perm_result$F[1],
            R_squared = perm_result$R2[1],
            P_value = perm_result$`Pr(>F)`[1]
          )

          permanova_summary <- rbind(permanova_summary, summary_row)
        }
      }
    }

    if (nrow(permanova_summary) > 0) {
      write.csv(permanova_summary,
        "results/diversity/tables/permanova_summary.csv",
        row.names = FALSE
      )
    }
  }

  cat("Statistical summaries saved.\n")
}

# Main function
main <- function() {
  cat("Starting Diversity Analysis...\n")
  cat("==========================================\n")

  # Create output directories
  create_directories()

  # Load data
  data <- load_data()
  abundance_data <- data$abundance_data
  metadata <- data$metadata

  # Alpha diversity analysis
  cat("\n--- Alpha Diversity Analysis ---\n")
  alpha_results <- calculate_alpha_diversity(abundance_data, metadata)
  plot_alpha_diversity(alpha_results)

  # Beta diversity analysis
  cat("\n--- Beta Diversity Analysis ---\n")
  beta_results <- calculate_beta_diversity(abundance_data, metadata)
  plot_beta_diversity(beta_results)

  # Statistical summaries
  cat("\n--- Statistical Analysis ---\n")
  generate_statistical_summary(alpha_results, beta_results)

  # Final summary
  cat("\n==========================================\n")
  cat("Diversity Analysis Complete!\n")
  cat("Results saved in: results/diversity/\n")
  cat("- Alpha diversity: results/diversity/alpha_diversity/\n")
  cat("- Beta diversity: results/diversity/beta_diversity/\n")
  cat("- Plots: results/diversity/plots/\n")
  cat("- Tables: results/diversity/tables/\n")
  cat("==========================================\n")
}

# Run main function
if (!interactive()) {
  main()
}
