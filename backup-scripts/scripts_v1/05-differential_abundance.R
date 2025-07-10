#!/usr/bin/env Rscript

# Differential Abundance Analysis Script for 16S rRNA Taxonomic Data
# Identifies taxa significantly different between sample categories
# Uses multiple statistical methods (DESeq2, LEfSe-style, ALDEx2)

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ALDEx2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
})

# Create output directories
create_directories <- function() {
  dirs <- c(
    "results/differential_abundance",
    "results/differential_abundance/plots",
    "results/differential_abundance/tables"
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
  cat("Loading taxonomic abundance data for differential analysis...\n")

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

  # Check if we have multiple categories for comparison
  n_categories <- length(unique(metadata$Category))
  if (n_categories < 2) {
    stop("Need at least 2 different categories for differential abundance analysis.")
  }
  cat("Found", n_categories, "categories:", paste(unique(metadata$Category), collapse = ", "), "\n")

  # Load abundance tables for different taxonomic levels
  levels <- c("species", "genus", "family", "order", "class", "phylum")
  abundance_data <- list()

  for (level in levels) {
    abundance_file <- file.path(taxonomy_dir, paste0(level, "_abundance_counts.csv"))
    relative_file <- file.path(taxonomy_dir, paste0(level, "_relative_abundance.csv"))

    if (file.exists(abundance_file) && file.exists(relative_file)) {
      counts <- read.csv(abundance_file, row.names = 1, check.names = FALSE)
      relative <- read.csv(relative_file, row.names = 1, check.names = FALSE)

      # Filter low abundance taxa (present in at least 10% of samples with >0.01% relative abundance)
      min_samples <- ceiling(ncol(relative) * 0.1)
      abundant_taxa <- rowSums(relative > 0.01) >= min_samples

      if (sum(abundant_taxa) > 0) {
        abundance_data[[level]] <- list(
          counts = counts[abundant_taxa, ],
          relative = relative[abundant_taxa, ]
        )
        cat(
          "Loaded", level, "level data:", sum(abundant_taxa), "taxa after filtering,",
          ncol(counts), "samples\n"
        )
      }
    }
  }

  if (length(abundance_data) == 0) {
    stop("No abundance data found after filtering. Please check taxonomy results.")
  }

  return(list(abundance_data = abundance_data, metadata = metadata))
}

# DESeq2 differential abundance analysis
run_deseq2_analysis <- function(count_data, metadata, level) {
  cat("Running DESeq2 analysis for", level, "level...\n")

  # Prepare count matrix (samples as columns, taxa as rows)
  count_matrix <- as.matrix(count_data)
  count_matrix <- round(count_matrix) # DESeq2 requires integer counts

  # Add pseudocount to handle zeros (common in microbiome data)
  count_matrix <- count_matrix + 1

  # Prepare metadata
  sample_data <- metadata[match(colnames(count_matrix), metadata$SampleID), ]
  rownames(sample_data) <- sample_data$SampleID
  sample_data$Category <- as.factor(sample_data$Category)

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~Category
  )

  # Filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  if (nrow(dds) == 0) {
    cat("Warning: No taxa remain after filtering for", level, "level\n")
    return(NULL)
  }

  # Run DESeq2
  dds <- DESeq(dds, quiet = TRUE)

  # Get results for all pairwise comparisons
  categories <- levels(sample_data$Category)
  deseq_results <- list()

  for (i in 1:(length(categories) - 1)) {
    for (j in (i + 1):length(categories)) {
      comparison_name <- paste0(categories[j], "_vs_", categories[i])

      res <- results(dds, contrast = c("Category", categories[j], categories[i]))
      res_df <- as.data.frame(res)
      res_df$taxon <- rownames(res_df)
      res_df$comparison <- comparison_name
      res_df$level <- level

      # Add significance categories
      res_df$significance <- "Not Significant"
      res_df$significance[res_df$padj < 0.05 & !is.na(res_df$padj)] <- "Significant (p < 0.05)"
      res_df$significance[res_df$padj < 0.01 & !is.na(res_df$padj)] <- "Highly Significant (p < 0.01)"

      # Add direction
      res_df$direction <- "No Change"
      res_df$direction[res_df$log2FoldChange > 1 & res_df$padj < 0.05 & !is.na(res_df$padj)] <- paste("Enriched in", categories[j])
      res_df$direction[res_df$log2FoldChange < -1 & res_df$padj < 0.05 & !is.na(res_df$padj)] <- paste("Enriched in", categories[i])

      deseq_results[[comparison_name]] <- res_df

      cat("  ", comparison_name, ": ", sum(res_df$padj < 0.05, na.rm = TRUE), " significant taxa\n")
    }
  }

  return(deseq_results)
}

# ALDEx2 differential abundance analysis - Fixed for multiple groups
run_aldex2_analysis <- function(count_data, metadata, level) {
  cat("Running ALDEx2 analysis for", level, "level...\n")

  # Prepare data
  count_matrix <- as.matrix(count_data)
  sample_data <- metadata[match(colnames(count_matrix), metadata$SampleID), ]
  conditions <- sample_data$Category

  # Get unique categories
  categories <- unique(conditions)

  # ALDEx2 works best with exactly 2 groups, so we'll do pairwise comparisons
  aldex_results <- list()

  for (i in 1:(length(categories) - 1)) {
    for (j in (i + 1):length(categories)) {
      comparison_name <- paste0(categories[j], "_vs_", categories[i])
      cat("  Running ALDEx2 for", comparison_name, "...\n")

      # Subset data for this comparison
      subset_samples <- conditions %in% c(categories[i], categories[j])
      subset_matrix <- count_matrix[, subset_samples]
      subset_conditions <- conditions[subset_samples]

      # Check if we have enough samples per group
      min_group_size <- min(table(subset_conditions))
      if (min_group_size < 3) {
        cat("    Warning: Small group sizes detected. Results may be unreliable.\n")
      }

      # Run ALDEx2
      tryCatch(
        {
          aldex_result <- aldex(subset_matrix, subset_conditions, mc.samples = 128, test = "t", effect = TRUE)

          # Prepare results dataframe
          aldex_df <- data.frame(
            taxon = rownames(aldex_result),
            effect_size = aldex_result$effect,
            p_value = aldex_result$we.ep,
            p_value_bh = aldex_result$we.eBH,
            comparison = comparison_name,
            level = level,
            stringsAsFactors = FALSE
          )

          # Add significance categories
          aldex_df$significance <- "Not Significant"
          aldex_df$significance[aldex_df$p_value_bh < 0.05] <- "Significant (p < 0.05)"
          aldex_df$significance[aldex_df$p_value_bh < 0.01] <- "Highly Significant (p < 0.01)"

          # Add direction based on effect size
          aldex_df$direction <- "No Change"
          aldex_df$direction[aldex_df$effect_size > 1 & aldex_df$p_value_bh < 0.05] <- paste("Enriched in", categories[j])
          aldex_df$direction[aldex_df$effect_size < -1 & aldex_df$p_value_bh < 0.05] <- paste("Enriched in", categories[i])

          aldex_results[[comparison_name]] <- aldex_df

          cat("    ", comparison_name, ":", sum(aldex_df$p_value_bh < 0.05, na.rm = TRUE), "significant taxa\n")
        },
        error = function(e) {
          cat("    Error in ALDEx2 for", comparison_name, ":", e$message, "\n")
        }
      )
    }
  }

  return(aldex_results)
}

# LEfSe-style analysis using Kruskal-Wallis and LDA
run_lefse_style_analysis <- function(relative_data, metadata, level) {
  cat("Running LEfSe-style analysis for", level, "level...\n")

  # Prepare data
  rel_matrix <- as.matrix(relative_data)
  sample_data <- metadata[match(colnames(rel_matrix), metadata$SampleID), ]

  lefse_results <- data.frame()

  for (taxon in rownames(rel_matrix)) {
    taxon_data <- rel_matrix[taxon, ]

    # Kruskal-Wallis test
    test_data <- data.frame(
      abundance = taxon_data,
      category = sample_data$Category
    )

    kw_test <- kruskal.test(abundance ~ category, data = test_data)

    if (kw_test$p.value < 0.05) {
      # Calculate LDA effect size (simplified)
      group_means <- aggregate(abundance ~ category, data = test_data, mean)
      max_mean <- max(group_means$abundance)
      max_group <- group_means$category[which.max(group_means$abundance)]

      # Simple effect size calculation
      effect_size <- log10(max_mean + 1e-6)

      result_row <- data.frame(
        taxon = taxon,
        p_value = kw_test$p.value,
        effect_size = effect_size,
        enriched_group = as.character(max_group),
        level = level,
        stringsAsFactors = FALSE
      )

      lefse_results <- rbind(lefse_results, result_row)
    }
  }

  # Apply multiple testing correction
  if (nrow(lefse_results) > 0) {
    lefse_results$p_value_bh <- p.adjust(lefse_results$p_value, method = "BH")

    # Add significance categories
    lefse_results$significance <- "Not Significant"
    lefse_results$significance[lefse_results$p_value_bh < 0.05] <- "Significant (p < 0.05)"
    lefse_results$significance[lefse_results$p_value_bh < 0.01] <- "Highly Significant (p < 0.01)"

    # Filter by significance and effect size
    lefse_results <- lefse_results[lefse_results$p_value_bh < 0.05 & lefse_results$effect_size > 2, ]

    cat("  LEfSe-style:", nrow(lefse_results), "significant taxa\n")
  } else {
    cat("  LEfSe-style: No significant taxa found\n")
  }

  return(lefse_results)
}

# Create volcano plots for DESeq2 results
plot_volcano <- function(deseq_results, level) {
  cat("Creating volcano plots for", level, "level...\n")

  for (comparison in names(deseq_results)) {
    res_df <- deseq_results[[comparison]]

    # Remove rows with NA values
    res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

    if (nrow(res_df) == 0) next

    # Create volcano plot
    p <- ggplot(res_df, aes(x = .data$log2FoldChange, y = -log10(.data$padj))) +
      geom_point(aes(color = .data$significance), alpha = 0.7, size = 2) +
      scale_color_manual(values = c(
        "Not Significant" = "grey70",
        "Significant (p < 0.05)" = "orange",
        "Highly Significant (p < 0.01)" = "red"
      )) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      theme_minimal() +
      labs(
        title = paste("Volcano Plot -", tools::toTitleCase(level), "Level"),
        subtitle = comparison,
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value",
        color = "Significance"
      ) +
      theme(legend.position = "bottom")

    ggsave(
      file.path(
        "results/differential_abundance/plots",
        paste0(level, "_", comparison, "_volcano.png")
      ),
      p,
      width = 10, height = 8, dpi = 300
    )
  }
}

# Create effect size plots for ALDEx2
plot_aldex2_effects <- function(aldex_results, level) {
  cat("Creating ALDEx2 effect plots for", level, "level...\n")

  for (comparison in names(aldex_results)) {
    aldex_df <- aldex_results[[comparison]]

    # Effect size plot
    p <- ggplot(aldex_df, aes(x = .data$effect_size, y = -log10(.data$p_value_bh))) +
      geom_point(aes(color = .data$significance), alpha = 0.7, size = 2) +
      scale_color_manual(values = c(
        "Not Significant" = "grey70",
        "Significant (p < 0.05)" = "orange",
        "Highly Significant (p < 0.01)" = "red"
      )) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      theme_minimal() +
      labs(
        title = paste("ALDEx2 Effect Plot -", tools::toTitleCase(level), "Level"),
        subtitle = comparison,
        x = "Effect Size",
        y = "-Log10 Adjusted P-value",
        color = "Significance"
      ) +
      theme(legend.position = "bottom")

    ggsave(
      file.path(
        "results/differential_abundance/plots",
        paste0(level, "_", comparison, "_aldex2_effect.png")
      ),
      p,
      width = 10, height = 8, dpi = 300
    )
  }
}

# Create heatmaps of significant taxa
create_heatmaps <- function(abundance_data, metadata, significant_taxa, level, method) {
  cat("Creating heatmap for", level, "level (", method, ")...\n")

  if (length(significant_taxa) == 0) {
    cat("No significant taxa to plot for", level, "level\n")
    return(NULL)
  }

  # Get relative abundance data for significant taxa
  rel_data <- abundance_data$relative[rownames(abundance_data$relative) %in% significant_taxa, ]

  if (nrow(rel_data) == 0) {
    return(NULL)
  }
  if (nrow(rel_data) == 0) {
    return(NULL)
  }

  # Log transform for better visualization
  log_data <- log10(rel_data + 1e-6)

  # Prepare sample annotation
  sample_data <- metadata[match(colnames(log_data), metadata$SampleID), ]
  annotation_df <- data.frame(
    Category = sample_data$Category,
    row.names = colnames(log_data)
  )

  # Color palette for categories
  n_categories <- length(unique(annotation_df$Category))
  colors <- RColorBrewer::brewer.pal(min(max(n_categories, 3), 8), "Set2")
  names(colors) <- unique(annotation_df$Category)
  annotation_colors <- list(Category = colors)

  # Create heatmap
  png(
    file.path(
      "results/differential_abundance/plots", # Changed from "png" to "plots"
      paste0(level, "_", method, "_significant_taxa_heatmap.png")
    ),
    width = 12, height = max(8, nrow(log_data) * 0.3), units = "in", res = 300
  )

  pheatmap(log_data,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    main = paste(
      "Significant Taxa -", tools::toTitleCase(level),
      "Level (", tools::toTitleCase(method), ")"
    ),
    fontsize_row = 8,
    fontsize_col = 8
  )

  dev.off()
}

# Save all results to files
save_results <- function(all_results) {
  cat("Saving differential abundance results...\n")

  for (level in names(all_results)) {
    level_results <- all_results[[level]]

    # Save DESeq2 results
    if (!is.null(level_results$deseq2)) {
      for (comparison in names(level_results$deseq2)) {
        write.csv(level_results$deseq2[[comparison]],
          file.path(
            "results/differential_abundance/tables",
            paste0(level, "_deseq2_", comparison, ".csv")
          ),
          row.names = FALSE
        )
      }
    }

    # Save ALDEx2 results
    if (!is.null(level_results$aldex2)) {
      for (comparison in names(level_results$aldex2)) {
        write.csv(level_results$aldex2[[comparison]],
          file.path(
            "results/differential_abundance/tables",
            paste0(level, "_aldex2_", comparison, ".csv")
          ),
          row.names = FALSE
        )
      }
    }

    # Save LEfSe-style results
    if (!is.null(level_results$lefse) && nrow(level_results$lefse) > 0) {
      write.csv(level_results$lefse,
        file.path(
          "results/differential_abundance/tables",
          paste0(level, "_lefse_style.csv")
        ),
        row.names = FALSE
      )
    }
  }

  cat("All results saved to results/differential_abundance/tables/\n")
}

# Main function
main <- function() {
  cat("Starting Differential Abundance Analysis...\n")
  cat("==========================================\n")

  # Create output directories
  create_directories()

  # Load data
  data <- load_data()
  abundance_data <- data$abundance_data
  metadata <- data$metadata

  # Run differential abundance analysis for each taxonomic level
  all_results <- list()

  for (level in names(abundance_data)) {
    cat("\n--- Processing", tools::toTitleCase(level), "Level ---\n")

    count_data <- abundance_data[[level]]$counts
    relative_data <- abundance_data[[level]]$relative

    level_results <- list()

    # DESeq2 analysis
    tryCatch(
      {
        deseq2_results <- run_deseq2_analysis(count_data, metadata, level)
        level_results$deseq2 <- deseq2_results
      },
      error = function(e) {
        cat("Error in DESeq2 analysis for", level, ":", e$message, "\n")
        level_results$deseq2 <- NULL
      }
    )

    # ALDEx2 analysis
    aldex2_results <- run_aldex2_analysis(count_data, metadata, level)
    level_results$aldex2 <- aldex2_results

    # LEfSe-style analysis
    tryCatch(
      {
        lefse_results <- run_lefse_style_analysis(relative_data, metadata, level)
        level_results$lefse <- lefse_results
      },
      error = function(e) {
        cat("Error in LEfSe analysis for", level, ":", e$message, "\n")
        level_results$lefse <- NULL
      }
    )

    all_results[[level]] <- level_results

    # Create plots
    if (!is.null(deseq2_results)) {
      plot_volcano(deseq2_results, level)
    }

    if (!is.null(aldex2_results)) {
      plot_aldex2_effects(aldex2_results, level)
    }

    # Create heatmaps for significant taxa
    if (!is.null(deseq2_results)) {
      for (comparison in names(deseq2_results)) {
        sig_taxa <- deseq2_results[[comparison]]$taxon[deseq2_results[[comparison]]$padj < 0.05 &
          !is.na(deseq2_results[[comparison]]$padj)]
        if (length(sig_taxa) > 0) {
          create_heatmaps(
            abundance_data[[level]], metadata, sig_taxa, level,
            paste0("deseq2_", comparison)
          )
        }
      }
    }

    if (!is.null(aldex2_results)) {
      for (comparison in names(aldex2_results)) {
        sig_taxa <- aldex2_results[[comparison]]$taxon[aldex2_results[[comparison]]$p_value_bh < 0.05]
        if (length(sig_taxa) > 0) {
          create_heatmaps(
            abundance_data[[level]], metadata, sig_taxa, level,
            paste0("aldex2_", comparison)
          )
        }
      }
    }

    if (!is.null(lefse_results) && nrow(lefse_results) > 0) {
      sig_taxa <- lefse_results$taxon
      if (length(sig_taxa) > 0) {
        create_heatmaps(abundance_data[[level]], metadata, sig_taxa, level, "lefse")
      }
    }
  }

  # Save all results
  save_results(all_results)

  # Final summary
  cat("\n==========================================\n")
  cat("Differential Abundance Analysis Complete!\n")
  cat("Results saved in: results/differential_abundance/\n")
  cat("- Plots: results/differential_abundance/plots/\n")
  # cat("- Heatmaps: results/differential_abundance/heatmaps/\n")
  cat("- Tables: results/differential_abundance/tables/\n")
  cat("==========================================\n")
}

# Run main function
if (!interactive()) {
  main()
}
