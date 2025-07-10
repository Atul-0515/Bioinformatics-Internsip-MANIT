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
    "results/differential_abundance/tables",
    "results/differential_abundance/reports"
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

  count_matrix <- as.matrix(count_data)
  count_matrix <- round(count_matrix)

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

run_aldex2_analysis <- function(count_data, metadata, level) {
  cat("Running ALDEx2 analysis for", level, "level...\n")

  # Prepare data
  count_matrix <- as.matrix(count_data)
  sample_data <- metadata[match(colnames(count_matrix), metadata$SampleID), ]
  conditions <- sample_data$Category

  # Get unique categories
  categories <- unique(conditions)

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
      # Calculate LDA effect size
      group_means <- aggregate(abundance ~ category, data = test_data, mean)
      max_mean <- max(group_means$abundance)
      max_group <- group_means$category[which.max(group_means$abundance)]

      # Simple effect size calculation
      effect_size <- log10(max_mean + 1e-6)

      result_row <- data.frame(
        taxon = taxon,
        p_value = kw_test$p.value,
        effect_size = effect_size,
        max_mean_abundance = max_mean,
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
    lefse_results <- lefse_results[lefse_results$p_value_bh < 0.05 & lefse_results$max_mean_abundance > 0.01, ]

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
      "results/differential_abundance/plots",
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

# Initialize summary report
initialize_summary_report <- function() {
  report_dir <- "results/differential_abundance/reports"
  if (!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE)
  }

  report_file <- file.path(report_dir, paste0(
    "differential_abundance_summary_",
    format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"
  ))

  # Write header
  writeLines(c(
    "DIFFERENTIAL ABUNDANCE ANALYSIS SUMMARY REPORT",
    "===============================================",
    paste("Analysis Date:", Sys.time()),
    paste("R Version:", R.version.string),
    "",
    "OVERVIEW",
    "--------"
  ), report_file)

  return(report_file)
}

# Log data summary
log_data_summary <- function(report_file, abundance_data, metadata) {
  summary_lines <- c(
    paste("Total Samples:", nrow(metadata)),
    paste("Sample Categories:", paste(unique(metadata$Category), collapse = ", ")),
    paste("Number of Categories:", length(unique(metadata$Category))),
    "",
    "TAXONOMIC LEVELS PROCESSED:",
    "-------------------------"
  )

  for (level in names(abundance_data)) {
    n_taxa <- nrow(abundance_data[[level]]$counts)
    summary_lines <- c(
      summary_lines,
      paste("  -", tools::toTitleCase(level), ":", n_taxa, "taxa")
    )
  }

  summary_lines <- c(summary_lines, "")

  cat(paste(summary_lines, collapse = "\n"), "\n", file = report_file, append = TRUE)
}

# Log method results summary
log_method_summary <- function(report_file, all_results) {
  method_summary <- c(
    "METHOD RESULTS SUMMARY",
    "=====================",
    ""
  )

  for (level in names(all_results)) {
    method_summary <- c(
      method_summary,
      paste("TAXONOMIC LEVEL:", toupper(level)),
      paste(rep("-", nchar(level) + 17), collapse = "")
    )

    level_results <- all_results[[level]]

    # DESeq2 summary
    if (!is.null(level_results$deseq2)) {
      method_summary <- c(method_summary, "DESeq2 Analysis:")
      for (comparison in names(level_results$deseq2)) {
        res <- level_results$deseq2[[comparison]]
        n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
        n_high_sig <- sum(res$padj < 0.01, na.rm = TRUE)
        method_summary <- c(
          method_summary,
          paste("  ", comparison, ":", n_sig, "significant taxa (", n_high_sig, "highly significant)")
        )
      }
    } else {
      method_summary <- c(method_summary, "DESeq2 Analysis: FAILED")
    }

    # ALDEx2 summary
    if (!is.null(level_results$aldex2)) {
      method_summary <- c(method_summary, "ALDEx2 Analysis:")
      for (comparison in names(level_results$aldex2)) {
        res <- level_results$aldex2[[comparison]]
        n_sig <- sum(res$p_value_bh < 0.05, na.rm = TRUE)
        n_high_sig <- sum(res$p_value_bh < 0.01, na.rm = TRUE)
        method_summary <- c(
          method_summary,
          paste("  ", comparison, ":", n_sig, "significant taxa (", n_high_sig, "highly significant)")
        )
      }
    } else {
      method_summary <- c(method_summary, "ALDEx2 Analysis: No results or failed")
    }

    # LEfSe summary
    if (!is.null(level_results$lefse) && nrow(level_results$lefse) > 0) {
      n_sig <- nrow(level_results$lefse)
      method_summary <- c(
        method_summary,
        paste("LEfSe-style Analysis:", n_sig, "significant biomarkers")
      )
    } else {
      method_summary <- c(method_summary, "LEfSe-style Analysis: No significant biomarkers")
    }

    method_summary <- c(method_summary, "")
  }

  cat(paste(method_summary, collapse = "\n"), "\n", file = report_file, append = TRUE)
}

# Log significant taxa details
log_significant_taxa <- function(report_file, all_results) {
  taxa_details <- c(
    "SIGNIFICANT TAXA DETAILS",
    "=======================",
    ""
  )

  for (level in names(all_results)) {
    level_results <- all_results[[level]]
    has_significant <- FALSE

    level_section <- c(
      paste("LEVEL:", toupper(level)),
      paste(rep("-", nchar(level) + 7), collapse = "")
    )

    # Top significant taxa from DESeq2
    if (!is.null(level_results$deseq2)) {
      for (comparison in names(level_results$deseq2)) {
        res <- level_results$deseq2[[comparison]]
        sig_taxa <- res[res$padj < 0.05 & !is.na(res$padj), ]

        if (nrow(sig_taxa) > 0) {
          has_significant <- TRUE
          sig_taxa <- sig_taxa[order(sig_taxa$padj), ]
          top_taxa <- head(sig_taxa, 5)

          level_section <- c(
            level_section,
            paste("DESeq2 -", comparison, "(Top 5):")
          )

          for (i in seq_len(nrow(top_taxa))) {
            direction <- ifelse(top_taxa$log2FoldChange[i] > 0, "↑", "↓")
            level_section <- c(
              level_section,
              sprintf(
                "  %s %s (LFC: %.2f, p-adj: %.2e)",
                direction, top_taxa$taxon[i],
                top_taxa$log2FoldChange[i], top_taxa$padj[i]
              )
            )
          }
          level_section <- c(level_section, "")
        }
      }
    }

    # Top significant taxa from ALDEx2
    if (!is.null(level_results$aldex2)) {
      for (comparison in names(level_results$aldex2)) {
        res <- level_results$aldex2[[comparison]]
        sig_taxa <- res[res$p_value_bh < 0.05, ]

        if (nrow(sig_taxa) > 0) {
          has_significant <- TRUE
          sig_taxa <- sig_taxa[order(sig_taxa$p_value_bh), ]
          top_taxa <- head(sig_taxa, 5)

          level_section <- c(
            level_section,
            paste("ALDEx2 -", comparison, "(Top 5):")
          )

          for (i in seq_len(nrow(top_taxa))) {
            direction <- ifelse(top_taxa$effect_size[i] > 0, "↑", "↓")
            level_section <- c(
              level_section,
              sprintf(
                "  %s %s (Effect: %.2f, p-adj: %.2e)",
                direction, top_taxa$taxon[i],
                top_taxa$effect_size[i], top_taxa$p_value_bh[i]
              )
            )
          }
          level_section <- c(level_section, "")
        }
      }
    }

    if (!is.null(level_results$lefse) && nrow(level_results$lefse) > 0) {
      res <- level_results$lefse
      has_significant <- TRUE
      res <- res[order(res$p_value_bh), ]
      top_taxa <- head(res, 5)

      level_section <- c(
        level_section,
        "LEfSe-style - Top 5 Biomarkers (by p-value):"
      )

      for (i in seq_len(nrow(top_taxa))) {
        level_section <- c(
          level_section,
          sprintf(
            "  - %s (Enriched in: %s, p-adj: %.2e)",
            top_taxa$taxon[i],
            top_taxa$enriched_group[i],
            top_taxa$p_value_bh[i]
          )
        )
      }
      level_section <- c(level_section, "")
    }
    if (has_significant) {
      taxa_details <- c(taxa_details, level_section)
    }
  }

  if (length(taxa_details) == 3) {
    taxa_details <- c(taxa_details, "No significant taxa found across all levels and methods.")
  }

  cat(paste(taxa_details, collapse = "\n"), "\n", file = report_file, append = TRUE)
}

# Finalize report
finalize_report <- function(report_file, all_results) {
  # Count total significant findings
  total_comparisons <- 0
  total_significant <- 0

  for (level in names(all_results)) {
    level_results <- all_results[[level]]

    if (!is.null(level_results$deseq2)) {
      for (comparison in names(level_results$deseq2)) {
        total_comparisons <- total_comparisons + 1
        res <- level_results$deseq2[[comparison]]
        total_significant <- total_significant + sum(res$padj < 0.05, na.rm = TRUE)
      }
    }
  }

  summary_footer <- c(
    "",
    "ANALYSIS SUMMARY",
    "===============",
    paste("Total taxonomic levels analyzed:", length(all_results)),
    paste("Total comparisons performed:", total_comparisons),
    paste("Total significant taxa found:", total_significant),
    "",
    "FILES GENERATED:",
    "- Volcano plots: results/differential_abundance/plots/*_volcano.png",
    "- Effect plots: results/differential_abundance/plots/*_aldex2_effect.png",
    "- Heatmaps: results/differential_abundance/plots/*_heatmap.png",
    "- Data tables: results/differential_abundance/tables/*.csv",
    "",
    paste("Analysis completed:", Sys.time()),
    "==============================================="
  )

  cat(paste(summary_footer, collapse = "\n"), "\n", file = report_file, append = TRUE)

  cat("Summary report saved to:", report_file, "\n")
}

# Main function
main <- function() {
  cat("Starting Differential Abundance Analysis...\n")
  cat("==========================================\n")

  # Create output directories
  create_directories()

  report_file <- initialize_summary_report()

  # Load data
  data <- load_data()
  abundance_data <- data$abundance_data
  metadata <- data$metadata

  log_data_summary(report_file, abundance_data, metadata)

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

  log_method_summary(report_file, all_results)

  log_significant_taxa(report_file, all_results)

  finalize_report(report_file, all_results)

  # Final summary
  cat("\n==========================================\n")
  cat("Differential Abundance Analysis Complete!\n")
  cat("Results saved in: results/differential_abundance/\n")
  cat("- Plots: results/differential_abundance/plots/\n")
  cat("- Reports: results/differential_abundance/reports/\n")
  cat("- Tables: results/differential_abundance/tables/\n")
  cat("==========================================\n")
}

if (!interactive()) {
  main()
}
