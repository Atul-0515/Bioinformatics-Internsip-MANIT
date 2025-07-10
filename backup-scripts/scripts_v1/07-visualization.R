#!/usr/bin/env Rscript

# ==============================================================================
# 07-visualization.R
# Comprehensive Visualization Script for Bioinformatics Project
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(RColorBrewer)
  library(pheatmap)
  library(vegan)
  library(plotly)
  library(gridExtra)
  library(ComplexHeatmap)
  library(ggpubr)
  library(scales)
  library(stringr)
  library(ggpubr)
  library(cowplot)
  library(scales)
  library(ggrepel)
  library(viridis)
  library(tibble)
})

# Set up plotting parameters
theme_set(theme_bw())
options(
  ggplot2.discrete.colour = function(...) scale_colour_brewer(..., palette = "Set2"),
  ggplot2.discrete.fill = function(...) scale_fill_brewer(..., palette = "Set2")
)

# Create output directory for plots
dir.create("results/final_visualizations", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# SECTION 1: QUALITY CONTROL SUMMARY PLOTS
# ==============================================================================

cat("Generating Quality Control summary plots...\n")

# Read trimming summary data
if (file.exists("results/qc/trimming_summary.csv")) {
  qc_data <- read.csv("results/qc/trimming_summary.csv")

  # 1. Read count comparison (before vs after trimming)
  read_comparison <- qc_data %>%
    select(sample, total_reads_before, total_reads_after) %>%
    pivot_longer(
      cols = c(total_reads_before, total_reads_after),
      names_to = "stage", values_to = "read_count"
    ) %>%
    mutate(stage = case_when(
      stage == "total_reads_before" ~ "Before Trimming",
      stage == "total_reads_after" ~ "After Trimming"
    ))

  p1 <- ggplot(read_comparison, aes(x = sample, y = read_count, fill = stage)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "Read Counts Before and After Quality Trimming",
      x = "Sample", y = "Number of Reads", fill = "Stage"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = scales::comma)

  # 2. Read retention percentage
  p2 <- ggplot(qc_data, aes(
    x = reorder(sample, reads_retained_percent),
    y = reads_retained_percent
  )) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
    labs(
      title = "Read Retention After Quality Filtering",
      x = "Sample", y = "Reads Retained (%)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 100)

  # 3. Q30 quality improvement
  q30_comparison <- qc_data %>%
    select(sample, q30_rate_before, q30_rate_after) %>%
    pivot_longer(
      cols = c(q30_rate_before, q30_rate_after),
      names_to = "stage", values_to = "q30_rate"
    ) %>%
    mutate(stage = case_when(
      stage == "q30_rate_before" ~ "Before Trimming",
      stage == "q30_rate_after" ~ "After Trimming"
    ))

  p3 <- ggplot(q30_comparison, aes(x = sample, y = q30_rate * 100, fill = stage)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "Q30 Quality Rate Improvement",
      x = "Sample", y = "Q30 Rate (%)", fill = "Stage"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # 4. Summary statistics table plot
  summary_stats <- qc_data %>%
    summarise(
      mean_retention = round(mean(reads_retained_percent), 2),
      median_retention = round(median(reads_retained_percent), 2),
      min_retention = round(min(reads_retained_percent), 2),
      max_retention = round(max(reads_retained_percent), 2),
      total_reads_input = sum(total_reads_before),
      total_reads_output = sum(total_reads_after)
    )

  # Create summary table visualization
  p4 <- ggplot() +
    theme_void() +
    labs(title = "Quality Control Summary Statistics") +
    annotation_custom(
      tableGrob(
        data.frame(
          Metric = c(
            "Mean Retention (%)", "Median Retention (%)",
            "Min Retention (%)", "Max Retention (%)",
            "Total Input Reads", "Total Output Reads"
          ),
          Value = c(
            summary_stats$mean_retention, summary_stats$median_retention,
            summary_stats$min_retention, summary_stats$max_retention,
            format(summary_stats$total_reads_input, big.mark = ","),
            format(summary_stats$total_reads_output, big.mark = ",")
          )
        ),
        theme = ttheme_default()
      )
    )

  # Combine plots
  qc_combined <- plot_grid(p1, p2, p3, p4,
    labels = c("A", "B", "C", "D"),
    ncol = 2, nrow = 2
  )

  # Save individual plots
  ggsave("results/final_visualizations/read_count_comparison.png", p1,
    width = 12, height = 6, dpi = 300
  )
  ggsave("results/final_visualizations/read_retention.png", p2,
    width = 10, height = 6, dpi = 300
  )
  ggsave("results/final_visualizations/q30_improvement.png", p3,
    width = 12, height = 6, dpi = 300
  )

  # Save combined plot
  ggsave("results/final_visualizations/qc_summary_combined.png", qc_combined,
    width = 16, height = 12, dpi = 300
  )

  cat("Quality Control plots saved successfully!\n")
} else {
  cat("Warning: trimming_summary.csv not found. Skipping QC visualization.\n")
}


# ==============================================================================
# SECTION 2: TAXONOMIC COMPOSITION PLOTS
# ==============================================================================

cat("Generating Taxonomic Composition plots...\n")

# Load sample categories
sample_categories <- read.csv("results/taxonomy/sample_categories.csv", stringsAsFactors = FALSE)

# Function to load and prepare taxonomic data
load_taxonomic_data <- function(level) {
  file_path <- paste0("results/taxonomy/merged_tables/", level, "_relative_abundance.csv")
  if (file.exists(file_path)) {
    data <- read.csv(file_path, row.names = 1, check.names = FALSE)
    return(data)
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
}

# Load data for different taxonomic levels
phylum_data <- load_taxonomic_data("phylum")
genus_data <- load_taxonomic_data("genus")
species_data <- load_taxonomic_data("species")

# 1. RELATIVE ABUNDANCE BAR PLOTS
if (!is.null(phylum_data)) {
  # Prepare data for plotting
  phylum_long <- phylum_data %>%
    tibble::rownames_to_column("Taxon") %>%
    pivot_longer(-Taxon, names_to = "SampleID", values_to = "Abundance") %>%
    left_join(sample_categories, by = "SampleID")

  # Keep only top 10 most abundant phyla, group others as "Other"
  top_phyla <- phylum_long %>%
    group_by(Taxon) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
    arrange(desc(mean_abundance)) %>%
    slice_head(n = 10) %>%
    pull(Taxon)

  phylum_plot_data <- phylum_long %>%
    mutate(Taxon_grouped = ifelse(Taxon %in% top_phyla, Taxon, "Other")) %>%
    group_by(SampleID, Category, Taxon_grouped) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")

  # Create stacked bar plot
  p1 <- ggplot(phylum_plot_data, aes(x = SampleID, y = Abundance, fill = Taxon_grouped)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Category, scales = "free_x") +
    labs(
      title = "Phylum-level Relative Abundance",
      x = "Sample ID", y = "Relative Abundance (%)",
      fill = "Phylum"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    scale_fill_brewer(type = "qual", palette = "Set3")

  ggsave("results/final_visualizations/phylum_abundance_barplot.png",
    plot = p1, width = 12, height = 8, dpi = 300
  )
}

# 2. TAXONOMIC COMPOSITION HEATMAP
if (!is.null(genus_data)) {
  # Filter for top 30 most abundant genera
  genus_means <- rowMeans(genus_data)
  top_genera <- names(sort(genus_means, decreasing = TRUE)[1:30])
  genus_subset <- as.matrix(genus_data[top_genera, ])

  # Add sample categories as annotation
  sample_annotation <- sample_categories %>%
    tibble::column_to_rownames("SampleID")

  # Create heatmap using base pheatmap (not ComplexHeatmap)
  png("results/final_visualizations/genus_heatmap.png", width = 12 * 100, height = 10 * 100, res = 300)
  pheatmap(log10(genus_subset + 0.1),
    annotation_col = sample_annotation,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    scale = "row",
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Top 30 Genera Heatmap (Log10 Relative Abundance)"
  )
  invisible(dev.off())
}

# 3. DIVERSITY OVERVIEW PLOTS
if (!is.null(species_data)) {
  # Calculate alpha diversity metrics
  diversity_metrics <- data.frame(
    SampleID = colnames(species_data),
    Shannon = apply(species_data, 2, function(x) diversity(x, index = "shannon")),
    Simpson = apply(species_data, 2, function(x) diversity(x, index = "simpson")),
    Richness = apply(species_data, 2, function(x) sum(x > 0))
  ) %>%
    left_join(sample_categories, by = "SampleID")

  # Shannon diversity boxplot
  p2 <- ggplot(diversity_metrics, aes(x = Category, y = Shannon, fill = Category)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(
      title = "Shannon Diversity by Category",
      x = "Category", y = "Shannon Diversity Index"
    ) +
    theme(legend.position = "none") +
    stat_compare_means()

  # Species richness boxplot
  p3 <- ggplot(diversity_metrics, aes(x = Category, y = Richness, fill = Category)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(
      title = "Species Richness by Category",
      x = "Category", y = "Number of Species"
    ) +
    theme(legend.position = "none") +
    stat_compare_means()

  # Combine diversity plots
  diversity_combined <- plot_grid(p2, p3, ncol = 2)
  ggsave("results/final_visualizations/diversity_overview.png",
    plot = diversity_combined, width = 12, height = 6, dpi = 300
  )

  # Beta diversity PCoA
  species_t <- t(species_data)
  bray_dist <- vegdist(species_t, method = "bray")
  pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)

  pcoa_data <- data.frame(
    SampleID = rownames(pcoa_result$points),
    PC1 = pcoa_result$points[, 1],
    PC2 = pcoa_result$points[, 2]
  ) %>%
    left_join(sample_categories, by = "SampleID")

  variance_explained <- round(pcoa_result$eig[1:2] / sum(pcoa_result$eig) * 100, 1)

  p4 <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = Category)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = "Beta Diversity (Bray-Curtis PCoA)",
      x = paste0("PC1 (", variance_explained[1], "%)"),
      y = paste0("PC2 (", variance_explained[2], "%)")
    ) +
    theme(legend.position = "bottom") +
    stat_ellipse(type = "norm", level = 0.95)

  ggsave("results/final_visualizations/beta_diversity_pcoa.png",
    plot = p4, width = 10, height = 8, dpi = 300
  )
}

cat("Taxonomic composition plots completed!\n")


# ==============================================================================
# SECTION 3: DIVERSITY ANALYSIS PLOTS
# ==============================================================================

cat("Generating Diversity Analysis plots...\n")

# Load diversity analysis results
load_diversity_results <- function() {
  # Load alpha diversity results
  alpha_files <- list.files("results/diversity/tables", pattern = "*_alpha_diversity.csv", full.names = TRUE)
  alpha_results <- list()

  for (file in alpha_files) {
    level <- gsub("_alpha_diversity.csv", "", basename(file))
    alpha_results[[level]] <- read.csv(file, stringsAsFactors = FALSE)
  }

  # Load beta diversity PCoA coordinates (we'll recreate from distance matrices)
  dist_files <- list.files("results/diversity/tables", pattern = "*_distance_matrix.csv", full.names = TRUE)
  beta_results <- list()

  for (file in dist_files) {
    filename <- basename(file)
    parts <- strsplit(gsub("_distance_matrix.csv", "", filename), "_")[[1]]
    level <- parts[1]
    dist_method <- paste(parts[-1], collapse = "_")

    if (!level %in% names(beta_results)) {
      beta_results[[level]] <- list()
    }

    dist_matrix <- read.csv(file, row.names = 1, check.names = FALSE)
    beta_results[[level]][[dist_method]] <- as.dist(dist_matrix)
  }

  return(list(alpha = alpha_results, beta = beta_results))
}

# 1. Alpha Diversity Comparative Boxplots
create_alpha_diversity_comparison <- function(alpha_results) {
  cat("Creating alpha diversity comparison plots...\n")

  # Combine all levels for comparison
  combined_data <- data.frame()

  for (level in names(alpha_results)) {
    if ("Category" %in% colnames(alpha_results[[level]])) {
      temp_data <- alpha_results[[level]] %>%
        select("SampleID", "Category", "Shannon", "Simpson", "Observed_Taxa") %>%
        mutate(Taxonomic_Level = str_to_title(level))

      combined_data <- rbind(combined_data, temp_data)
    }
  }

  if (nrow(combined_data) > 0) {
    # Shannon diversity comparison across levels
    p_shannon <- ggplot(combined_data, aes(x = .data$Taxonomic_Level, y = .data$Shannon, fill = .data$Category)) +
      geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
      geom_point(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.6) +
      scale_fill_brewer(type = "qual", palette = "Set2") +
      labs(
        title = "Shannon Diversity Across Taxonomic Levels",
        x = "Taxonomic Level", y = "Shannon Index"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Observed taxa comparison
    p_obs <- ggplot(combined_data, aes(x = .data$Taxonomic_Level, y = .data$Observed_Taxa, fill = .data$Category)) +
      geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
      geom_point(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.6) +
      scale_fill_brewer(type = "qual", palette = "Set2") +
      labs(
        title = "Observed Taxa Across Taxonomic Levels",
        x = "Taxonomic Level", y = "Number of Observed Taxa"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Combined plot
    combined_plot <- plot_grid(p_shannon, p_obs, ncol = 1, align = "v")

    ggsave("results/final_visualizations/alpha_diversity_comparison.png",
      combined_plot,
      width = 12, height = 10, dpi = 300
    )
  }
}

# 2. Beta Diversity PCoA Summary Plot
create_beta_diversity_summary <- function(beta_results, alpha_results) {
  cat("Creating beta diversity summary plots...\n")

  # Load sample metadata
  metadata_file <- "results/taxonomy/sample_categories.csv"
  if (file.exists(metadata_file)) {
    metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

    # Create PCoA plots for genus level (most informative)
    if ("genus" %in% names(beta_results)) {
      pcoa_plots <- list()

      for (dist_method in names(beta_results$genus)) {
        dist_matrix <- beta_results$genus[[dist_method]]

        # Perform PCoA
        pcoa <- cmdscale(dist_matrix, k = 2, eig = TRUE)
        eigenvals <- pcoa$eig[pcoa$eig > 0]
        var_explained <- eigenvals / sum(eigenvals) * 100

        # Create data frame
        pcoa_df <- data.frame(
          SampleID = rownames(pcoa$points),
          PC1 = pcoa$points[, 1],
          PC2 = pcoa$points[, 2],
          stringsAsFactors = FALSE
        )

        # Add metadata
        pcoa_df <- merge(pcoa_df, metadata, by = "SampleID", all.x = TRUE)

        # Create plot
        p <- ggplot(pcoa_df, aes(x = .data$PC1, y = .data$PC2, color = .data$Category)) +
          geom_point(size = 3, alpha = 0.7) +
          stat_ellipse(type = "norm", level = 0.68, alpha = 0.3) +
          scale_color_brewer(type = "qual", palette = "Set1") +
          labs(
            title = paste("PCoA -", str_to_title(gsub("_", " ", dist_method))),
            x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
            y = paste0("PC2 (", round(var_explained[2], 1), "%)")
          ) +
          theme_minimal() +
          theme(legend.position = "bottom")

        pcoa_plots[[dist_method]] <- p
      }

      # Combine plots
      if (length(pcoa_plots) >= 2) {
        combined_pcoa <- plot_grid(
          plotlist = pcoa_plots[seq_len(min(4, length(pcoa_plots)))],
          ncol = 2, align = "hv"
        )

        ggsave("results/final_visualizations/beta_diversity_pcoa_summary.png",
          combined_pcoa,
          width = 14, height = 10, dpi = 300
        )
      }
    }
  }
}

# 3. Rarefaction Curves
create_rarefaction_curves <- function() {
  cat("Creating rarefaction curves...\n")

  # Load count data
  count_files <- list.files("results/taxonomy/merged_tables",
    pattern = "*_abundance_counts.csv", full.names = TRUE
  )

  if (length(count_files) > 0) {
    # Use genus level for rarefaction
    genus_file <- count_files[grepl("genus", count_files)][1]

    if (!is.null(genus_file) && file.exists(genus_file)) {
      count_data <- read.csv(genus_file, row.names = 1, check.names = FALSE)
      count_matrix <- t(count_data) # Transpose for vegan

      # Remove samples with zero counts
      count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

      if (nrow(count_matrix) > 0) {
        # Calculate rarefaction curves
        min_sample_size <- min(rowSums(count_matrix))

        if (min_sample_size > 100) {
          # Use a more reasonable step size
          step_size <- max(100, min_sample_size %/% 50)
          rare_curves <- rarecurve(count_matrix, step = step_size, sample = min_sample_size)

          # Extract data for plotting
          rare_data <- data.frame()

          for (i in seq_len(length(rare_curves))) {
            sample_name <- names(rare_curves)[i]

            # Check if sample_name exists
            if (is.null(sample_name) || sample_name == "") {
              sample_name <- rownames(count_matrix)[i]
            }

            curve_values <- rare_curves[[i]]

            # Check if curve_values has proper names
            if (is.null(names(curve_values)) || length(names(curve_values)) == 0) {
              # Create sequential read counts
              reads_seq <- seq(step_size, min_sample_size, by = step_size)
              if (length(reads_seq) != length(curve_values)) {
                reads_seq <- seq_len(length(curve_values)) * step_size
              }
            } else {
              reads_seq <- as.numeric(names(curve_values))
              # Remove NAs
              valid_indices <- !is.na(reads_seq)
              reads_seq <- reads_seq[valid_indices]
              curve_values <- curve_values[valid_indices]
            }

            # Only proceed if we have valid data
            if (length(reads_seq) > 0 && length(curve_values) > 0 &&
              length(reads_seq) == length(curve_values)) {
              curve_data <- data.frame(
                SampleID = rep(sample_name, length(reads_seq)),
                Reads = reads_seq,
                Taxa = as.numeric(curve_values),
                stringsAsFactors = FALSE
              )

              rare_data <- rbind(rare_data, curve_data)
            }
          }

          # Only proceed if we have data
          if (nrow(rare_data) > 0) {
            # Load metadata for coloring
            metadata_file <- "results/taxonomy/sample_categories.csv"
            if (file.exists(metadata_file)) {
              metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
              rare_data <- merge(rare_data, metadata, by = "SampleID", all.x = TRUE)

              # Set default category for samples without metadata
              rare_data$Category[is.na(rare_data$Category)] <- "Unknown"

              # Create plot
              p_rare <- ggplot(rare_data, aes(
                x = .data$Reads, y = .data$Taxa,
                group = .data$SampleID, color = .data$Category
              )) +
                geom_line(alpha = 0.7) +
                scale_color_brewer(type = "qual", palette = "Set1") +
                labs(
                  title = "Rarefaction Curves (Genus Level)",
                  x = "Number of Reads", y = "Number of Taxa"
                ) +
                theme_minimal() +
                theme(legend.position = "bottom")

              ggsave("results/final_visualizations/rarefaction_curves.png",
                p_rare,
                width = 10, height = 8, dpi = 300
              )

              cat("Rarefaction curves plot saved successfully.\n")
            }
          } else {
            cat("Warning: No valid rarefaction data to plot.\n")
          }
        } else {
          cat("Warning: Sample sizes too small for meaningful rarefaction curves.\n")
        }
      } else {
        cat("Warning: No samples with non-zero counts found.\n")
      }
    }
  }
}

# 4. Statistical Comparison Summary
create_statistical_summary_plot <- function(alpha_results) {
  cat("Creating statistical summary plots...\n")

  # Load statistical test results if available
  stat_files <- list.files("results/diversity/tables",
    pattern = "*_statistical_tests.csv", full.names = TRUE
  )

  if (length(stat_files) > 0) {
    stat_summary <- data.frame()

    for (file in stat_files) {
      level <- gsub("_statistical_tests.csv", "", basename(file))
      stats <- read.csv(file, stringsAsFactors = FALSE)
      stats$Level <- str_to_title(level)
      stat_summary <- rbind(stat_summary, stats)
    }

    if (nrow(stat_summary) > 0) {
      # P-value comparison plot
      p_stats <- ggplot(stat_summary, aes(x = .data$Level, y = -log10(.data$P_value), fill = .data$Level)) +
        geom_col(alpha = 0.7) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        scale_fill_brewer(type = "qual", palette = "Set3") +
        labs(
          title = "Statistical Significance Across Taxonomic Levels",
          subtitle = "Kruskal-Wallis test for Shannon diversity",
          x = "Taxonomic Level", y = "-log10(P-value)"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )

      ggsave("results/final_visualizations/statistical_summary.png",
        p_stats,
        width = 10, height = 6, dpi = 300
      )
    }
  }

  # Create effect size comparison if possible
  create_effect_size_plot(alpha_results)
}

# 5. Effect Size Visualization
create_effect_size_plot <- function(alpha_results) {
  effect_data <- data.frame()

  for (level in names(alpha_results)) {
    alpha_data <- alpha_results[[level]]

    if ("Category" %in% colnames(alpha_data) && length(unique(alpha_data$Category)) == 2) {
      categories <- unique(alpha_data$Category)

      # Calculate effect sizes (Cohen's d) for Shannon diversity
      group1 <- alpha_data$Shannon[alpha_data$Category == categories[1]]
      group2 <- alpha_data$Shannon[alpha_data$Category == categories[2]]

      # Check if we have valid data
      if (length(group1) > 0 && length(group2) > 0 &&
        !all(is.na(group1)) && !all(is.na(group2))) {
        # Cohen's d calculation
        pooled_sd <- sqrt(((length(group1) - 1) * var(group1, na.rm = TRUE) +
          (length(group2) - 1) * var(group2, na.rm = TRUE)) /
          (length(group1) + length(group2) - 2))

        if (pooled_sd > 0) {
          cohens_d <- (mean(group1, na.rm = TRUE) - mean(group2, na.rm = TRUE)) / pooled_sd

          effect_data <- rbind(effect_data, data.frame(
            Level = str_to_title(level),
            Effect_Size = abs(cohens_d),
            Comparison = paste(categories[1], "vs", categories[2]),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  if (nrow(effect_data) > 0) {
    p_effect <- ggplot(effect_data, aes(x = .data$Level, y = .data$Effect_Size, fill = .data$Level)) +
      geom_col(alpha = 0.7) +
      geom_hline(yintercept = c(0.2, 0.5, 0.8), linetype = "dashed", alpha = 0.5) +
      scale_fill_brewer(type = "qual", palette = "Set3") +
      labs(
        title = "Effect Sizes (Cohen's d) for Shannon Diversity",
        subtitle = "Dashed lines: small (0.2), medium (0.5), large (0.8) effects",
        x = "Taxonomic Level", y = "Effect Size (|Cohen's d|)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )

    ggsave("results/final_visualizations/effect_sizes.png",
      p_effect,
      width = 10, height = 6, dpi = 300
    )
  }
}

# Main diversity visualization function
generate_diversity_visualizations <- function() {
  # Load all diversity results
  diversity_data <- load_diversity_results()

  # Generate all plots
  create_alpha_diversity_comparison(diversity_data$alpha)
  create_beta_diversity_summary(diversity_data$beta, diversity_data$alpha)
  create_rarefaction_curves()
  create_statistical_summary_plot(diversity_data$alpha)

  cat("Diversity analysis visualizations completed!\n")
}

# Execute diversity visualizations
generate_diversity_visualizations()


# ==============================================================================
# SECTION 4: DIFFERENTIAL ABUNDANCE PLOTS
# ==============================================================================

cat("Generating Differential Abundance plots...\n")

# Load differential abundance results
load_differential_results <- function() {
  results_dir <- "results/differential_abundance/tables"

  if (!dir.exists(results_dir)) {
    cat("Warning: Differential abundance results not found. Skipping DA plots.\n")
    return(NULL)
  }

  # Get all result files
  deseq_files <- list.files(results_dir, pattern = "deseq2.*\\.csv$", full.names = TRUE)
  aldex_files <- list.files(results_dir, pattern = "aldex2.*\\.csv$", full.names = TRUE)
  lefse_files <- list.files(results_dir, pattern = "lefse.*\\.csv$", full.names = TRUE)

  results <- list(
    deseq2 = lapply(deseq_files, read.csv),
    aldex2 = lapply(aldex_files, read.csv),
    lefse = lapply(lefse_files, read.csv)
  )

  # Add file names for reference
  names(results$deseq2) <- basename(deseq_files)
  names(results$aldex2) <- basename(aldex_files)
  names(results$lefse) <- basename(lefse_files)

  return(results)
}

# Filter to most important comparisons
filter_key_comparisons <- function(results, max_comparisons = 6) {
  if (is.null(results) || length(results$deseq2) == 0) {
    return(results)
  }

  # Calculate significance scores for each comparison
  comparison_scores <- data.frame()

  for (i in seq_along(results$deseq2)) {
    file_name <- names(results$deseq2)[i]
    data <- results$deseq2[[i]]

    if (nrow(data) == 0) next

    # Score based on number of significant results and effect sizes
    sig_count <- sum(data$padj < 0.05, na.rm = TRUE)
    high_sig_count <- sum(data$padj < 0.01, na.rm = TRUE)
    max_effect <- max(abs(data$log2FoldChange), na.rm = TRUE)

    score <- sig_count + (high_sig_count * 2) + (max_effect * 0.5)

    comparison_scores <- rbind(comparison_scores, data.frame(
      file_name = file_name,
      index = i,
      score = score,
      sig_count = sig_count
    ))
  }

  # Select top comparisons
  if (nrow(comparison_scores) > max_comparisons) {
    top_indices <- head(comparison_scores[order(comparison_scores$score, decreasing = TRUE), "index"], max_comparisons)

    # Filter results to top comparisons
    results$deseq2 <- results$deseq2[top_indices]

    cat(sprintf(
      "Filtered to top %d most significant comparisons out of %d total\n",
      length(top_indices), nrow(comparison_scores)
    ))
  }

  return(results)
}

# Combined volcano plot for multiple comparisons
create_combined_volcano <- function(deseq_results, max_plots = 6) {
  if (length(deseq_results) == 0) {
    return(NULL)
  }

  combined_data <- data.frame()

  for (i in seq_along(deseq_results)) {
    if (i > max_plots) break

    data <- deseq_results[[i]]
    if (nrow(data) == 0) next

    # Clean data
    plot_data <- data[!is.na(data$padj) & !is.na(data$log2FoldChange), ]

    # Add comparison identifier
    file_name <- names(deseq_results)[i]
    comparison_name <- str_remove(file_name, "^.*?_deseq2_") %>%
      str_remove("\\.csv$") %>%
      str_trunc(20)

    plot_data$comparison <- comparison_name
    combined_data <- rbind(combined_data, plot_data)
  }

  if (nrow(combined_data) == 0) {
    return(NULL)
  }

  p <- ggplot(combined_data, aes(x = .data$log2FoldChange, y = -log10(.data$padj))) +
    geom_point(aes(color = .data$significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Not Significant" = "grey70",
      "Significant (p < 0.05)" = "#F39C12",
      "Highly Significant (p < 0.01)" = "#E74C3C"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5, color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5, color = "red") +
    facet_wrap(~ .data$comparison, scales = "free", ncol = 3) +
    theme_minimal() +
    labs(
      title = "Volcano Plots - Key Comparisons",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significance"
    ) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 8)
    )

  return(p)
}

# Summary heatmap of significant taxa across comparisons
create_significance_heatmap <- function(deseq_results) {
  if (length(deseq_results) == 0) {
    return(NULL)
  }

  # Collect all significant taxa across comparisons
  sig_matrix <- data.frame()

  for (i in seq_along(deseq_results)) {
    data <- deseq_results[[i]]
    if (nrow(data) == 0) next

    file_name <- names(deseq_results)[i]
    comparison_name <- str_remove(file_name, "^.*?_deseq2_") %>%
      str_remove("\\.csv$") %>%
      str_trunc(15)

    # Get significant taxa
    sig_taxa <- data[data$padj < 0.05 & !is.na(data$padj), ]

    if (nrow(sig_taxa) > 0) {
      sig_taxa$comparison <- comparison_name
      sig_taxa$log2FC_signed <- ifelse(sig_taxa$log2FoldChange > 0, 1, -1)

      sig_matrix <- rbind(sig_matrix, sig_taxa[, c("taxon", "comparison", "log2FC_signed", "padj")])
    }
  }

  if (nrow(sig_matrix) == 0) {
    return(NULL)
  }

  # Create matrix for heatmap (limit to top taxa)
  taxa_counts <- table(sig_matrix$taxon)
  top_taxa <- names(head(sort(taxa_counts, decreasing = TRUE), 20))

  heatmap_data <- sig_matrix[sig_matrix$taxon %in% top_taxa, ] %>%
    select("taxon", "comparison", "log2FC_signed") %>%
    distinct() %>%
    pivot_wider(names_from = "comparison", values_from = "log2FC_signed", values_fill = 0)

  # Convert to matrix for plotting
  taxa_names <- heatmap_data$taxon
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- str_trunc(taxa_names, 40)

  # Convert back to long format for ggplot
  plot_data <- expand.grid(
    Taxon = factor(rownames(heatmap_matrix), levels = rev(rownames(heatmap_matrix))),
    Comparison = factor(colnames(heatmap_matrix))
  ) %>%
    mutate(
      Value = as.vector(heatmap_matrix),
      Significance = case_when(
        .data$Value == 1 ~ "Enriched",
        .data$Value == -1 ~ "Depleted",
        TRUE ~ "Not Significant"
      )
    )

  p <- ggplot(plot_data, aes(x = .data$Comparison, y = .data$Taxon, fill = .data$Significance)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = c(
      "Enriched" = "#E74C3C",
      "Depleted" = "#3498DB",
      "Not Significant" = "grey90"
    )) +
    theme_minimal() +
    labs(
      title = "Differential Abundance Heatmap - Top 20 Taxa",
      x = "Comparison",
      y = "Taxon",
      fill = "Direction"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )

  return(p)
}

# Enhanced Significant Taxa Bar Plot with levels
create_enhanced_significant_taxa_plot <- function(results_data, method = "deseq2") {
  if (is.null(results_data) || length(results_data) == 0) {
    return(NULL)
  }

  # Count significant taxa by level and comparison
  sig_counts <- data.frame()

  for (file_name in names(results_data)) {
    data <- results_data[[file_name]]
    if (nrow(data) == 0) next

    # Extract level and comparison from filename
    parts <- str_split(str_remove(file_name, "\\.csv$"), "_")[[1]]
    level <- parts[2] # Assuming format like "genus_deseq2_comparison.csv"
    comparison <- paste(parts[3:length(parts)], collapse = "_")
    comparison <- str_trunc(comparison, 20)

    if (method == "deseq2") {
      sig_05 <- sum(data$padj < 0.05, na.rm = TRUE)
      sig_01 <- sum(data$padj < 0.01, na.rm = TRUE)
    } else if (method == "aldex2") {
      sig_05 <- sum(data$p_value_bh < 0.05, na.rm = TRUE)
      sig_01 <- sum(data$p_value_bh < 0.01, na.rm = TRUE)
    }

    sig_counts <- rbind(sig_counts, data.frame(
      Level = level,
      Comparison = comparison,
      `p < 0.05` = sig_05,
      `p < 0.01` = sig_01,
      check.names = FALSE
    ))
  }

  if (nrow(sig_counts) == 0) {
    return(NULL)
  }

  # Reshape for plotting
  plot_data <- sig_counts %>%
    pivot_longer(
      cols = c("p < 0.05", "p < 0.01"),
      names_to = "Significance",
      values_to = "Count"
    )

  p <- ggplot(plot_data, aes(x = .data$Level, y = .data$Count, fill = .data$Significance)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("p < 0.05" = "#F39C12", "p < 0.01" = "#E74C3C")) +
    theme_minimal() +
    labs(
      title = paste("Significant Taxa Count Across All Comparisons -", str_to_title(method)),
      x = "Taxonomic Level",
      y = "Total Number of Significant Taxa",
      fill = "Significance Level"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}

# Main plotting function - streamlined
plot_differential_abundance <- function() {
  results <- load_differential_results()
  if (is.null(results)) {
    return()
  }

  # Filter to key comparisons to reduce plot overload
  results <- filter_key_comparisons(results, max_comparisons = 6)

  plot_list <- list()

  # 1. Combined volcano plot for key comparisons
  if (length(results$deseq2) > 0) {
    combined_volcano <- create_combined_volcano(results$deseq2)
    if (!is.null(combined_volcano)) {
      plot_list[["combined_volcano"]] <- combined_volcano
    }
  }

  # 2. Enhanced significant taxa summary plots
  if (length(results$deseq2) > 0) {
    sig_bar_deseq2 <- create_enhanced_significant_taxa_plot(results$deseq2, "deseq2")
    if (!is.null(sig_bar_deseq2)) {
      plot_list[["sig_summary_deseq2"]] <- sig_bar_deseq2
    }
  }

  if (length(results$aldex2) > 0) {
    sig_bar_aldex2 <- create_enhanced_significant_taxa_plot(results$aldex2, "aldex2")
    if (!is.null(sig_bar_aldex2)) {
      plot_list[["sig_summary_aldex2"]] <- sig_bar_aldex2
    }
  }

  # 3. Method comparison plot (if both methods available)
  if (length(results$deseq2) > 0 && length(results$aldex2) > 0) {
    # Create one representative comparison plot
    deseq_file <- names(results$deseq2)[1]
    corresponding_aldex <- NULL

    for (j in seq_along(results$aldex2)) {
      aldex_file <- names(results$aldex2)[j]
      if (str_detect(deseq_file, str_extract(aldex_file, "(?<=aldex2_).*(?=\\.csv)"))) {
        corresponding_aldex <- results$aldex2[[j]]
        break
      }
    }

    if (!is.null(corresponding_aldex)) {
      # Merge data for method comparison
      merged_data <- merge(
        results$deseq2[[1]][, c("taxon", "log2FoldChange", "padj")],
        corresponding_aldex[, c("taxon", "effect_size", "p_value_bh")],
        by = "taxon",
        suffix = c("_deseq2", "_aldex2")
      )

      if (nrow(merged_data) > 0) {
        merged_data$category <- "Not Significant"
        merged_data$category[merged_data$padj < 0.05 & merged_data$p_value_bh < 0.05] <- "Both Significant"
        merged_data$category[merged_data$padj < 0.05 & merged_data$p_value_bh >= 0.05] <- "DESeq2 Only"
        merged_data$category[merged_data$padj >= 0.05 & merged_data$p_value_bh < 0.05] <- "ALDEx2 Only"

        method_comparison <- ggplot(merged_data, aes(x = .data$log2FoldChange, y = .data$effect_size)) +
          geom_point(aes(color = .data$category), alpha = 0.7, size = 2) +
          scale_color_manual(values = c(
            "Not Significant" = "grey70",
            "Both Significant" = "#E74C3C",
            "DESeq2 Only" = "#3498DB",
            "ALDEx2 Only" = "#2ECC71"
          )) +
          geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
          theme_minimal() +
          labs(
            title = "DESeq2 vs ALDEx2 Method Comparison",
            x = "DESeq2 Log2 Fold Change",
            y = "ALDEx2 Effect Size",
            color = "Significance"
          ) +
          theme(legend.position = "bottom")

        plot_list[["method_comparison"]] <- method_comparison
      }
    }
  }

  # Save plots with descriptive names
  for (plot_name in names(plot_list)) {
    ggsave(
      filename = file.path("results/final_visualizations", paste0("DA_", plot_name, ".png")),
      plot = plot_list[[plot_name]],
      width = 14, height = 10, dpi = 300
    )
  }

  cat(sprintf(
    "Generated %d summary differential abundance plots instead of %d individual plots\n",
    length(plot_list), length(results$deseq2) * 2
  ))
  cat("Differential abundance plots saved to results/final_visualizations/\n")

  return(plot_list)
}

# Run the optimized plotting function
plots <- plot_differential_abundance()


# ==============================================================================
# SECTION 5: MACHINE LEARNING & BIOMARKER PLOTS
# ==============================================================================
cat("Generating ML and Biomarker Discovery plots...\n")

# Function to read ML results
read_ml_results <- function(results_dir = "results") {
  ml_dir <- file.path(results_dir, "ml_biomarker_discovery")

  # Read biomarker rankings for each taxonomic level
  biomarker_files <- list.files(file.path(ml_dir, "biomarker_panels"),
    pattern = "_biomarker_ranking.csv",
    full.names = TRUE
  )

  biomarkers_list <- list()
  for (file in biomarker_files) {
    level <- gsub(".*/(\\w+)_biomarker_ranking.csv", "\\1", file)
    biomarkers_list[[level]] <- read.csv(file, row.names = 1) %>%
      mutate(taxonomic_level = level) %>%
      rownames_to_column("feature")
  }

  # Read combined biomarker results
  combined_file <- file.path(ml_dir, "biomarker_panels", "combined_top_biomarkers.csv")
  combined_biomarkers <- NULL
  if (file.exists(combined_file)) {
    combined_biomarkers <- read.csv(combined_file)
  }

  # Read classification reports (parse text file)
  report_files <- list.files(file.path(ml_dir, "classification_reports"),
    pattern = "_classification_report.txt",
    full.names = TRUE
  )

  model_performance <- data.frame()
  for (file in report_files) {
    level <- gsub(".*/(\\w+)_classification_report.txt", "\\1", file)
    lines <- readLines(file)

    # Extract accuracy values
    accuracy_lines <- lines[grepl("CV Accuracy:", lines)]
    if (length(accuracy_lines) > 0) {
      models <- c("Random Forest", "Gradient Boosting", "SVM", "Logistic Regression")
      for (i in seq_along(accuracy_lines)) {
        if (i <= length(models)) {
          acc_value <- as.numeric(gsub(".*CV Accuracy: ([0-9.]+).*", "\\1", accuracy_lines[i]))
          model_performance <- rbind(
            model_performance,
            data.frame(
              level = level,
              model = models[i],
              accuracy = acc_value
            )
          )
        }
      }
    }
  }

  return(list(
    biomarkers = biomarkers_list,
    combined = combined_biomarkers,
    performance = model_performance
  ))
}

# Read ML results
ml_results <- read_ml_results()

# 1. Feature Importance Heatmap
if (length(ml_results$biomarkers) > 0) {
  # Combine top biomarkers from each level
  top_biomarkers <- do.call(rbind, lapply(ml_results$biomarkers, function(x) {
    x %>%
      head(15) %>%
      select(feature, mean_importance, taxonomic_level)
  }))

  # Create feature importance heatmap
  importance_matrix <- top_biomarkers %>%
    select(feature, taxonomic_level, mean_importance) %>%
    pivot_wider(names_from = taxonomic_level, values_from = mean_importance, values_fill = 0) %>%
    column_to_rownames("feature") %>%
    as.matrix()

  # Generate and save heatmap directly
  png("results/final_visualizations/feature_importance_heatmap.png",
    width = 12, height = 10, units = "in", res = 300
  )

  pheatmap(importance_matrix,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color = colorRampPalette(c("white", "red", "darkred"))(50),
    main = "Feature Importance Across Taxonomic Levels",
    fontsize = 8
  )

  dev.off()
}

# 2. Model Performance Comparison
if (nrow(ml_results$performance) > 0) {
  p_performance <- ggplot(ml_results$performance, aes(x = model, y = accuracy, fill = level)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = round(accuracy, 3)),
      position = position_dodge(width = 0.9),
      vjust = -0.3, size = 3
    ) +
    labs(
      title = "Model Performance Across Taxonomic Levels",
      x = "Model Type",
      y = "Cross-Validation Accuracy",
      fill = "Taxonomic Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    scale_y_continuous(limits = c(0, max(ml_results$performance$accuracy) * 1.1))

  ggsave("results/final_visualizations/model_performance_comparison.png",
    plot = p_performance, width = 12, height = 8, dpi = 300
  )
}

# 3. Top Biomarkers Bar Plot
if (!is.null(ml_results$combined)) {
  top_n <- min(20, nrow(ml_results$combined))

  p_biomarkers <- ml_results$combined %>%
    head(top_n) %>%
    mutate(feature = factor(feature, levels = rev(feature))) %>%
    ggplot(aes(x = mean_importance, y = feature, fill = level)) +
    geom_col(alpha = 0.8) +
    labs(
      title = "Top Biomarkers Across All Taxonomic Levels",
      x = "Mean Feature Importance",
      y = "Biomarker",
      fill = "Taxonomic Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set3")

  ggsave("results/final_visualizations/top_biomarkers.png",
    plot = p_biomarkers, width = 14, height = 10, dpi = 300
  )
}

# 4. Biomarker Distribution by Taxonomic Level
if (length(ml_results$biomarkers) > 0) {
  biomarker_counts <- sapply(ml_results$biomarkers, nrow)

  p_distribution <- data.frame(
    level = names(biomarker_counts),
    count = biomarker_counts
  ) %>%
    ggplot(aes(x = reorder(level, count), y = count, fill = level)) +
    geom_col(alpha = 0.8, show.legend = FALSE) +
    geom_text(aes(label = count), hjust = -0.1) +
    coord_flip() +
    labs(
      title = "Number of Biomarkers by Taxonomic Level",
      x = "Taxonomic Level",
      y = "Number of Biomarkers"
    ) +
    theme_minimal() +
    scale_fill_viridis_d()

  ggsave("results/final_visualizations/biomarker_distribution.png",
    plot = p_distribution, width = 10, height = 6, dpi = 300
  )
}

cat("ML visualization plots completed!\n")


# ==============================================================================
# SECTION 6: INTEGRATED SUMMARY DASHBOARD
# ==============================================================================

cat("Creating integrated summary dashboard...\n")

# Function to create summary dashboard
create_summary_dashboard <- function() {
  # Read key data files
  qc_data <- tryCatch(read.csv("results/qc/trimming_summary.csv"), error = function(e) NULL)
  alpha_data <- tryCatch(read.csv("results/diversity/tables/genus_alpha_diversity.csv"), error = function(e) NULL)

  # Create summary plots list
  plots <- list()

  # 1. QC Summary mini plot
  if (!is.null(qc_data)) {
    p1 <- ggplot(qc_data, aes(
      x = reorder(sample, .data$reads_retained_percent),
      y = .data$reads_retained_percent
    )) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      labs(title = "Read Retention (%)", x = "", y = "%") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    plots$qc <- p1
  }

  # 2. Diversity summary
  if (!is.null(alpha_data) && "Category" %in% colnames(alpha_data)) {
    p2 <- ggplot(alpha_data, aes(x = .data$Category, y = .data$Shannon, fill = .data$Category)) +
      geom_boxplot(alpha = 0.7) +
      labs(title = "Shannon Diversity", x = "", y = "Index") +
      theme_minimal() +
      theme(legend.position = "none")
    plots$diversity <- p2
  }

  # 3. Taxa richness
  if (!is.null(alpha_data)) {
    p3 <- ggplot(alpha_data, aes(x = .data$Category, y = .data$Observed_Taxa, fill = .data$Category)) +
      geom_boxplot(alpha = 0.7) +
      labs(title = "Taxa Richness", x = "", y = "Count") +
      theme_minimal() +
      theme(legend.position = "none")
    plots$richness <- p3
  }

  # 4. Sample size summary
  if (!is.null(qc_data)) {
    p4 <- ggplot(qc_data, aes(x = 1, y = .data$total_reads_after / 1000)) +
      geom_boxplot(fill = "lightgreen", alpha = 0.7) +
      labs(title = "Final Read Counts", x = "", y = "Reads (K)") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    plots$reads <- p4
  }

  # Combine plots if available
  if (length(plots) >= 2) {
    combined_dashboard <- plot_grid(
      plotlist = plots,
      labels = LETTERS[seq_along(plots)],
      ncol = 2
    )

    ggsave("results/final_visualizations/summary_dashboard.png",
      combined_dashboard,
      width = 12, height = 8, dpi = 300
    )
  }
}

# Multi-panel publication figure
create_publication_figure <- function() {
  # Load key plots if they exist
  plot_files <- c(
    "results/final_visualizations/phylum_abundance_barplot.png",
    "results/final_visualizations/diversity_overview.png",
    "results/final_visualizations/beta_diversity_pcoa.png"
  )

  existing_plots <- plot_files[file.exists(plot_files)]

  if (length(existing_plots) > 0) {
    cat("Key visualization files found. Publication figure components available.\n")
  }
}

# Execute dashboard creation
create_summary_dashboard()
create_publication_figure()


# ==============================================================================
# SAVE SUMMARY STATISTICS
# ==============================================================================

cat("Compiling summary statistics...\n")

# Initialize summary list
summary_stats <- list()

# QC Statistics
if (file.exists("results/qc/trimming_summary.csv")) {
  qc_data <- read.csv("results/qc/trimming_summary.csv")
  summary_stats$qc <- list(
    total_samples = nrow(qc_data),
    mean_retention = round(mean(qc_data$reads_retained_percent), 2),
    total_input_reads = sum(qc_data$total_reads_before),
    total_output_reads = sum(qc_data$total_reads_after)
  )
}

# Diversity Statistics
alpha_files <- list.files("results/diversity/tables", pattern = "*alpha_diversity.csv", full.names = TRUE)
if (length(alpha_files) > 0) {
  alpha_data <- read.csv(alpha_files[1])
  if ("Shannon" %in% colnames(alpha_data)) {
    summary_stats$diversity <- list(
      mean_shannon = round(mean(alpha_data$Shannon, na.rm = TRUE), 3),
      mean_simpson = round(mean(alpha_data$Simpson, na.rm = TRUE), 3),
      mean_richness = round(mean(alpha_data$Observed_Taxa, na.rm = TRUE), 1)
    )
  }
}

# Taxonomic Statistics
tax_files <- list.files("results/taxonomy/merged_tables", pattern = "*relative_abundance.csv", full.names = TRUE)
if (length(tax_files) > 0) {
  # Get sample info from first file
  tax_data <- read.csv(tax_files[1], row.names = 1)
  summary_stats$taxonomy <- list(
    total_taxa_detected = nrow(tax_data),
    samples_analyzed = ncol(tax_data)
  )
}

# ML Results
if (file.exists("results/ml_biomarker_discovery/biomarker_panels/combined_top_biomarkers.csv")) {
  biomarkers <- read.csv("results/ml_biomarker_discovery/biomarker_panels/combined_top_biomarkers.csv")
  summary_stats$ml <- list(
    total_biomarkers = nrow(biomarkers),
    top_biomarker = biomarkers$feature[1],
    max_importance = round(max(biomarkers$mean_importance), 3)
  )
}

# Save summary statistics
summary_df <- data.frame(
  Analysis = character(),
  Metric = character(),
  Value = character(),
  stringsAsFactors = FALSE
)

for (analysis in names(summary_stats)) {
  for (metric in names(summary_stats[[analysis]])) {
    summary_df <- rbind(summary_df, data.frame(
      Analysis = toupper(analysis),
      Metric = gsub("_", " ", str_to_title(metric)),
      Value = as.character(summary_stats[[analysis]][[metric]])
    ))
  }
}

# Write summary table
write.csv(summary_df, "results/final_visualizations/analysis_summary_statistics.csv", row.names = FALSE)

# Create summary table plot
if (nrow(summary_df) > 0) {
  p_summary <- ggplot() +
    theme_void() +
    labs(title = "Analysis Summary Statistics") +
    annotation_custom(
      tableGrob(summary_df, theme = ttheme_default(base_size = 10))
    )

  ggsave("results/final_visualizations/summary_statistics_table.png",
    p_summary,
    width = 10, height = 8, dpi = 300
  )
}

cat("Summary statistics saved to: results/final_visualizations/analysis_summary_statistics.csv\n")
cat("Dashboard and summary visualizations completed!\n")
