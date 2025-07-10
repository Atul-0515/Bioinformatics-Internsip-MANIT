# 16S rRNA Diversity Analysis Script Breakdown

## 1. Script Purpose
This script analyzes **microbial diversity** in your gut microbiome samples. It measures two main types of diversity:
- **Alpha diversity**: How many different microbes are in each individual sample and how evenly distributed they are
- **Beta diversity**: How similar or different the microbial communities are between different samples/groups

Think of it like analyzing the biodiversity in different forests - alpha diversity tells you how many species are in each forest, while beta diversity tells you how similar the forests are to each other.

## 2. Key Functions

### Data Loading and Preparation
- Loads the taxonomic abundance data from the previous analysis step
- Reads sample information (which group each sample belongs to: normal, obesity, diabetes, or both)
- Organizes data at different taxonomic levels (species, genus, family, etc.)

### Alpha Diversity Calculations
- **Observed Taxa**: Simply counts how many different microbes are present
- **Shannon Index**: Measures both richness (number of species) and evenness (how equally distributed they are)
- **Simpson Index**: Focuses more on the most abundant species
- **Chao1**: Estimates the total number of species, including rare ones that might have been missed

### Beta Diversity Analysis
- **Distance Calculations**: Measures how different samples are from each other using three methods:
  - Bray-Curtis: Good for abundance data
  - Jaccard: Focuses on presence/absence
  - Euclidean: Standard mathematical distance
- **PCoA (Principal Coordinates Analysis)**: Creates 2D plots to visualize how samples cluster together
- **PERMANOVA**: Statistical test to see if groups are significantly different

## 3. Input/Output

### Input Required:
- Taxonomic abundance tables (from the previous taxonomic profiling step)
- Sample metadata file showing which group each sample belongs to
- Must have completed the taxonomic profiling analysis first

### Output Produced:
- **Tables**: CSV files with diversity metrics for each sample
- **Plots**: Box plots, violin plots, PCoA plots, and heatmaps
- **Statistical Results**: PERMANOVA test results and summary statistics
- All saved in organized folders under `results/diversity/`

## 4. Methods Used

### Statistical Methods:
- **Kruskal-Wallis Test**: Non-parametric test to compare diversity between groups (doesn't assume normal distribution)
- **PERMANOVA**: Tests if microbial community composition differs significantly between groups
- **Principal Coordinates Analysis (PCoA)**: Dimensionality reduction to visualize sample relationships

### R Packages Used:
- **vegan**: Main package for ecological diversity analysis
- **ggplot2**: Creates publication-quality plots
- **dplyr/tidyr**: Data manipulation and organization
- **pheatmap**: Creates heatmaps for distance matrices

## 5. Results Generated

### Visual Outputs:
- **Alpha Diversity Plots**: Box plots and violin plots showing diversity differences between your 4 groups
- **PCoA Plots**: Scatter plots showing how samples cluster by group
- **Distance Heatmaps**: Color-coded matrices showing similarity between all sample pairs

### Statistical Outputs:
- **Summary Tables**: Mean and standard deviation of diversity metrics for each group
- **Statistical Test Results**: P-values showing if groups are significantly different
- **Distance Matrices**: Numerical values showing how different each sample is from every other sample