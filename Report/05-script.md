# Differential Abundance Analysis Script Summary

## 1. Script Purpose
This script identifies which bacterial species/taxa are significantly different between your four sample groups (normal, obesity, diabetes, and combined T2D+obesity). Essentially, it answers the question: "Which microbes are more or less abundant in diseased patients compared to healthy controls?"

## 2. Key Functions

### Data Loading and Preparation
- Loads taxonomic abundance data from previous analysis steps
- Filters out rare bacteria (keeps only those present in at least 10% of samples)
- Organizes data by taxonomic levels (species, genus, family, etc.)

### Statistical Analysis Methods
The script uses three different statistical approaches to ensure robust results:

**DESeq2 Analysis:**
- Originally designed for gene expression data, adapted for microbiome
- Good for detecting fold-changes between groups
- Handles the compositional nature of microbiome data

**ALDEx2 Analysis:**
- Specifically designed for compositional data like microbiome
- Uses centered log-ratio transformation
- More conservative approach, reduces false positives

**LEfSe-style Analysis:**
- Uses Kruskal-Wallis test (non-parametric)
- Identifies biomarker taxa that characterize each group
- Calculates effect sizes to measure biological significance

### Visualization Creation
- Volcano plots showing statistical significance vs fold-change
- Heatmaps displaying abundance patterns of significant taxa
- Effect size plots for ALDEx2 results

## 3. Input/Output

### Input Requirements:
- Taxonomic abundance tables (counts and relative abundances) from previous steps
- Sample metadata with group categories
- Results from the taxonomic profiling script (03-taxonomic_profiling.py)

### Output Generated:
- **Tables:** CSV files with statistical results for each comparison and method
- **Plots:** Volcano plots, effect plots, and heatmaps
- **Directories:** Organized results in `results/differential_abundance/`

## 4. Methods Used

### Statistical Tests:
- **DESeq2:** Negative binomial distribution modeling with size factor normalization
- **Kruskal-Wallis test:** Non-parametric test for multiple group comparisons
- **Welch's t-test:** Used within ALDEx2 for pairwise comparisons
- **Multiple testing correction:** Benjamini-Hochberg false discovery rate (FDR)

### Data Processing:
- Pseudocount addition (+1) for DESeq2 to handle zero values
- Log-transformation for visualization
- Centered log-ratio transformation in ALDEx2
- Filtering based on prevalence and abundance thresholds

### Effect Size Calculations:
- Log2 fold-change (DESeq2)
- Effect size based on between/within group variance (ALDEx2)
- LDA-style effect size (LEfSe approach)

## 5. Results Generated

### Statistical Tables:
- Lists of significantly different taxa between each group comparison
- P-values, adjusted p-values, fold-changes, and effect sizes
- Direction of change (enriched in which group)

### Visualizations:
- **Volcano Plots:** Show statistical significance vs biological effect size
- **Heatmaps:** Display abundance patterns of significant taxa across samples
- **Effect Plots:** ALDEx2-specific visualization of effect sizes

### Key Comparisons:
All pairwise comparisons between your four groups:
- Normal vs Obesity
- Normal vs Diabetes  
- Normal vs T2D+Obesity
- Obesity vs Diabetes
- Obesity vs T2D+Obesity
- Diabetes vs T2D+Obesity
