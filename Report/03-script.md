# Taxonomic Profiling Script Analysis

## 1. Script Purpose
This script performs **taxonomic identification and classification** of bacterial species in your 16S rRNA microbiome samples. Essentially, it takes your cleaned sequencing data and identifies which bacteria are present in each sample and how abundant they are. Think of it as creating a detailed census of the bacterial community in each gut microbiome sample.

## 2. Key Functions

### Primary Analytical Steps:
- **Database Setup**: Checks for and extracts the reference database needed to identify bacteria
- **Taxonomic Classification**: Uses Kraken2 to identify bacterial species in each sample by comparing sequences to known bacterial genomes
- **Abundance Estimation**: Uses Bracken to calculate how many reads belong to each bacterial species (correcting for biases)
- **Data Consolidation**: Combines results from all samples into organized tables
- **Metadata Organization**: Creates a mapping file linking each sample to its group (normal, obesity, diabetes, or both)

### Core Processing Functions:
- `process_individual_kraken2()`: Identifies bacteria in each sample
- `run_bracken()`: Calculates bacterial abundances at different taxonomic levels
- `create_abundance_tables()`: Merges data from all samples into comparative tables
- `generate_taxonomy_summary()`: Creates summary statistics

## 3. Input/Output

### Input Required:
- **Trimmed FASTQ files**: Clean sequencing data from the quality control step
- **Kraken2 database**: Reference database containing known bacterial genomes
- **Sample organization**: Samples organized by group (normal, obesity, diabetes, both)

### Output Generated:
- **Abundance tables**: CSV files showing bacterial counts and percentages for each sample
  - Raw counts (actual number of sequences)
  - Relative abundances (percentages)
- **Taxonomic reports**: Detailed classification results for each sample
- **Sample metadata**: File mapping each sample to its clinical group
- **Summary statistics**: Overall processing results and bacterial diversity counts

## 4. Methods Used

### Bioinformatics Tools:
- **Kraken2**: Fast taxonomic classifier that assigns bacterial identities to DNA sequences
- **Bracken**: Abundance estimation tool that corrects for genome size biases

### Classification Approach:
- **Paired-end processing**: Analyzes both forward and reverse reads together for better accuracy
- **Multi-level taxonomy**: Classifies bacteria at 6 different levels:
  - Species (most specific)
  - Genus
  - Family
  - Order
  - Class
  - Phylum (most general)

### Quality Controls:
- File existence and size validation
- Processing success tracking
- Error logging for troubleshooting

## 5. Results Generated

### Main Output Files:
1. **Abundance Tables** (for each taxonomic level):
   - `species_abundance_counts.csv`: Raw bacterial counts
   - `species_relative_abundance.csv`: Bacterial percentages
   - Similar files for genus, family, order, class, and phylum levels

2. **Sample Metadata**:
   - `sample_categories.csv`: Links each sample to its clinical group

3. **Processing Reports**:
   - Individual Kraken2 and Bracken reports for each sample
   - Summary statistics (total samples, taxa detected, processing success rates)

### Data Structure:
- Rows represent different bacterial taxa
- Columns represent individual samples
- Values show either raw counts or relative percentages
