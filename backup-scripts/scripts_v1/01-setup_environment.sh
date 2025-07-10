#!/usr/bin/env bash

# # 16S rRNA Analysis Environment Setup for macOS ARM (M1/M2)
# # This script installs all necessary tools that work natively on ARM processors

# echo "Setting up 16S rRNA analysis environment for macOS ARM..."

# # # Update mamba
# mamba update -n base -c defaults mamba

# # Create conda environment for 16S analysis
# echo "Creating conda environment: microbiome_analysis"
# mamba create -n microbiome_analysis python=3.9 -y

# # Activate environment
# mamba activate microbiome_analysis

# # Install quality control tools
# echo "Installing quality control tools..."
# mamba install -c bioconda fastqc -y
# mamba install -c bioconda fastp -y
# mamba install -c bioconda multiqc -y

# # Install taxonomic profiling tools (ARM compatible)
# echo "Installing taxonomic profiling tools..."
# mamba install -c bioconda kraken2 -y
# mamba install -c bioconda bracken -y
# mamba install -c bioconda metaphlan -y

# # Install sequence processing tools
# echo "Installing sequence processing tools..."
# mamba install -c bioconda seqtk -y
# mamba install -c bioconda bbmap -y

# # Install R and necessary packages
# echo "Installing R and packages..."
# mamba install -c conda-forge r-base=4.3 -y
# mamba install -c conda-forge r-essentials -y
# mamba install -c bioconda bioconductor-phyloseq -y
# mamba install -c conda-forge r-vegan -y
# mamba install -c conda-forge r-ggplot2 -y
# mamba install -c conda-forge r-dplyr -y
# mamba install -c conda-forge r-tidyr -y
# mamba install -c conda-forge r-readr -y
# mamba install -c bioconda bioconductor-deseq2 -y
# mamba install -c conda-forge r-randomforest -y
# mamba install -c conda-forge r-caret -y

# # Install Python packages for ML and analysis
# echo "Installing Python packages..."
# pip install pandas numpy scipy matplotlib seaborn
# pip install scikit-learn xgboost lightgbm
# pip install biopython plotly
# pip install jupyter notebook

# # Install additional analysis tools
# echo "Installing additional tools..."
# mamba install -c bioconda vsearch -y
# mamba install -c bioconda mafft -y
# mamba install -c bioconda fasttree -y

# Create directory structure
echo "Creating directory structure..."
mkdir -p results/{qc,taxonomy,diversity,differential,ml,visualization}
# mkdir -p databases
mkdir -p logs

# echo "Environment setup complete!"
# echo "To activate the environment, run: conda activate microbiome_analysis"

echo "Setup completed successfully!"
