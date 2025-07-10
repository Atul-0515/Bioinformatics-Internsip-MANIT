#!/usr/bin/env bash

# Quality Control Script for 16S rRNA FASTQ files
# Performs initial QC, trimming, and post-trim QC

# Activate environment
# mamba activate microbiome_analysis


# Set variables
FASTQ_DIR="fastq_files"
QC_DIR="results/qc"
TRIM_DIR="results/qc/trimmed"
THREADS=8

echo "Starting Quality Control Analysis..."

# Create output directories
mkdir -p ${QC_DIR}/raw_fastqc
mkdir -p ${QC_DIR}/trimmed_fastqc
mkdir -p ${TRIM_DIR}
mkdir -p logs

# Function to process each sample directory
process_sample() {
    local sample_dir=$1
    local sample_name=$(basename "$sample_dir")
    
    echo "Processing sample: $sample_name"
    
    # Create sample-specific directories
    mkdir -p ${QC_DIR}/raw_fastqc/${sample_name}
    mkdir -p ${QC_DIR}/trimmed_fastqc/${sample_name}
    mkdir -p ${TRIM_DIR}/${sample_name}
    
    # Run FastQC on raw files
    echo "Running FastQC on raw files for $sample_name..."
    fastqc ${sample_dir}/*.fastq.gz \
        -o ${QC_DIR}/raw_fastqc/${sample_name} \
        -t ${THREADS} \
        2> logs/fastqc_raw_${sample_name}.log
    
    # Process each pair of files
    for r1_file in ${sample_dir}/*_1.fastq.gz; do
        if [[ -f "$r1_file" ]]; then
            # Get corresponding R2 file
            r2_file=${r1_file/_1.fastq.gz/_2.fastq.gz}
            
            # Extract base name for output
            base_name=$(basename "$r1_file" _1.fastq.gz)
            
            echo "Trimming: $base_name"
            
            # Run fastp for quality trimming (gentle parameters for 16S)
            fastp \
                -i "$r1_file" \
                -I "$r2_file" \
                -o "${TRIM_DIR}/${sample_name}/${base_name}_1_trimmed.fastq.gz" \
                -O "${TRIM_DIR}/${sample_name}/${base_name}_2_trimmed.fastq.gz" \
                --thread ${THREADS} \
                --qualified_quality_phred 20 \
                --unqualified_percent_limit 20 \
                --length_required 100 \
                --detect_adapter_for_pe \
                --correction \
                --json "${QC_DIR}/${sample_name}_${base_name}_fastp.json" \
                --html "${QC_DIR}/${sample_name}_${base_name}_fastp.html" \
                2> logs/fastp_${sample_name}_${base_name}.log
        fi
    done
    
    # Run FastQC on trimmed files
    echo "Running FastQC on trimmed files for $sample_name..."
    fastqc ${TRIM_DIR}/${sample_name}/*.fastq.gz \
        -o ${QC_DIR}/trimmed_fastqc/${sample_name} \
        -t ${THREADS} \
        2> logs/fastqc_trimmed_${sample_name}.log
}

# Process all sample directories
for sample_dir in ${FASTQ_DIR}/*/; do
    if [[ -d "$sample_dir" ]]; then
        process_sample "$sample_dir"
    fi
done

# Generate MultiQC report
echo "Generating MultiQC report..."
multiqc ${QC_DIR} -o ${QC_DIR} -n multiqc_report.html

# Generate summary statistics
echo "Generating summary statistics..."
python3 << 'EOF'
import os
import json
import pandas as pd
from pathlib import Path

# Collect fastp statistics
fastp_stats = []
qc_dir = Path("results/qc")

for json_file in qc_dir.glob("*_fastp.json"):
    with open(json_file) as f:
        data = json.load(f)
    
    sample_info = {
        'sample': json_file.stem.replace('_fastp', ''),
        'total_reads_before': data['summary']['before_filtering']['total_reads'],
        'total_reads_after': data['summary']['after_filtering']['total_reads'],
        'q30_rate_before': data['summary']['before_filtering']['q30_rate'],
        'q30_rate_after': data['summary']['after_filtering']['q30_rate'],
        'adapter_trimmed_reads': data['adapter_cutting']['adapter_trimmed_reads'],
        'reads_passed_filter': data['filtering_result']['passed_filter_reads']
    }
    fastp_stats.append(sample_info)

# Create summary DataFrame
if fastp_stats:
    df = pd.DataFrame(fastp_stats)
    df['reads_retained_percent'] = (df['total_reads_after'] / df['total_reads_before'] * 100).round(2)
    df.to_csv('results/qc/trimming_summary.csv', index=False)
    print("Quality control summary saved to results/qc/trimming_summary.csv")
    print("\nSample processing summary:")
    print(df[['sample', 'total_reads_before', 'total_reads_after', 'reads_retained_percent']])
else:
    print("No fastp JSON files found for summary generation")
EOF

echo "Quality Control Analysis Complete!"
echo "Results saved in: ${QC_DIR}"
echo "Trimmed files saved in: ${TRIM_DIR}"
echo "View the MultiQC report: ${QC_DIR}/multiqc_report.html"