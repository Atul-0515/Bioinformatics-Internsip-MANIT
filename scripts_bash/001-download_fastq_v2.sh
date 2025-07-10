#!/usr/bin/env bash

# Simple script to download FASTQ files from SRR/ERR numbers using ENA
# Usage: ./download_fastq_ena.sh

# Create output directory
mkdir -p fastq_files

# Function to get ENA download URLs and download files
download_fastq() {
    local txt_file=$1
    local condition=$(basename "$txt_file" .txt)
    
    # Count total accessions in file
    local total=$(grep -v '^$' "$txt_file" | grep -v '^#' | wc -l)
    local completed=0
    
    echo "Processing $txt_file for condition: $condition"
    echo "Total accessions: $total"
    mkdir -p "fastq_files/$condition"
    
    while IFS= read -r accession; do
        # Skip empty lines and comments
        [[ -z "$accession" || "$accession" =~ ^#.*$ ]] && continue
        
        # Clean accession
        accession=$(echo "$accession" | tr -d '[:space:]')
        
        ((completed++))
        local remaining=$((total - completed))
        
        echo "[$completed/$total] Downloading $accession... ($remaining remaining)"
        
        # Get ENA URLs
        local urls_file="temp_${accession}.txt"
        wget -q -O "$urls_file" "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"
        
        # Extract and download each FASTQ file
        tail -n +2 "$urls_file" | cut -f2 | tr ';' '\n' | grep -v '^$' | while read url; do
            if [ -n "$url" ]; then
                echo "  -> Downloading $(basename $url)..."
                wget -c -P "fastq_files/$condition" "ftp://$url"
            fi
        done
        
        rm -f "$urls_file"
        echo "  ✓ Completed $accession"
        
    done < "$txt_file"
    
    echo "Finished $condition: $completed/$total accessions completed"
    echo ""
}

# Process all txt files in runs directory
for txt_file in runs/*.txt; do
    if [ -f "$txt_file" ]; then
        download_fastq "$txt_file"
    fi
done

echo "Download complete! Files organized in fastq_files/ by condition."