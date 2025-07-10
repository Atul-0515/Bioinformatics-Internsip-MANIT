#!/usr/bin/env python3

"""
Taxonomic Profiling Script for 16S rRNA data
Uses Kraken2 and Bracken for taxonomy assignment and abundance estimation

FIXED: Now correctly handles multiple trimmed files per sample directory
       by treating EACH INDIVIDUAL TRIMMED FASTQ PAIR as a separate biological sample,
       as clarified by user.
       Generates a sample_categories.csv for downstream R script.
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
import json
import tarfile


def run_command(cmd, log_file=None):
    """Run shell command and handle errors"""
    if log_file:
        with open(log_file, "w") as f:
            result = subprocess.run(cmd, shell=True, stdout=f, stderr=f)
    else:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error running command: {cmd}")
        if not log_file:
            print(f"Error output: {result.stderr}")
    return result.returncode == 0


def setup_directories():
    """Create necessary directories"""
    dirs = [
        "results/taxonomy/kraken2_output",
        "results/taxonomy/bracken_output",
        "results/taxonomy/reports",
        "results/taxonomy/merged_tables",
        "logs",
        "databases",
    ]
    for dir_path in dirs:
        Path(dir_path).mkdir(parents=True, exist_ok=True)


def extract_manual_database():
    """Extract manually downloaded Kraken2 database"""
    db_dir = Path("databases")
    tar_files = list(db_dir.glob("*.tar.gz"))

    if not tar_files:
        return None

    print(f"Found database archive: {tar_files[0]}")
    print("Extracting database...")

    try:
        with tarfile.open(tar_files[0], "r:gz") as tar:
            members = tar.getnames()
            print(f"Archive contains {len(members)} files/directories")
            tar.extractall("databases/")

        print("Extraction complete, looking for database files...")

        for root, dirs, files in os.walk(db_dir):
            root_path = Path(root)
            if "hash.k2d" in files:
                print(f"Found Kraken2 database files in: {root_path}")
                return str(root_path)

        if (db_dir / "hash.k2d").exists():
            print(f"Found Kraken2 database files in: {db_dir}")
            return str(db_dir)

        print("Database files not found after extraction")
        return None

    except Exception as e:
        print(f"Error extracting database: {e}")
        return None


def check_database():
    """Check if Kraken2 database exists or extract from manual download"""
    db_dir = Path("databases")
    print("Checking for existing database...")

    for root, dirs, files in os.walk(db_dir):
        root_path = Path(root)
        if "hash.k2d" in files and "opts.k2d" in files and "taxo.k2d" in files:
            print(f"Found complete Kraken2 database at: {root_path}")
            return str(root_path)

    print("No extracted database found, checking for archive files...")
    extracted_path = extract_manual_database()

    if extracted_path:
        return extracted_path

    print("No Kraken2 database found!")
    print(
        "Please download a database and place the .tar.gz file in the 'databases' directory"
    )
    return None


def get_available_read_lengths(db_path):
    """Check which read lengths are available in the database"""
    db_path = Path(db_path)
    kmer_files = list(db_path.glob("database*mers.kmer_distrib"))

    available_lengths = []
    for file in kmer_files:
        filename = file.name
        if "database" in filename and "mers.kmer_distrib" in filename:
            length_str = filename.replace("database", "").replace(
                "mers.kmer_distrib", ""
            )
            try:
                length = int(length_str)
                available_lengths.append(length)
            except ValueError:
                continue

    return sorted(available_lengths)


def process_individual_kraken2(r1_file, r2_file, srr_id, db_path, threads=8):
    """Run Kraken2 on an individual read pair"""
    output_dir = Path("results/taxonomy/kraken2_output/")
    output_dir.mkdir(parents=True, exist_ok=True)

    kraken_output_file = output_dir / f"{srr_id}_kraken2.out"
    kraken_report_file = output_dir / f"{srr_id}_kraken2_report.txt"

    # Check if input files exist and have content
    if not Path(r1_file).exists() or not Path(r2_file).exists():
        print(f"Input files not found for {srr_id}: {r1_file}, {r2_file}")
        return False
    if Path(r1_file).stat().st_size == 0 or Path(r2_file).stat().st_size == 0:
        print(f"Empty input files for {srr_id}: {r1_file}, {r2_file}")
        return False

    cmd = f"""kraken2 --db {db_path} \
             --threads {threads} \
             --paired {r1_file} {r2_file} \
             --output {kraken_output_file} \
             --report {kraken_report_file} \
             --use-names"""

    print(f"Running Kraken2 for sample: {srr_id}...")
    success = run_command(cmd, f"logs/kraken2_{srr_id}.log")

    if success and kraken_output_file.exists() and kraken_report_file.exists():
        # Quick validation
        with open(kraken_report_file) as f:
            lines = f.readlines()
            classified_lines = [line for line in lines if not line.startswith("U\t")]
            print(
                f"Kraken2 report has {len(classified_lines)} classified taxa for {srr_id}"
            )
        return True

    print(f"Kraken2 failed for {srr_id}")
    return False


def run_bracken(srr_id, db_path, read_length=100):
    """Run Bracken for abundance estimation for an individual sample"""
    kraken_report = f"results/taxonomy/kraken2_output/{srr_id}_kraken2_report.txt"
    bracken_output_base = f"results/taxonomy/bracken_output/{srr_id}_bracken"

    if not Path(kraken_report).exists():
        print(f"Kraken report not found: {kraken_report}")
        return False

    # Check if report has classified reads
    with open(kraken_report) as f:
        lines = f.readlines()
        if len(lines) < 2:
            print(f"Kraken report appears empty: {kraken_report}")
            return False

    levels = ["S", "G", "F", "O", "C", "P"]  # Species to Phylum
    success_count = 0

    for level in levels:
        cmd = f"""bracken -d {db_path} \
                 -i {kraken_report} \
                 -o {bracken_output_base}_{level}.bracken \
                 -w {bracken_output_base}_{level}_report.txt \
                 -r {read_length} \
                 -l {level} \
                 -t 10"""

        success = run_command(cmd, f"logs/bracken_{srr_id}_{level}.log")
        if success:
            success_count += 1

    return success_count > 0


def create_abundance_tables():
    """Create merged abundance tables from Bracken outputs"""
    bracken_dir = Path("results/taxonomy/bracken_output")
    output_dir = Path("results/taxonomy/merged_tables")
    output_dir.mkdir(parents=True, exist_ok=True)

    levels = {
        "S": "species",
        "G": "genus",
        "F": "family",
        "O": "order",
        "C": "class",
        "P": "phylum",
    }

    for level_code, level_name in levels.items():
        print(f"Creating merged abundance table for {level_name} level...")

        bracken_files = list(bracken_dir.glob(f"*_bracken_{level_code}.bracken"))
        if not bracken_files:
            print(f"No Bracken files found for {level_name} level")
            continue

        merged_data = {}
        sample_srr_ids = []

        for bracken_file in bracken_files:
            srr_id = bracken_file.stem.replace(f"_bracken_{level_code}", "")
            sample_srr_ids.append(srr_id)

            try:
                df = pd.read_csv(bracken_file, sep="\t")
                for _, row in df.iterrows():
                    taxon_name = row["name"]
                    abundance = row["new_est_reads"]

                    if taxon_name not in merged_data:
                        merged_data[taxon_name] = {}
                    merged_data[taxon_name][srr_id] = abundance

            except Exception as e:
                print(f"Error processing {bracken_file}: {e}")

        if merged_data:
            # Create abundance table
            abundance_df = pd.DataFrame(merged_data).T.fillna(0)
            abundance_df = abundance_df.reindex(columns=sorted(sample_srr_ids))

            # Save raw counts
            abundance_df.to_csv(output_dir / f"{level_name}_abundance_counts.csv")

            # Calculate relative abundances
            relative_abundance = abundance_df.copy()
            for col in abundance_df.columns:
                total = abundance_df[col].sum()
                if total > 0:
                    relative_abundance[col] = (abundance_df[col] / total) * 100
                else:
                    relative_abundance[col] = 0

            relative_abundance.to_csv(
                output_dir / f"{level_name}_relative_abundance.csv"
            )
            print(
                f"Saved {level_name} tables with {len(abundance_df)} taxa and {len(sample_srr_ids)} samples"
            )


def generate_taxonomy_summary():
    """Generate summary statistics for taxonomy results"""
    print("Generating taxonomy summary...")

    summary_stats = {
        "total_samples_processed": 0,
        "total_taxa_detected": {},
        "processing_details": {},
    }

    levels = ["species", "genus", "family", "order", "class", "phylum"]

    for level in levels:
        abundance_file = (
            f"results/taxonomy/merged_tables/{level}_relative_abundance.csv"
        )
        if Path(abundance_file).exists():
            try:
                df = pd.read_csv(abundance_file, index_col=0)
                sample_cols = list(df.columns)

                summary_stats["total_samples_processed"] = len(sample_cols)
                summary_stats["total_taxa_detected"][level] = len(df)
                summary_stats["processing_details"][level] = {
                    "samples": len(sample_cols),
                    "taxa": len(df),
                    "sample_srr_ids": sample_cols,
                }

            except Exception as e:
                print(f"Error processing {level} level summary: {e}")

    with open("results/taxonomy/taxonomy_summary.json", "w") as f:
        json.dump(summary_stats, f, indent=2)

    print("Taxonomy Summary:")
    for key, value in summary_stats.items():
        if key != "processing_details":
            print(f"  {key}: {value}")


def write_sample_category_metadata(sample_category_map):
    """Write CSV file mapping SRR_ID to category for R script"""
    output_path = Path("results/taxonomy/sample_categories.csv")
    df = pd.DataFrame(sample_category_map.items(), columns=["SampleID", "Category"])
    df = df.sort_values(by="SampleID")
    df.to_csv(output_path, index=False)
    print(f"Sample category metadata saved to {output_path}")


def main():
    """Main function to run taxonomic profiling"""
    print("Starting Taxonomic Profiling with Kraken2 Database...")

    setup_directories()

    # Check for database
    db_path = check_database()
    if not db_path:
        print("Database setup failed. Exiting.")
        return

    # Determine read length
    available_lengths = get_available_read_lengths(db_path)
    if available_lengths:
        read_length = 100 if 100 in available_lengths else available_lengths[0]
    else:
        read_length = 100
    print(f"Using read length: {read_length}")

    # Find trimmed samples
    trimmed_base_dir = Path("results/qc/trimmed")
    if not trimmed_base_dir.exists():
        print("No trimmed results directory found. Please run quality control first.")
        return

    category_dirs = [d for d in trimmed_base_dir.iterdir() if d.is_dir()]
    if not category_dirs:
        print("No category directories found in results/qc/trimmed/")
        return

    samples_to_process = []
    sample_category_map = {}

    # Collect individual FASTQ pairs
    for category_dir in category_dirs:
        category_name = category_dir.name
        r1_files = sorted(category_dir.glob("*_1_trimmed.fastq.gz"))

        for r1_file in r1_files:
            r2_file = r1_file.parent / r1_file.name.replace(
                "_1_trimmed.fastq.gz", "_2_trimmed.fastq.gz"
            )
            srr_id = r1_file.name.replace("_1_trimmed.fastq.gz", "")

            if r2_file.exists():
                samples_to_process.append(
                    {
                        "srr_id": srr_id,
                        "r1_path": str(r1_file),
                        "r2_path": str(r2_file),
                        "category": category_name,
                    }
                )
                sample_category_map[srr_id] = category_name

    if not samples_to_process:
        print("No trimmed FASTQ pairs found. Exiting.")
        return

    print(f"Found {len(samples_to_process)} samples for taxonomic profiling.")

    # Process samples
    successful_samples = 0
    for sample_info in samples_to_process:
        srr_id = sample_info["srr_id"]
        print(f"\nProcessing sample: {srr_id} (Category: {sample_info['category']})")

        # Run Kraken2 and Bracken
        if process_individual_kraken2(
            sample_info["r1_path"], sample_info["r2_path"], srr_id, db_path
        ):
            if run_bracken(srr_id, db_path, read_length):
                successful_samples += 1
                print(f"✓ Successfully processed {srr_id}")

    print(
        f"\nProcessing complete: {successful_samples}/{len(samples_to_process)} samples successful"
    )

    if successful_samples > 0:
        create_abundance_tables()
        write_sample_category_metadata(sample_category_map)
        generate_taxonomy_summary()
        print("\nTaxonomic profiling complete! Results in results/taxonomy/")
    else:
        print("\nNo samples were successfully processed. Check logs for errors.")


if __name__ == "__main__":
    main()
