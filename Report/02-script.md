# Quality Control Script Analysis

## 1. Script Purpose
This script checks the quality of DNA sequencing data and cleans it up before analysis. It removes poor-quality sequences and adapters (leftover sequencing chemicals) that could interfere with identifying bacteria.

## 2. Key Functions

### What the script does step-by-step:
1. **Quality Check**: Examines raw sequencing files to see how good the data is
2. **Data Cleaning**: Removes low-quality sequences and unwanted adapter sequences
3. **Second Quality Check**: Verifies that the cleaning worked properly
4. **Summary Reports**: Creates easy-to-read reports showing what was done

### Processing approach:
- Works on one sample at a time
- Handles paired DNA reads (forward and reverse sequences)
- Uses 8 computer processors simultaneously for faster processing

## 3. Input/Output

### What goes in:
- DNA sequencing files stored in `fastq_files/` folder
- Each sample has its own subfolder
- Files are compressed and named with `_1` and `_2` for paired reads

### What comes out:
- Quality reports for original data
- Cleaned DNA sequence files
- Quality reports for cleaned data  
- Individual processing reports for each sample
- Combined summary report
- Statistics table showing how much data was kept
- Log files tracking any problems

## 4. Methods Used

### Software tools:
- **FastQC**: Creates quality reports with graphs and statistics
- **fastp**: Does the actual cleaning and trimming of sequences
- **MultiQC**: Combines all reports into one easy-to-view dashboard

### Cleaning standards applied:
- Keeps sequences with quality score ≥20 (99% accuracy)
- Removes reads if more than 20% of bases are low quality
- Discards sequences shorter than 100 base pairs
- Automatically finds and removes adapter contamination
- Fixes sequencing errors where possible

## 5. Results Generated

### Key measurements:
- How many DNA sequences started with vs. kept after cleaning
- Percentage of high-quality bases (Q30 rate)
- How many adapter sequences were removed
- What percentage of data was retained for each sample

### Reports created:
- Visual quality reports showing graphs of sequence quality
- Before/after comparison charts
- Summary table with statistics for all samples
- Interactive web-based dashboard for easy viewing
