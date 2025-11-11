![cover](images/00-coverpage.png)

## 📑 Table of contents

1. [📝 Overview](#-overview)
2. [📦 Installation](#-installation)
3. [🏃 How to run](#-how-to-run)
4. [📂 Input file formats & expected structure](#-input-file-formats--expected-structure)
5. [🔄 App workflows and usage](#-app-workflows-and-usage)
   - [🧬 Genome Indexing Workflow](#-genome-indexing-workflow)
   - [🔬 Bulk RNA-Seq Analysis Workflow](#-bulk-rna-seq-analysis-workflow)
   - [📊 Data Visualization Workflow](#-data-visualization-workflow)
6. [📚 Citation](#-citation)
7. [🤝 Acknowledgments](#-acknowledgments)
8. [💡 Inspiration](#-inspiration)


## 📝 Overview
<p align="justify"> <b>Bulk RNAlyzer</b> is a comprehensive web-based application built with Streamlit that provides a complete workflow for analyzing and visualizing bulk RNA-seq data. This user-friendly interface allows researchers to perform end-to-end RNA-seq analysis without requiring extensive command-line expertise. </p>

#### 🛠️ Pipeline Overview

| Description | Tools Used |
|--------------|-------------|
| Quality Control | FastQC |
| Reads Trimming | Trimmomatic |
| Alignment to Reference Genome | HISAT2 |
| Quantification | StringTie |
| Differential Expression Analysis | DESeq2 |
| Data Visualizations | ComplexHeatmap |

#### ✨ Key Features

- User-Friendly Interface: No command-line expertise required.
- Flexible Starting Points: Begin bulk RNA-seq analysis from any step in the pipeline.
- Interactive Visualizations: Generate publication-quality heatmaps, boxplots, volcano plots, and PCA plots.
- Comprehensive Output: Provides detailed reports and summary statistics.
- Multi-Condition Comparison: Supports comparison of multiple combinations of conditions when more than two conditions are present.
- Normalized Data Output: Produces normalized gene and transcript count matrices.
- Detailed DEG Results: Outputs overall DEG CSV, gene summary CSV, and gene biotype-specific DEG CSV files, along with .rds objects and associated visualizations.

## 📦 Installation

#### Create conda environment

conda create -n bulkrnalyzer python=3.11

conda activate bulkrnalyzer

#### Install tools for bulk RNA-seq analysis

#### 1. Fastqc:

```
sudo apt install fastqc
```

#### 2. Trimmomatic

Download binary from official website

#### 3. HISAT2

- Using the `HISAT2` officical [link](https://daehwankimlab.github.io/hisat2/download/#version-hisat2-221), download the binary or install from [source code](https://github.com/DaehwanKimLab/hisat2)
- After downloading/installing add `HISAT2` directory to System's PATH:

1. Open nano:

"""
nano ~/.bashrc
"""

2. Add path of `HISAT2` directory:

"""
export PATH=$PATH:/home/yourusername/tools/hisat2-2.2.1
"""

3. Save changes and exit:

- `Ctrl + O` Save
- Will be prompted for filename, press `Enter`
- `Ctrl + X` Exit

4. Reload shell:

"""
source ~/.bashrc
"""

5. Verify by running:

```
hisat2-build
```

#### 4. Samtools

"""
sudo apt install samtools
"""

#### 5. Sambamba

"""
sudo apt install sambamba
"""

#### 6. gffread

```
sudo apt install gffread
```

#### 7. StringTie

"""
sudo apt install stringtie
"""

#### 🐍 Install Python Dependencies

```
pip install streamlit pandas numpy pathlib subprocess -y
```

#### Install R Packages

- Open terminal and run:

```
Rscript install_rpackages.R
```

## 🏃 How to run

**Navigate to the directory containing the `main.py` file and run in terminal:**

```
streamlit run main.py
```

The application will open in your default web browser at `http://localhost:8501`

## 📂 Input file formats & expected structure

#### 🔬 Bulk RNA-Seq Analysis Inputs

🧬 Raw Sequencing Data

- Format: FASTQ files (compressed .fastq.gz preferred)

🏷️ Sample Metadata

- Format: CSV
- Required Columns: Sample identifiers (matching FASTQ file names), Condition/treatment groups etc
- All column names should be lowercase
- Example Strcuture:

| sample | condition | time_point | batch |
|-----------|-----------|------------|-------|
| sample1 | treated | 24h | batch1 |
| sample2 | control | 24h | batch1 |

📚 Reference Files

- Genome Assembly: FASTA format (.fa, .fasta)
- Gene Annotations: GTF or GFF3 format
- Adapter Sequences: FASTA format

#### 🎯 Visualization Inputs

📊 Count Matrices

- Gene Count Matrix: CSV with genes as rows, samples as columns
- Transcript Count Matrix: CSV with transcripts as rows, samples as columns
- Required Columns: gene_id, gene_symbol (for gene matrix) or transcript_id, transcript_name (for transcript matrix)

🎨 Visualization Parameters

- Gene Lists: Comma-separated gene names
- Plot Dimensions: Customizable width and height in pixels
- Color Palettes: Multiple predefined color schemes available for heatmap

## 🔄 App Workflows and Usage

**Add image of streamlit app**

The application provides three main analysis types:

- 🧬 Index Genome - Build reference genome indexes
- 🔬 Bulk RNA-Seq Analysis - Complete pipeline from raw reads to differential expression
- 📊 Data Visualizations - Create interactive plots from analysis results

1. 🧬 Genome Indexing Workflow

- Create custom HISAT2 reference genome indexes for alignment, including support for transcript-aware and SNP-aware indexing.

Select Index Type:

- HFM (basic genome index)
- HGFM with transcripts
- HGFM with SNPs
- HGFM with SNPs and transcripts

Input Files:

- Genome FASTA file
- Annotation file (GTF/GFF3) for transcript-aware indexing
- VCF file for SNP-aware indexing

Parameters:

- Number of threads
- Large index option (for genomes >4 billion bases)
- Offrate value (controls index density)
- Output: HISAT2 index files (.ht2)


2. 🔬 Bulk RNA-Seq Analysis Workflow

Start analysis from any step in the pipeline:

| Steps | Description           | Input Requirements |
|------|-----------------------|--------------------|
| 1️⃣ | **QC Before Trimming** | Raw FASTQ files |
| 2️⃣ | **Reads Trimming** | Raw FASTQ files, Adapter sequences |
| 3️⃣ | **QC After Trimming** | Trimmed FASTQ files |
| 4️⃣ | **Alignment** | Trimmed FASTQ files, HISAT2 index |
| 5️⃣ | **Quantification** | BAM files, Annotation file |
| 6️⃣ | **DESeq2 Analysis** | Count matrices, Metadata file |

Step 1: Quality Control (Before Trimming) 🔍

- Tool: FastQC
- Output: HTML reports with quality metrics
- Assess raw read quality before processing

Step 2: Read Trimming ✂️

- Tool: Trimmomatic
- Parameters:
   - PHRED encoding (33/64)
   - Leading/Trailing quality thresholds
   - Optional: Sliding window trimming, minimum length
   - Adapter Selection: Built-in Trimmomatic adapters or custom adapter sequences
- Output: Trimmed FASTQ files (*_paired.fastq.gz, *_unpaired.fastq.gz)

Step 3: Quality Control (After Trimming) 🔍

- Tool: FastQC
- Verify trimming effectiveness and final read quality
- Output: Post-trimming quality reports

Step 4: Alignment 🧬

- Tool: HISAT2 with SAMtools and Sambamba
- Parameters:
   - Number of threads
   - Maximum alignments per read
   - Intron length boundaries
   - Keep intermediate files option

- Process:
   - HISAT2 alignment to generate SAM files
   - SAM to sorted BAM conversion (SAMtools)
   - Duplicate marking (Sambamba)
   - BAM indexing (SAMtools)
- Output: Deduplicated, sorted BAM files with indexes

Step 5: Quantification 📊

- Tool: StringTie with prepDE.py
- Parameters:
   - Minimum transcript coverage
   - Minimum isoform abundance
   - Junction coverage
   - Annotation Handling: Automatic GFF3 to GTF conversion if needed
- Output:
   - Gene count matrix (gene_count_matrix.csv)
   - Transcript count matrix (transcript_count_matrix.csv)

Step 6: Differential Expression Analysis 📈

- Tool: DESeq2 (R package)
- Parameters:
   - Statistical test (Wald or LRT)
   - Log fold change threshold
   - Significance level (alpha)
   - Fit type (parametric, local, mean)
- Comparison Setup:
   - Wald test: Define condition contrasts
   - LRT test: Select model combinations
- Output:
   - Normalized count matrices
   - Differential expression results
   - PCA plots, volcano plots, heatmaps
   - Biotype-specific results

3. 📊 Data Visualization Workflow

🔥 Heatmaps

- Features:
   - Z-score normalized expression
   - Pre-defined color schemes
   - Row and column annotations
   - Publication-ready resolution

📦 Boxplots

- Features:
   - Individual gene/transcript expression
   - Condition-wise comparisons
   - Pre-defined color palettes
   - Multiple genes in separate plots

**Add visualization images**

## 📚 Citation

If you use Bulk RNAlyzer in your research, please cite:

## 🙏 Acknowledgements

- Bulk RNAlyzer builds upon the work of open-source bioinformatics tools and libraries.
- Please also cite the tools used in the analysis pipeline:
   - HISAT2 - For efficient and sensitive alignment of RNA-seq reads
   - FastQC - For comprehensive quality control of sequencing data
   - Trimmomatic - For flexible trimming of sequence data
   - StringTie - For accurate transcript quantification and assembly
   - DESeq2 - For robust differential expression analysis
   - SAMtools - For handling high-throughput sequencing data
   - ComplexHeatmap - For creating sophisticated heatmaps
   - ggplot2 - For elegant data visualization
   - EnhancedVolcano - For publication-ready volcano plots
   - Streamlit - For creating interactive web applications with Python
   - Pandas - For data manipulation and analysis
   - Other open-source tools

## 💡 Inspiration

This application was inspired by the need for:

- Accessible bioinformatics for researchers without computational backgrounds
- Reproducible workflows in transcriptomics analysis
- Integrated solutions that bridge multiple analysis steps
- Visual analytics for better data interpretation
