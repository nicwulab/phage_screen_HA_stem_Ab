# MiniHA Antibody Phage Display Screening Analysis Pipeline

## Abstract

This repository contains a comprehensive computational pipeline for analyzing PacBio sequencing data from phage display antibody library screening experiments. The pipeline processes raw sequencing reads, identifies and annotates antibody variable region sequences, quantifies enrichment through multiple panning rounds, and generates visualizations and statistical analyses of antibody variant frequencies. The workflow is specifically optimized for analysis of heavy chain (IGH) sequences, with particular focus on IGHV1-69 and related variable region families.

## Overview

The pipeline performs end-to-end analysis of phage display screening data, including:
- Quality filtering and adapter trimming of sequencing reads
- Sequence deduplication and read counting
- Immunoglobulin gene annotation using BLAST and pyir
- Enrichment ratio calculation across panning rounds
- Sequence logo generation using Kabat numbering
- Family-level and variant-level frequency analysis
- Machine learning-based sequence evaluation (optional)

## Methodology

### Sequencing Data Processing

Circular consensus sequences (CCSs) were generated from raw PacBio subreads using SMRTLink v13.0 with the following parameters:
- Minimum accuracy requirement: 99.9%
- Minimum number of passes: 3

FASTQ format sequences were processed using BioPython's SeqIO module. Quality filtering was performed to remove reads with more than five nucleotides having Phred quality scores <40. Adapter sequences were identified and trimmed from Fab sequences, and reads lacking complete adapter sequences were excluded from further analysis.

### Sequence Processing and Annotation

Filtered reads were converted to FASTA format and deduplicated using seqkit v2.7.0. Unique sequences were aligned to reference immunoglobulin gene sequences using BLAST+ v2.12.0+ and pyir for comprehensive annotation of variable (V), diversity (D), and joining (J) gene segments.

### Frequency and Enrichment Calculations

The frequency (F) of a Fab variant *i* in sample *s* was computed for each replicate as follows:

$$F_{i,s} = \frac{readcount_{i,s} + 1}{\sum_s (readcount_{i,s} + 1)}$$

A pseudocount of 1 was added to each variant to avoid division by zero in subsequent calculations. Enrichment (E) of a Fab variant *i* in replicate *k* after phage display panning was calculated as:

$$E_{i,k} = \frac{F_{post-selection,i,k}}{F_{pre-selection,i,k}}$$

## Pipeline Workflow

### 1. Data Integration and Preprocessing

**Script:** `scripts/Fq_Integration.sh`

- Integrates multiple FASTQ files from different samples and panning rounds
- Adds sample-specific prefixes to read identifiers
- Trims flanking adapter sequences:
  - Forward adapter: `AGCTATGACCCACTCTTTCAACAGTCTTATCGTCATCG`
  - Reverse adapter: `TTCAGttcaggaggaatttaaaatgaaaaagac`
- Removes duplicate sequences using seqkit
- Generates read count statistics for total, unique, and duplicated reads

**Input:** Multiple FASTQ files specified in `Files_loc`
**Output:** 
- `PacBio/HeadTail_out.fq`: Trimmed, deduplicated sequences
- `PacBio/clean.fastqc.gz`: Compressed deduplicated sequences
- `Result/Counts_total.csv`: Total read counts per sample
- `Result/Counts_uniq.csv`: Unique read counts per sample
- `Result/Counts_duplciate.csv`: Duplicated read counts per sample

### 2. Linker Sequence Identification and Read Splitting

**Scripts:** `scripts/split_fq.py`, `scripts/split_fq2.py`

- Identifies linker sequence: `ttctagataattaattaggaggaatttaaaatgaaatacctattgcctacggcagccgctggattgttattactcgctgcccaaccagccatggcc`
- Splits reads into upstream and downstream sequences when linker is detected
- Preserves full sequences when linker is absent

**Output:** `PacBio/clean_split.fa`: Split sequences with `_up` and `_dn` suffixes

### 3. Immunoglobulin Gene Annotation

**Tool:** pyir (parallelized immunoglobulin repertoire analysis)

- Annotates sequences with V, D, J gene assignments
- Identifies CDR regions and framework sequences
- Generates amino acid translations
- Performs germline alignment

**Output:** 
- `PacBio/clean.tsv.gz`: Annotated unique sequences (TSV format)
- `PacBio/Duplicated.tsv.gz`: Annotated duplicated sequences

### 4. Heavy Chain Filtering and Quality Control

**Script:** `scripts/Remove_stopCondon_reads.py`

Filters sequences based on the following criteria:
- Locus: IGH (heavy chain) only
- Minimum amino acid length: 60 residues
- Maximum germline start position: 30 (ensures adequate 5' coverage)
- No stop codons in variable domain
- Alignment length matches germline length

**Output:** 
- `PacBio/clean_IGH.tsv.gz`: Filtered heavy chain sequences
- `PacBio/duplicated_IGH.tsv.gz`: Filtered duplicated heavy chain sequences

### 5. Amino Acid Sequence Deduplication

- Extracts amino acid sequences from annotated data
- Removes duplicate amino acid sequences
- Maintains mapping between unique sequences and original read IDs

**Output:** 
- `PacBio/all_AA.fa`: All amino acid sequences
- `PacBio/duplicated_IGH_AA.fa`: Deduplicated amino acid sequences
- `Result/duplicates_IGH_AA.txt`: Mapping of unique sequences to read IDs

### 6. Enrichment Analysis

**Scripts:** `scripts/counts_cal.py`, `scripts/duplciated_counts.py`

- Calculates frequency ratios across samples and panning rounds
- Performs exponential regression analysis for enrichment trends
- Identifies top enriched variants based on multiple criteria:
  - Exponential regression coefficients
  - Frequency ratios
  - Wild-type (wt) vs. mutant (cm) comparisons

**Output:**
- `Result/Ratio_matrix.csv`: Enrichment ratios for all variants
- `Result/ExpRe_100_list.txt`: Top 100 variants by exponential regression
- `Result/Ratio_100_list.txt`: Top 100 variants by ratio
- `Result/Wt_100_list.txt`: Top 100 wild-type variants
- `Result/Cm_100_list.txt`: Top 100 mutant variants

### 7. Family-Level Analysis

**Scripts:** `scripts/Family_seq.sh`, `scripts/Family_cal.R`

- Groups sequences by immunoglobulin variable region family (e.g., IGHV1-69, IGHV3-15, IGHV3-23)
- Calculates family-specific enrichment statistics
- Generates family-specific sequence logos and visualizations

**Output:** Family-specific directories in `Result/` containing:
- Top variant lists
- Count statistics
- Sequence logos
- Enrichment plots

### 8. Sequence Logo Generation

**Scripts:** `scripts/KABAT_SEQLOGO.py`, `scripts/Seq_log_from_tsv.py`, `scripts/I_Can_Do_all.py`

- Converts sequences to Kabat numbering scheme
- Generates sequence logos using logomaker/weblogo
- Creates position-specific amino acid frequency matrices
- Supports filtering by variable region family, sample type, and panning stage

**Usage:**
```bash
python scripts/I_Can_Do_all.py \
  -i PacBio/clean.tsv.gz PacBio/Duplicated.tsv.gz \
  -O Result/SeqLogo \
  -s wt cm \
  -p P2 P3 \
  -t 100 \
  -fv IGHV1-69 IGHV3-21
```

**Output:** Sequence logo images in PNG/SVG format

### 9. Visualization and Statistical Analysis

**Scripts:** Multiple R scripts in `scripts/` directory

- `V_class_plot.R`: Variable region family distribution plots
- `Ratio_matrix_plot.R`: Enrichment ratio heatmaps
- `Figure_1_barplot.R`: Bar plots for top variants
- `Family_cal.R`: Family-specific statistical analysis
- `co_association.R`: Co-occurrence analysis
- `Paired_aa_counts_pair.R`: Paired amino acid frequency analysis

**Output:** Visualization files in `Picture/` and `plot/` directories

## Dependencies

### Software Requirements

- **Python 3.x** with packages:
  - pandas
  - numpy
  - BioPython (SeqIO)
  - abnumber (for Kabat numbering)
  - logomaker
  - matplotlib
  - xgboost (for ML models)
  - multiprocessing

- **R** with packages:
  - ggplot2
  - dplyr
  - Additional packages as required by individual scripts

- **Command-line tools:**
  - seqkit v2.7.0
  - BLAST+ v2.12.0+
  - pyir
  - cutadapt
  - cd-hit (for clustering)
  - weblogo (optional, for sequence logos)

### System Requirements

- Multi-core processor (pipeline utilizes parallel processing)
- Sufficient RAM for large sequence datasets (recommended: 64GB+)
- Disk space: ~10-20GB for intermediate and output files

## Input Data Format

### FASTQ Files

Input FASTQ files should be named according to the pattern:
`{PanningRound}_{SampleType}_{AdditionalInfo}.fastq`

Example naming convention:
- `P0_Ab_input_F3R3-bc01_fw-bc01_rev.fastq`: Pre-selection (P0) antibody library
- `P2_cm_F4R4-bc02_fw-bc02_rev.fastq`: Panning round 2, mutant sample
- `P2_wt_F5R5-bc03_fw-bc03_rev.fastq`: Panning round 2, wild-type sample
- `P3_cm_F6R6-bc04_fw-bc04_rev.fastq`: Panning round 3, mutant sample
- `P3_wt_F7R7-bc05_fw-bc05_rev.fastq`: Panning round 3, wild-type sample

### Configuration File

`Files_loc` should contain:
- Line starting with `Fastq`: Space-separated list of input FASTQ file paths
- Line starting with `OUTPUT`: Output directory path

## Output Structure

```
Result/
├── Counts_total.csv          # Total read counts per sample
├── Counts_uniq.csv            # Unique read counts per sample
├── Counts_duplciate.csv       # Duplicated read counts per sample
├── Counts_repeat.csv          # Repeat read statistics
├── Ratio_matrix.csv           # Enrichment ratio matrix
├── ID_reads.txt               # Mapping of read IDs to sequence IDs
├── duplicates_IGH_AA.txt     # Amino acid sequence deduplication mapping
├── Du_AA_counts.csv           # Amino acid sequence counts
├── Du_AA_KABAT.csv            # Kabat-numbered amino acid sequences
├── V_class_family_count.tsv   # Variable region family counts
├── ExpRe_100_list.txt         # Top 100 by exponential regression
├── Ratio_100_list.txt         # Top 100 by enrichment ratio
├── Wt_100_list.txt            # Top 100 wild-type variants
├── Cm_100_list.txt            # Top 100 mutant variants
├── Top100ForUpload.tsv        # Top 100 sequences for external analysis
├── {Family}/                  # Family-specific analysis results
│   ├── {Family}_Alltop{N}.list
│   ├── {Family}_ReadsID_counts.csv
│   └── {Family}.tsv
└── SeqLogo/                   # Sequence logo outputs

PacBio/
├── All.fq                     # Integrated FASTQ file
├── HeadTail_out.fq            # Trimmed sequences
├── clean.fastqc.gz            # Deduplicated sequences (compressed)
├── clean.fa                   # Deduplicated sequences (FASTA)
├── clean_split.fa             # Split sequences
├── clean.tsv.gz               # Annotated unique sequences
├── clean_IGH.tsv.gz           # Filtered heavy chain sequences
├── Duplicated.tsv.gz          # Annotated duplicated sequences
├── duplicated_IGH.tsv.gz      # Filtered duplicated heavy chain sequences
├── all_AA.fa                  # All amino acid sequences
└── duplicated_IGH_AA.fa       # Deduplicated amino acid sequences
```

## Usage Examples

### Complete Pipeline Execution

Execute the main integration script:
```bash
bash scripts/Fq_Integration.sh
```

### Family-Specific Analysis

Analyze a specific variable region family:
```bash
bash scripts/Family_seq.sh
```

Modify the `FAMILY` variable in the script (e.g., `FAMILY=HV1-69`).

### Custom Sequence Logo Generation

Generate sequence logos with specific filters:
```bash
python scripts/I_Can_Do_all.py \
  -i PacBio/clean.tsv.gz PacBio/Duplicated.tsv.gz \
  -O Result/SeqLogo \
  -s Ab \
  -t 100 \
  -fv IGHV1-69
```

### Extract Top Variants

Extract top 100 sequences for external analysis:
```bash
python scripts/Top100ForUpload.py
```

## Key Statistics

### Sample Read Counts

After integration and quality filtering:

| Sample   | Total Reads | Unique Reads |
|----------|-------------|--------------|
| P0_Ab    | 375,990     | -            |
| P2_cm    | 329,090     | -            |
| P2_wt    | 307,941     | -            |
| P3_cm    | 395,761     | -            |
| P3_wt    | 355,928     | -            |
| **Total**| **1,764,710** | **1,156,148** |

### Processing Statistics

- **Adapter trimming success rate:** 99.53% (1,764,710 / 1,773,204 reads)
- **Duplication rate:** 608,562 duplicated records removed
- **Linker detection rate:** 99.87% (1,154,584 / 1,156,148 reads)
- **Sequence splitting:** 2,310,732 sequences in final library
- **Heavy chain filtering:** 942,051 unique sequences after quality control

## Additional Tools

### Machine Learning Evaluation

The pipeline includes optional machine learning models for sequence evaluation:
- **XGBoost model:** `Result/xgb_model.json`
- **Usage:** See `fasta_to_KABAT.py` for example implementation

### Sequence Clustering

CD-HIT clustering for sequence similarity analysis:
```bash
cd-hit -i Result/Dupli_aa.fa \
  -o Result/cluster/Dupli_aa0.8.fa \
  -c 0.8 \
  -M 32000 \
  -T 8
```

## Citation

If using this pipeline, please cite:
- seqkit: Shen et al., PLoS ONE 2016
- BLAST+: Camacho et al., BMC Bioinformatics 2009
- pyir: (please cite appropriate reference)

## Contact and Support

For questions or issues regarding this pipeline, please contact the development team or refer to the individual script documentation.
