# User Guide

This guide provides detailed instructions for using the MfsFlow toolkit to process high-sensitivity full-length transcriptome data.

## Table of Contents

1. [Workflow Overview](#workflow-overview)
2. [Input Data Preparation](#input-data-preparation)
3. [Running the Pipeline](#running-the-pipeline)
4. [Understanding Barcode Modes](#understanding-barcode-modes)
5. [Interpreting Reports](#interpreting-reports)
6. [Advanced Configuration](#advanced-configuration)

## Workflow Overview

MfsFlow implements a four-stage pipeline:

```
FASTQ Files → Filtering → Mapping → Counting → Summarising → Reports
```

### Stage 1: Filtering
- Quality control of reads based on Phred scores
- Barcode and UMI extraction
- Read filtering based on quality cutoffs
- FASTQ splitting for parallel processing

### Stage 2: Mapping
- STAR alignment to reference genome
- Supports both UMI and internal read mapping
- Streaming barcode correction (optional)
- BAM file sorting and indexing

### Stage 3: Counting
- FeatureCounts for read quantification
- UMI deduplication (Hamming distance-based)
- Support for exon and intron counting
- Downsampling for saturation analysis

### Stage 4: Summarising
- Quality control metrics calculation
- Saturation curve generation
- Gene body coverage analysis
- HTML report generation (auto and manual modes)

## Input Data Preparation

### Supported Data Types

MfsFlow supports paired-end FASTQ files from MGI sequencing platforms.

### Read Structure Requirements

The pipeline supports two read structure modes:

#### Mode 1: R2 Longer than R1 (Automatic Barcode Detection)
- R2 is 20bp longer than R1
- Last 20bp of R2 are treated as the read barcode
- No samplesheet required

Example:
```
R1: 100bp
R2: 120bp (100bp read + 20bp barcode)
```

#### Mode 2: Equal-Length R1 and R2 (Samplesheet Required)
- R1 and R2 have the same length
- Barcode must be provided in a samplesheet
- Each FASTQ pair is assigned a specific barcode

Samplesheet format (`samplesheet.csv`):
```csv
read1,read2,barcode
S1_R1.fq.gz,S1_R2.fq.gz,ACGTACGTACGTACGTACGT
S2_R1.fq.gz,S2_R2.fq.gz,TGCATGCATGCATGCATGCA
```

### Reference Genome Preparation

#### Directory Structure
```
/path/to/reference/
├── star/                   # STAR index directory
│   ├── Genome              # STAR genome index
│   ├── SA                  # STAR suffix array
│   └── ...
└── genes/
    └── genes.gtf           # Or genes.gtf.gz; gene annotation file
```

#### Building STAR Index

Download reference files from [GENCODE](https://www.gencodegenes.org):

```bash
# Download FASTA and GTF files
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

# Decompress
gunzip gencode.v44.primary_assembly.fa.gz
gunzip gencode.v44.annotation.gtf.gz

# Create directories
mkdir -p /path/to/reference/star
mkdir -p /path/to/reference/genes

# Copy GTF file
cp gencode.v44.annotation.gtf /path/to/reference/genes/genes.gtf

# Build STAR index
STAR --runMode genomeGenerate \
     --genomeDir /path/to/reference/star \
     --genomeFastaFiles /path/to/reference/gencode.v44.primary_assembly.fa \
     --sjdbGTFfile /path/to/reference/genes/genes.gtf \
     --runThreadN 20
```

## Running the Pipeline

### Basic Command

```bash
mfsflow \
  --fastqs /path/to/fastq_dir \
  --genomeDir /path/to/reference \
  --sample SAMPLE_NAME \
  --plate 1 \
  --outdir /path/to/output \
  --threads 20
```

### Command-Line Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--fastqs` | Yes | Directory containing input R1/R2 FASTQ files |
| `--genomeDir` | Yes | Reference directory containing `star/` and `genes/genes.gtf` or `genes/genes.gtf.gz` |
| `--sample` | Yes | Sample name (used for output naming) |
| `--outdir` | No | Output directory (default: `./<sample_name>`) |
| `--threads` | No | Number of threads (default: 30) |
| `--tmpRoot` | No | Temporary root for intermediate files (e.g., `/dev/shm`) |
| `--stage` | No | Start from specific stage (Filtering/Mapping/Counting/Summarising) |
| `--resume` | No | Reuse the saved run configuration; requires `--outdir` and an explicit `--stage` |
| `--plate` | No | Plate ID for automatic barcode mode (mutually exclusive with `--manual`, `--expectBarcode`) |
| `--manual` | No | Manual barcode IDs (comma-separated, e.g., `"20,21,22"`) (mutually exclusive with `--plate`, `--expectBarcode`) |
| `--expectBarcode` | No | Path to custom barcode file (mutually exclusive with `--plate`, `--manual`) |
| `--samplesheet` | No** | CSV samplesheet for equal-length R1/R2 data |

When none of `--plate`, `--manual`, or `--expectBarcode` is specified, the pipeline automatically runs barcode discovery mode, inferring the barcode set from observed reads.
**Required for equal-length R1/R2 data.

### Example Commands

#### Example 1: Plate Mode (Automatic Barcodes)
```bash
mfsflow \
  --fastqs /data/fastq/Sample01 \
  --genomeDir /ref/human_GRCh38 \
  --sample Sample01 \
  --plate 1 \
  --outdir /output/Sample01 \
  --threads 40
```

#### Example 2: Manual Mode (Custom Barcode IDs)
```bash
mfsflow \
  --fastqs /data/fastq/Sample02 \
  --genomeDir /ref/mouse_GRCm39 \
  --sample Sample02 \
  --manual "20,21,22,23,24" \
  --outdir /output/Sample02 \
  --threads 20
```

#### Example 3: Custom Barcode File
```bash
mfsflow \
  --fastqs /data/fastq/Sample03 \
  --genomeDir /ref/human_GRCh38 \
  --sample Sample03 \
  --expectBarcode /data/barcodes/custom_barcodes.tsv \
  --outdir /output/Sample03 \
  --threads 20
```

Custom barcode file format (`custom_barcodes.tsv`):

**Format**: TSV (tab-separated)

**Columns**:
| Column | Name | Description |
|--------|------|-------------|
| 1 | `wellID` | Well ID (e.g., A1, B2) |
| 2 | `umi_barcode` | UMI barcode sequence(s), multiple separated by commas |
| 3 | `internal_barcode` | Internal barcode sequence(s), multiple separated by commas |

**Example 1: With header row (recommended)**
```tsv
wellID	umi_barcode	internal_barcode
A1	ACGTACGTACGTACGTACGT	TGCAACGTACGTACGTACGT
A2	TGCATGCATGCATGCATGCA	ACGTTGCATGCATGCATGCA
```

**Example 2: Without header row**
```tsv
A1	ACGTACGTACGTACGTACGT	TGCAACGTACGTACGTACGT
A2	TGCATGCATGCATGCATGCA	ACGTTGCATGCATGCATGCA
```

**Example 3: Multiple barcodes per well (comma-separated)**
```tsv
wellID	umi_barcode	internal_barcode
A1	ACGTACGTACGTACGTACGT,TGCAACGTACGTACGTACGT	TGCAACGTACGTACGTACGT,ACGTTGCATGCATGCATGCA
A2	TGCATGCATGCATGCATGCA	ACGTTGCATGCATGCATGCATGCA
```

**Notes**:
- Use **tab** as delimiter (not spaces or commas)
- Barcode sequences are case-insensitive (automatically converted to uppercase)
- Multiple barcodes in one column are separated by commas
- Header row is optional (if present, it will be automatically skipped)
- File extension can be `.tsv`, `.txt`, or any text format

**Creating the file**:
```bash
# Create using echo (note the -e flag and \t for tab)
echo -e "wellID\tumi_barcode\tinternal_barcode" > custom_barcodes.tsv
echo -e "A1\tACGTACGTACGTACGTACGT\tTGCAACGTACGTACGTACGT" >> custom_barcodes.tsv
echo -e "A2\tTGCATGCATGCATGCATGCA\tACGTTGCATGCATGCATGCA" >> custom_barcodes.tsv

# Verify the file
cat custom_barcodes.tsv
```

#### Example 4: Barcode Discovery Mode
```bash
mfsflow \
  --fastqs /data/fastq/Sample04 \
  --genomeDir /ref/human_GRCh38 \
  --sample Sample04 \
  --outdir /output/Sample04 \
  --threads 20
```

When no barcode mode is specified, the pipeline automatically runs barcode discovery —
it infers the barcode set from observed reads, writes a discovery report, then continues.

#### Example 5: Equal-Length R1/R2 with Samplesheet
```bash
mfsflow \
  --fastqs /data/fastq/Sample05 \
  --samplesheet /data/samplesheet/Sample05.csv \
  --genomeDir /ref/human_GRCh38 \
  --sample Sample05 \
  --plate 1 \
  --outdir /output/Sample05 \
  --threads 20
```

#### Example 6: Resume from Specific Stage
```bash
# Resume from Mapping stage (skip filtering)
mfsflow \
  --outdir /output/Sample01 \
  --resume \
  --stage Mapping
```

#### Example 7: Using Memory-Backed Temporary Directory
```bash
mfsflow \
  --fastqs /data/fastq/Sample01 \
  --genomeDir /ref/human_GRCh38 \
  --sample Sample01 \
  --plate 1 \
  --outdir /output/Sample01 \
  --threads 20 \
  --tmpRoot /dev/shm
```

**Note**: Using `/dev/shm` can significantly speed up I/O-intensive operations but requires sufficient RAM.

## Understanding Barcode Modes

MfsFlow supports four barcode selection modes:

### 1. Plate Mode (`--plate`)

- Uses pre-defined barcode lists for standard plates
- Barcode IDs correspond to well positions (1-12 for auto, 1-24 for manual)
- Strict barcode detection: fails if expected barcodes are not observed

Use case: Standard plate-based experiments with known barcode sets.

### 2. Manual Mode (`--manual`)

- Specify custom barcode IDs (comma-separated)
- Uses manual barcode lists
- Strict barcode detection

Use case: Custom experiments with specific barcode subsets.

### 3. Custom Barcode File (`--expectBarcode`)

- Provide a custom TSV file with barcode sequences
- Most flexible option
- Strict barcode detection

Use case: Completely custom barcode designs.

### 4. Barcode Discovery Mode (default)

- Automatically infers barcode set from observed reads (triggered when no `--plate`/`--manual`/`--expectBarcode` is specified)
- Compares observed barcodes against bundled lists (exact/Hamming-1 matching)
- Writes discovery report before proceeding
- Less strict than other modes

Use case: Unknown or complex barcode configurations, or when you are not sure which
plate/manual ID to use. This is the default behavior and is often the best first-pass
mode for new datasets because it gives a transparent barcode inference report instead
of requiring you to guess the plate/manual ID upfront.

## Interpreting Reports

MfsFlow generates two types of HTML reports:

### 1. Automatic Report (`<project>_auto_plate_report.html`)

- Overview statistics
- Read distribution plots
- QC metrics
- Suitable for quick quality assessment

### 2. Manual Report (`<project>_manual_report.html`)

- Detailed per-well statistics
- UMI-based vs read-based statistics
- Per-well QC distribution
- Saturation curves
- Gene body coverage plots
- Suitable for in-depth analysis

### Key Metrics

#### Sample Overview
- **Total Reads**: Number of read pairs for PE data, or reads for SE data
- **Mapped Reads**: Number of read pairs/reads successfully mapped to reference
- **Mapping Ratio**: Percentage of mapped read pairs/reads
- **UMI Fraction**: Fraction of read pairs/reads containing valid UMIs
- **Exon+Intron Ratio**: Fraction of read pairs/reads mapping to exons and introns

#### UMI-Based vs Read-Based Statistics
- **Genes Detected**: Number of genes with at least one read/UMI
- **UMIs per Gene**: Median number of UMIs per gene
- **Reads per UMI**: Median number of read pairs/reads per UMI (indicates amplification bias)

#### Read Distribution
- **Exonic**: Reads mapping to exons
- **Intronic**: Reads mapping to introns
- **Intergenic**: Reads mapping to intergenic regions
- **Unmapped**: Reads that failed to map

#### Saturation Analysis
- Shows gene detection saturation as a function of sequencing depth
- Helps determine if more sequencing is needed

#### Gene Body Coverage
- Shows 5' to 3' coverage bias
- Full-length transcripts should show uniform coverage

## Advanced Configuration

For advanced configuration options, see the [Configuration Reference](configuration.md).

### Using the Generated YAML Configuration

MfsFlow writes the complete run configuration to
`XPRESS_PROCESSING/config/run_config.yaml` before execution. Use the normal
CLI arguments to create a run; the generated YAML preserves the exact paths
and options used for that run. The current CLI does not accept a separate
`--config` argument.

### Performance Optimization

1. **Increase threads**: Use `--threads` to specify more CPU cores
2. **Use fast storage**: SSD or RAM-backed storage for intermediate files
3. **Enable streaming**: Set `stream_bc_correction: true` in config
4. **Adjust downsampling**: Modify downsampling levels based on your needs

### Troubleshooting

For common issues and solutions, see the [Troubleshooting Guide](troubleshooting.md).

## Next Steps

After running MfsFlow:

1. Review the HTML reports in `outs/` directory
2. Load the `.h5ad` file into Scanpy or Seurat for downstream analysis
3. Use MEX matrices with your preferred single-cell analysis tools

Example (Python/Scanpy):
```python
import scanpy as sc

# Load the h5ad file
adata = sc.read_h5ad("output/XPRESS_PROCESSING/expression/SAMPLE.h5ad")

# Basic QC
sc.pp.calculate_qc_metrics(adata)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
```
