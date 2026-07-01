# MfsFlow

**MfsFlow** is a comprehensive and flexible pipeline designed for processing high-sensitivity full-length transcriptome data. It is optimized for MGI sequencing chemistry and provides an end-to-end workflow from raw FASTQ data to gene expression quantification.

## 🌟 Features

- **End-to-End Pipeline**: Complete workflow including filtering, mapping, counting, and statistical analysis
- **Dual Report Modes**: Automatic and manual report generation for different analysis scenarios
- **Flexible Barcode Modes**: Support for plate-based, manual, and discovery barcode modes
- **Efficient Processing**: Streamed BAM I/O between stages, single-pass FASTQ splitting, and in-stream barcode correction
- **Detailed Reporting**: Comprehensive HTML reports with interactive plots and QC metrics
- **MEX & H5AD Output**: Standard MEX format and AnnData (H5AD) for seamless downstream analysis

## 📋 Prerequisites

Ensure you have the following installed on your system:

### System Requirements
- **Python 3.8+**
- **Linux/Unix** environment (recommended for production)

### Bundled Tools
MfsFlow comes with bundled bioinformatics tools in the `software/` directory:
- **STAR**: For genome alignment
- **Samtools**: For BAM file manipulation
- **Pigz**: For parallel gzip compression
- **Seqkit**: For FASTQ manipulation
- **FeatureCounts** (from Subread): For read counting

The pipeline automatically resolves these tools from the `software/` directory on Linux systems. No manual installation required.

### Python Dependencies
Install via pip:
```bash
pip install -r requirements.txt
```

Required Python packages:
- PyYAML>=5.1
- pysam>=0.21
- numpy>=1.20
- pandas>=1.3
- scipy>=1.7
- matplotlib>=3.5
- anndata>=0.8
- h5py>=3.1
- Pillow>=9.0

**Note**: All packages use open version ranges for maximum compatibility. If you encounter dependency conflicts, you can adjust versions as needed.

## 🚀 Installation

### 1. Clone the Repository
```bash
git clone https://github.com/lishuangshuang0616/MfsFlow.git
cd MfsFlow
```

### 2. Install Python Dependencies
```bash
pip install -r requirements.txt
```

### 3. Install in Development Mode (Optional)
```bash
pip install .
```

This will install the `mfsflow` command-line tool.

### 4. Set Up External Tools
On **Linux**, the bundled tools under `software/` are prepended to `PATH` automatically — no manual installation is needed. 

Reference databases can be built with STAR directly, or prepared from compatible
DNBelab C Series (`dnbc4tools`) and Cell Ranger references. The expected layout is
`<genomeDir>/star` for the STAR index and `<genomeDir>/genes/genes.gtf` or
`<genomeDir>/genes/genes.gtf.gz` for the GTF.

## 📚 Quick Start

### Basic Usage

Run the analysis pipeline with command-line parameters:

```bash
python3 run_analysis_pipeline.py \
  --fastqs /path/to/fastq_dir \
  --genomeDir /path/to/reference \
  --sample SAMPLE_NAME \
  --outdir /path/to/output \
  --threads 20
```

Or use the installed CLI:

```bash
mfsflow \
  --fastqs /path/to/fastq_dir \
  --genomeDir /path/to/reference \
  --sample SAMPLE_NAME \
  --outdir /path/to/output \
  --threads 20
```

### Input Modes

The pipeline supports two input modes based on read structure:

1. **R2 is 20bp longer than R1**: The last 20bp of R2 are treated as the read barcode
2. **R1 and R2 have the same length**: Provide `--samplesheet`; barcode is assigned per FASTQ pair

#### Samplesheet Format (Equal-Length R1/R2)

```csv
read1,read2,barcode
A_R1.fq.gz,A_R2.fq.gz,ACGTACGTACGTACGTACGT
B_R1.fq.gz,B_R2.fq.gz,TGCATGCATGCATGCATGCA
```

### Barcode Modes

Choose one of the following barcode selection modes:

```bash
# Plate mode (automatic barcode list)
--plate 1

# Manual mode (specify barcode IDs)
--manual 20,21,22,23,24

# Custom barcode file
--expectBarcode barcodes.tsv

# No mode specified → barcode discovery (default)
# Automatically infers plate/manual barcode set from observed reads
```

**Note**: When `--manual`, `--plate`, or `--expectBarcode` is specified, barcode detection is strict — if none of the expected barcodes are observed, the run stops. When no mode is specified, the pipeline automatically runs barcode discovery: it compares observed read barcodes against the bundled manual/plate barcode lists with exact/Hamming-1 matching, writes a discovery report, then continues with the inferred barcode set.

## 📖 Documentation

Comprehensive documentation is available in the `docs/` directory:

- [**User Guide**](docs/user_guide.md): Detailed usage instructions and workflow explanation
- [**Installation Guide**](docs/installation.md): Step-by-step installation instructions
- [**Configuration Reference**](docs/configuration.md): Complete configuration options
- [**Output Files**](docs/output.md): Description of all output files
- [**Result Usage**](docs/output.md#expression-files): How to read H5AD and MEX outputs
- [**Troubleshooting**](docs/troubleshooting.md): Common issues and solutions
- [**Architecture**](docs/architecture.md): Technical architecture and design

## 📁 Output Structure

Results are saved under your specified `out_dir`:

```
output_directory/
├── XPRESS_PROCESSING/
│   ├── config/              # Generated run configuration and expected barcode tables
│   ├── logs/               # Pipeline logs
│   ├── barcodes/           # Kept barcodes, barcode binning, and discovery report
│   ├── expression/         # MEX matrices and <project>.h5ad
│   ├── stats/              # QC tables, Q30 statistics, and plots
│   ├── intermediate/       # Temporary BAM and split FASTQ working files
│   └── *.bam              # Final UB-corrected sorted BAM
└── outs/                   # Final reports and customer-facing outputs
    ├── auto_report.html    # Automatic report
    ├── manual_report.html  # Manual report
    ├── expression/         # H5AD and MEX matrices
    ├── stats/              # QC tables and plots
    ├── bam/                # Final UB-corrected BAM and index
    └── config/             # Run config and expected barcode table
```

MfsFlow supports resuming interrupted runs from intermediate stages. After a run
finishes and reports are generated, final deliverables are moved into `outs/`;
use `outs/` as the completed result directory.

Before each run, MfsFlow checks Python dependencies, required external tools,
disk space, reference integrity, and resume-stage inputs. Successful stages
write artifact manifests and `.success` markers under
`XPRESS_PROCESSING/logs/stages/`.

## 🔧 Advanced Usage

### Resume from a Specific Stage

```bash
mfsflow \
  --fastqs /path/to/fastq_dir \
  --genomeDir /path/to/reference \
  --sample SAMPLE_NAME \
  --plate 1 \
  --outdir /path/to/output \
  --threads 20 \
  --stage Mapping
```

Available stages: `Filtering`, `Mapping`, `Counting`, `Summarising`

### Use Temporary Directory in Memory

```bash
mfsflow \
  --fastqs /path/to/fastq_dir \
  --genomeDir /path/to/reference \
  --sample SAMPLE_NAME \
  --plate 1 \
  --outdir /path/to/output \
  --threads 20 \
  --tmpRoot /dev/shm
```


## 📊 Reference Genome Preparation

### Build STAR Index

Download the genome sequence (FASTA) and gene annotation (GTF) files from [GENCODE](https://www.gencodegenes.org).

```bash
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles /path/to/genome.fa \
     --sjdbGTFfile /path/to/genes.gtf \
     --runThreadN 20
```

The pipeline expects STAR index at `<genomeDir>/star` and GTF at `<genomeDir>/genes/genes.gtf` or `<genomeDir>/genes/genes.gtf.gz`.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📧 Contact

For questions, bug reports, or feature requests, please open an issue on GitHub.

---

**Version**: 1.0.0  
**Last Updated**: 2026-07-01
