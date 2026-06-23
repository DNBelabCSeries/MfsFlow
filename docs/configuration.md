# Configuration Reference

This document provides a complete reference for all configuration options in MfsFlow. Configuration can be specified via command-line arguments or a YAML configuration file.

## Table of Contents

1. [Configuration Methods](#configuration-methods)
2. [Project Configuration](#project-configuration)
3. [Sample Configuration](#sample-configuration)
4. [Sequence Files Configuration](#sequence-files-configuration)
5. [Reference Configuration](#reference-configuration)
6. [Filtering Options](#filtering-options)
7. [Barcode Options](#barcode-options)
8. [Counting Options](#counting-options)
9. [Performance Options](#performance-options)
10. [Output Options](#output-options)
11. [Example Configurations](#example-configurations)

## Configuration Methods

### Method 1: Command-Line Arguments

```bash
mfsflow \
  --fastqs /path/to/fastqs \
  --genomeDir /path/to/reference \
  --sample SAMPLE_NAME \
  --plate 1 \
  --outdir /path/to/output \
  --threads 20
```

### Method 2: YAML Configuration File

Create a YAML file and pass it to MfsFlow:

```bash
mfsflow --config /path/to/config.yaml
```

See `mfsflow/yaml/example_config.yaml` (in source repo) for a complete example.

## Project Configuration

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `project` | string | (required) | Project name, used for naming output files |

### Command-Line Equivalent
- Set via `--sample` (used as project name)

### Example
```yaml
project: "Sample01"
```

## Sample Configuration

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sample.sample_type` | string | "auto" | Sample type: "auto", "manual", or "custom" |
| `sample.sample_id` | string/int | "1" | Sample ID(s): single number for auto, comma-separated for manual |
| `barcodes.barcode_file` | string | "" | Path to custom barcode file (for `--expectBarcode` mode) |

### Command-Line Equivalent
- `--plate 1` → `sample_type: "auto"`, `sample_id: "1"`
- `--manual "20,21,22"` → `sample_type: "manual"`, `sample_id: "20,21,22"`
- `--expectBarcode /path/to/barcodes.tsv` → `sample_type: "custom"`, `barcodes.barcode_file: /path/to/barcodes.tsv`

### Custom Barcode File Format (`--expectBarcode`)

When using `--expectBarcode`, you need to provide a TSV file with custom barcode sequences.

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
A2	TGCATGCATGCATGCATGCA	ACGTTGCATGCATGCATGCA
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

### Example
```yaml
sample:
  sample_type: auto  # or "manual" or "custom"
  sample_id: 1       # or "20,21,22" for manual
```

## Sequence Files Configuration

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sequence_files.file1.name` | string | (required) | Path to first FASTQ file (R1) |
| `sequence_files.file1.base_definition` | list | [] | Base definition for file1 (e.g., `["cDNA(11-100)", "UMI(1-10)"]`) |
| `sequence_files.file2.name` | string | (required) | Path to second FASTQ file (R2) |
| `sequence_files.file2.base_definition` | list | [] | Base definition for file2 (e.g., `["cDNA(1-100)", "BC(101-120)"]`) |

### Command-Line Equivalent
- Set via `--fastqs` (directory containing FASTQ files)

### Base Definition Vocabulary
- `BC(n)`: Barcode of length n
- `UMI(n)`: UMI of length n
- `cDNA(n)`: cDNA of length n
- Ranges can be specified: `cDNA(11-100)`

### Example
```yaml
sequence_files:
  file1:
    name: /path/to/Sample01_R1.fq.gz
    base_definition:
      - "cDNA(11-100)"
      - "UMI(1-10)"
  file2:
    name: /path/to/Sample01_R2.fq.gz
    base_definition:
      - "cDNA(1-100)"
      - "BC(101-120)"
```

**Note**: Base definitions are typically auto-detected from read lengths. Manual specification is rarely needed.

## Reference Configuration

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reference.STAR_index` | string | (required) | Path to STAR genome index directory |
| `reference.GTF_file` | string | (required) | Path to gene annotation file (`.gtf` or `.gtf.gz`) |
| `reference.additional_files` | string | null | Additional reference files (optional) |
| `reference.additional_STAR_params` | string | "--clip3pAdapterSeq CTGTCTCTTATACACATCT" | Additional STAR parameters |

### Command-Line Equivalent
- Set via `--genomeDir` (expects `star/` and `genes/genes.gtf` or `genes/genes.gtf.gz`)

### Example
```yaml
reference:
  STAR_index: /path/to/reference/star
  GTF_file: /path/to/reference/genes/genes.gtf
  additional_files: null
  additional_STAR_params: "--clip3pAdapterSeq CTGTCTCTTATACACATCT"
```

### Reference Database Layout

`--genomeDir` should point to a reference root directory with this layout:

```text
/path/to/reference/
├── star/              # STAR genome index directory
└── genes/
    └── genes.gtf      # Or genes.gtf.gz; gene annotation used for feature counting
```

The pipeline then resolves:

```yaml
reference:
  STAR_index: /path/to/reference/star
  GTF_file: /path/to/reference/genes/genes.gtf
```

Compatible reference databases include:

- DNBelab C Series references generated by `dnbc4tools`, as long as they contain a STAR index and a compatible GTF.
- 10x Genomics Cell Ranger references, as long as the STAR index directory and GTF are exposed in the expected paths or copied/symlinked to the layout above.
- Custom references built with STAR from a genome FASTA and GTF.

For Cell Ranger-style references, the annotation is usually under `genes/genes.gtf`
or `genes/genes.gtf.gz`; both are supported.
If the STAR index is stored in a different subdirectory, create a symlink named `star`
or set `reference.STAR_index` directly in the YAML config.

The report displays the transcriptome name from the reference path. For example,
`/path/to/reference/star` is shown as `reference`. You can override this by adding
one of these optional fields:

```yaml
reference:
  transcriptome_name: GRCh38_GENCODE_v44
```

## Filtering Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `filter_cutoffs.BC_filter.num_bases` | int | 4 | Number of bases under quality cutoff to filter BC |
| `filter_cutoffs.BC_filter.phred` | int | 5 | Phred score cutoff for BC filtering |
| `filter_cutoffs.UMI_filter.num_bases` | int | 3 | Number of bases under quality cutoff to filter UMI |
| `filter_cutoffs.UMI_filter.phred` | int | 5 | Phred score cutoff for UMI filtering |

### Example
```yaml
filter_cutoffs:
  BC_filter:
    num_bases: 4
    phred: 5
  UMI_filter:
    num_bases: 3
    phred: 5
```

### Recommendations
- **Stringent filtering**: `num_bases: 2`, `phred: 10`
- **Relaxed filtering**: `num_bases: 5`, `phred: 3`
- **Default**: Suitable for most MGI data

## Barcode Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `barcodes.barcode_num` | int | null | Number of top barcodes to use (null = auto-detect) |
| `barcodes.automatic` | bool | false | Enable automatic barcode selection |
| `barcodes.BarcodeBinning` | int | 1 | Hamming distance for barcode binning |
| `barcodes.nReadsperCell` | int | 1 | Minimum reads per cell barcode |

### Example
```yaml
barcodes:
  barcode_num: null      # Auto-detect
  automatic: no          # Manual barcode selection
  BarcodeBinning: 1     # Hamming distance 1
  nReadsperCell: 1      # At least 1 read per cell
```

### BarcodeBinning Explanation
- `0`: No binning (exact match only)
- `1`: Bin barcodes within Hamming distance 1
- `2`: Bin barcodes within Hamming distance 2

**Recommendation**: Use `1` for most datasets. Use `0` only if you have perfect barcode sequences.

### Automatic Barcode Detection
If `barcode_num` is null and `automatic` is false, the pipeline will:
1. Count barcode frequencies
2. Identify knee point in cumulative fraction plot
3. Select barcodes above knee point

## Counting Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `counting_opts.introns` | bool | true | Count intronic reads (yes/no) |
| `counting_opts.strand` | int | 1 | Strand specificity: 0=unstranded, 1=positively stranded, 2=negatively stranded |
| `counting_opts.internal_strand` | int | 0 | Strand specificity for internal reads (override) |
| `counting_opts.Ham_Dist` | int | 1 | Hamming distance for UMI collapsing |
| `counting_opts.twoPass` | bool | false | Perform STAR two-pass mapping |
| `counting_opts.gene_body_max_reads` | int | 10000000 | Target primary mapped reads for gene body coverage per read source. Reads are selected by stable hash sampling across the whole BAM, so coordinate-sorted BAMs are not biased toward early genomic regions. Use 0 to include all reads. |
| `counting_opts.gene_body_sample_seed` | int | 42 | Seed for deterministic gene body coverage sampling |

### Example
```yaml
counting_opts:
  introns: yes
  strand: 1
  internal_strand: 0
  Ham_Dist: 1
  twoPass: no
```

### Strand Specificity
- `0`: Unstranded (e.g., Smart-seq2)
- `1`: Positively stranded (read matches gene strand)
- `2`: Negatively stranded (read is opposite to gene strand)

**MGI data**: Typically use `strand: 1`.

### Hamming Distance for UMI Collapsing
- `0`: No collapsing (exact UMI match only)
- `1`: Collapse UMIs within Hamming distance 1
- `2`: Collapse UMIs within Hamming distance 2

**Recommendation**: Use `1` for most datasets.

## Performance Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `performance_opts.stream_bc_correction` | bool | true | Enable streaming barcode correction |
| `performance_opts.tmp_root` | string | null | Temporary root directory (e.g., `/dev/shm`) |
| `performance_opts.min_free_gb` | number | 5 | Minimum free disk space required before execution |
| `performance_opts.disk_space_multiplier` | number | 4.0 | Estimated workspace requirement relative to FASTQ size |
| `num_threads` | int | 30 | Number of threads to use |

### Command-Line Equivalent
- `--threads 20` → `num_threads: 20`
- `--tmpRoot /dev/shm` → `performance_opts.tmp_root: /dev/shm`

### Example
```yaml
performance_opts:
  stream_bc_correction: true
  tmp_root: null          # or "/dev/shm"
num_threads: 30
```

### Streaming Barcode Correction
When enabled (`true`):
- Reduces intermediate file I/O
- Uses Unix pipes for data streaming
- Significantly faster for large datasets

**Recommendation**: Always enable unless debugging.

### Using `/dev/shm` for Temporary Files
If you have sufficient RAM:
```yaml
performance_opts:
  tmp_root: /dev/shm
```

This will store temporary files in memory, significantly reducing I/O latency.

**Warning**: Ensure you have enough RAM. 64GB+ recommended for large datasets.

## Output Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `out_dir` | string | (required) | Output directory path |
| `make_stats` | bool | true | Generate stats files and plots |
| `make_h5ad` | bool | true | Generate .h5ad file |
| `make_sorted_bam` | bool | true | Generate sorted BAM file |
| `make_ub_bam` | bool | false | Generate UB-corrected BAM file |
| `which_Stage` | string | "Filtering" | Start from this stage |

### Runtime Preflight and Stage State

Before execution, MfsFlow validates required Python packages, external tools,
available disk space, GTF readability, STAR index completeness, and artifacts
required by the selected resume stage.

Disk checks can be tuned in `performance_opts`:

```yaml
performance_opts:
  min_free_gb: 5
  disk_space_multiplier: 4.0
```

Each successful stage writes an artifact manifest and success marker under:

```text
XPRESS_PROCESSING/logs/stages/<Stage>.manifest.json
XPRESS_PROCESSING/logs/stages/<Stage>.success
```

Starting from `Mapping`, `Counting`, or `Summarising` fails early when required
upstream artifacts are missing or empty.

### Command-Line Equivalent
- `--outdir /path/to/output` → `out_dir: /path/to/output`
- `--stage Mapping` → `which_Stage: "Mapping"`

### Example
```yaml
out_dir: /path/to/output
make_stats: yes
make_h5ad: yes
make_sorted_bam: yes
make_ub_bam: no
which_Stage: Filtering
```

### Stage Selection
Available stages:
- `Filtering`: Start from beginning (default)
- `Mapping`: Skip filtering, start from mapping
- `Counting`: Skip filtering and mapping, start from counting
- `Summarising`: Skip to summarising stage

**Use case**: Resume failed runs from the last successful stage.

## Tool Paths

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `toolkit_directory` | string | "." | Path to MfsFlow toolkit directory |
| `read_layout` | string | "PE" | Read layout: "PE" (paired-end) or "SE" (single-end) |

### Example
```yaml
toolkit_directory: /path/to/MfsFlow
read_layout: PE
```

## Complete Example Configuration

See `mfsflow/yaml/example_config.yaml` (in source repo) for a complete example.

### Minimal Configuration (Auto Mode)
```yaml
project: "Sample01"
sample:
  sample_type: auto
  sample_id: 1
reference:
  STAR_index: /path/to/reference/star
  GTF_file: /path/to/reference/genes/genes.gtf
out_dir: /path/to/output
num_threads: 20
```

### Full Configuration (Manual Mode)
```yaml
project: "Sample01"
sample:
  sample_type: manual
  sample_id: "20,21,22,23,24"
reference:
  STAR_index: /path/to/reference/star
  GTF_file: /path/to/reference/genes/genes.gtf
  additional_STAR_params: "--clip3pAdapterSeq CTGTCTCTTATACACATCT"
out_dir: /path/to/output
num_threads: 30
filter_cutoffs:
  BC_filter:
    num_bases: 4
    phred: 5
  UMI_filter:
    num_bases: 3
    phred: 5
barcodes:
  barcode_num: null
  automatic: no
  BarcodeBinning: 1
  nReadsperCell: 1
counting_opts:
  introns: yes
  strand: 1
  internal_strand: 0
  Ham_Dist: 1
  twoPass: no
performance_opts:
  stream_bc_correction: true
  tmp_root: null
make_stats: yes
make_h5ad: yes
make_sorted_bam: yes
make_ub_bam: no
which_Stage: Filtering
toolkit_directory: /path/to/MfsFlow
read_layout: PE
```

## Next Steps

After configuring:
1. Run the pipeline (see [User Guide](user_guide.md))
2. Check output files (see [Output Files](output.md))
3. Interpret results (see [User Guide - Interpreting Reports](user_guide.md#interpreting-reports))
