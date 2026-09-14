# mfs-stitcher

This directory contains an optional Rust implementation of the molecule
stitching and isoform-indexing algorithms. It is deliberately kept separate
from the Python package while result-parity and large-data benchmarks are being
completed.

The production `mfsflow` command continues to use the Python implementation.
The Rust crate is currently validated by CI and can be built or tested with:

```bash
cargo test --locked --manifest-path crates/mfs-stitcher/Cargo.toml
cargo run --release --manifest-path crates/mfs-stitcher/Cargo.toml -- --help
```

The crate accepts BAM/GTF inputs directly and supports optional isoform index,
MEX, counts TSV, and molecule TSV output. Do not substitute it for the Python
stage in production without comparing matrix values, read classifications, and
resource usage on a representative dataset.

## Command-line usage

The normal run only needs an indexed, coordinate-sorted BAM, an output BAM, and
a GTF. Isoform indexes are built in memory automatically:

```bash
mfs-stitcher \
  --input input.bam \
  --output stitched.bam \
  --gtf genes.gtf \
  --matrix-out mex_matrix/ \
  --counts-tsv isoform_counts.tsv.gz \
  --molecules-tsv molecules.tsv.gz \
  --threads 20
```

The `--matrix-out`, `--counts-tsv`, and `--molecules-tsv` flags are optional. Use `--skip-iso` when transcript compatibility is not needed. The `--isoform` and `--junction` flags must be supplied together when reusing precomputed indexes. To build indexes without processing a BAM, use:

```bash
mfs-stitcher --gtf genes.gtf --index-only reference/prefix
```

Advanced filtering and BAM tag options are available through `--cells`,
`--genes`, `--cell-tag`, `--umi-tag`, `--contig`, and `--gene-identifier`.

## Input Parameters & CLI Options

The parameters are organized into four functional groups:

### 1. Core Inputs
| Option | Type | Required | Default | Description |
| :--- | :--- | :---: | :--- | :--- |
| `-i`, `--input`, `--input-bam` | `Path` | Yes* | - | Coordinate-sorted and indexed input BAM file (must have `.bai` in same directory). |
| `-g`, `--gtf` | `Path` | Yes | - | GTF annotation file (plain `.gtf` or gzip-compressed `.gtf.gz`). Used to define gene regions and exon/splice models. |

*\*Not required when running with `--index-only`.*

### 2. Outputs
| Option | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-o`, `--output`, `--output-bam` | `Path` | - | Target path for the stitched consensus BAM file. A `.bai` index will be automatically generated alongside it. |
| `--matrix-out` | `Path` | - | Directory path for exporting standard 10x-compatible MEX format sparse count matrices (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`). |
| `--counts-tsv` | `Path` | - | Path to save the aggregated cell-by-isoform quantification table (e.g. `isoform_counts.tsv.gz`). |
| `--molecules-tsv` | `Path` | - | Path to save the single-molecule audit ledger detailing every stitched molecule (e.g. `molecules.tsv.gz`). |

### 3. Annotation and Indexing
| Option | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--isoform` | `Path` | - | Path to a precomputed interval index JSON (`.intervals.json.gz`). Must be used together with `--junction`. |
| `--junction` | `Path` | - | Path to a precomputed junction skip index JSON (`.refskip.json.gz`). Must be used together with `--isoform`. |
| `--dump-index` | `String` | - | Write the in-memory GTF isoform indexes to `<PREFIX>.intervals.json.gz` and `<PREFIX>.refskip.json.gz` while processing. |
| `--index-only` | `String` | - | Build and save the GTF isoform index only without processing any BAM file. |
| `--skip-iso` | `Flag` | `false` | Skip transcript isoform calling. Skips building/querying transcript tries, producing stitched BAMs faster with lower memory. Cannot be used with `--matrix-out`. |
| `--contig` | `String` | - | Restrict execution or indexing to a single chromosome/contig (e.g. `--contig chr1`). Ideal for distributed cluster jobs. |
| `--gene-identifier`| `String` | `gene_id` | Attribute key to extract from GTF attribute column (field 9). Typical options: `gene_id` or `gene_name`. |

### 4. Read & Runtime Options
| Option | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-t`, `--threads` | `Int` | `4` | Number of worker threads for parallel chromosome and gene processing. |
| `--single-end` | `Flag` | `false` | Process reads as single-end. By default (paired-end), only properly paired and mapped mate pairs are used. |
| `--cell-tag` | `String` | `CB` | BAM tag containing the error-corrected cell barcode. |
| `--umi-tag`, `--UMI-tag` | `String` | `UB` | BAM tag containing the molecular identifier (UMI). |
| `--cells` | `Path` | - | Text file with cell whitelist barcodes (one per line). Reads with barcodes not in this list will be ignored. |
| `--genes` | `Path` | - | Text file containing target gene IDs/names (one per line). Restricts processing to only these genes. |


## Output Files & Specifications

`mfs-stitcher` produces four main categories of outputs:

### 1. Stitched Consensus BAM (`--output stitched.bam`) & Index (`.bai`)

Each BAM record represents an assembled full-length consensus transcript molecule for a unique `(Cell, Gene, UMI)` triplet. Output records are coordinate-sorted and indexed with an accompanying `.bai` file.

- **Read Name (`QNAME`)**: `<Cell>:<Gene>:<UMI>` (e.g. `CellA:ENSG00000123456:AGCTAGCT`)
- **Flags**: `0` (forward strand) or `16` (reverse strand)
- **MAPQ**: `255`
- **CIGAR**: Consensus alignment operations (`M` for aligned bases, `N` for spliced introns, `D` for deletions)
- **SAM Tags**:
  | Tag | Type | Description |
  | :--- | :--- | :--- |
  | `CB` / `BC` | `Z` | Cell barcode |
  | `UB` / `UM` | `Z` | Molecule UMI sequence (configured via `--umi-tag`) |
  | `GX` / `XT` | `Z` | Assigned gene ID |
  | `CT` | `Z` | Compatible transcript/isoform IDs (comma-separated, e.g. `ENST00000380152,ENST00000543210`) |
  | `NR` | `i` | Number of raw sequencing reads supporting this consensus molecule |
  | `ER` | `i` | Number of exonic reads contributing to this molecule |
  | `IR` | `i` | Number of intronic reads contributing to this molecule |
  | `NC` | `i` | Number of conflicting bases between spliced and continuous alignments (if conflict occurred) |
  | `IL` | `B,I` | Genomic interval coordinates where alignment conflicts were detected |

### 2. Isoform Counts Summary (`--counts-tsv isoform_counts.tsv.gz`)

Aggregated per-cell, per-isoform expression quantification table (gzipped TSV).

| Column | Description |
| :--- | :--- |
| `cell` | Cell barcode identifier |
| `gene_id` | Ensembl/custom gene ID |
| `transcript_id` | Ensembl/custom transcript (isoform) ID |
| `unique_counts` | Number of UMIs unambiguously assigned to only this transcript |
| `fractional_counts` | Fractional UMI counts allocating ambiguous multi-isoform molecules evenly ($1/N$) |
| `total_reads` | Total supporting sequencing reads for this cell and transcript |

### 3. Molecule-level Audit Detail (`--molecules-tsv molecules.tsv.gz`)

Comprehensive single-molecule level ledger detailing every stitched molecule (gzipped TSV).

| Column | Description |
| :--- | :--- |
| `cell` | Cell barcode identifier |
| `gene_id` | Gene ID |
| `umi` | Molecule UMI sequence |
| `transcripts` | Comma-separated list of compatible transcript IDs (or `Unassigned`) |
| `n_transcripts` | Number of compatible transcripts |
| `reads` | Number of raw sequencing reads supporting this molecule |
| `is_unique` | `1` if uniquely assigned to a single isoform; `0` if ambiguous ($N > 1$) or unassigned |

### 4. 10x-Compatible MEX Matrix (`--matrix-out <DIR>`)

A standard sparse MatrixMarket directory for seamless loading into Seurat (`Read10X`) or Scanpy (`sc.read_10x_mtx`):

- `matrix.mtx.gz`: Sparse count matrix (rows: transcripts/features, columns: cell barcodes)
- `features.tsv.gz`: Feature annotations (`<transcript_id>\t<transcript_id>\tGene Expression`)
- `barcodes.tsv.gz`: Cell barcodes list

