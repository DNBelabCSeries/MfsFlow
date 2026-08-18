# Troubleshooting Guide

This guide helps you diagnose and resolve common issues with the MfsFlow pipeline.

## Table of Contents

1. [General Troubleshooting Steps](#general-troubleshooting-steps)
2. [Installation Issues](#installation-issues)
3. [Configuration Issues](#configuration-issues)
4. [Runtime Errors](#runtime-errors)
5. [Output Issues](#output-issues)
6. [Performance Issues](#performance-issues)
7. [Report Issues](#report-issues)
8. [Getting Help](#getting-help)

## General Troubleshooting Steps

Before diving into specific issues, follow these general steps:

### 1. Check Log Files
```bash
# View main pipeline log
less XPRESS_PROCESSING/logs/pipeline.log

# Check for errors
grep -i "error\|fail\|exception" XPRESS_PROCESSING/logs/pipeline.log
```

### 2. Verify Input Files
```bash
# Check FASTQ files exist and are readable
ls -lh /path/to/fastqs/
zcat /path/to/fastqs/Sample_R1.fq.gz | head -4

# Check reference files exist
ls -lh /path/to/reference/star/
ls -lh /path/to/reference/genes/genes.gtf
```

### 3. Run with Fewer Threads
Sometimes reducing thread count helps isolate issues:
```bash
mfsflow ... --threads 4
```

### 4. Check System Resources
```bash
# Check available RAM
free -h

# Check disk space
df -h

# Check CPU load
top
```

### 5. Enable Verbose Output
Check if scripts support verbose output (add `set -x` to bash scripts).

## Installation Issues

### Issue: Bundled tools not found

**Error message**:
```
Error: STAR command not found
```
or
```
Error: Cannot find bundled STAR in software/ directory
```

**Cause**: The bundled tools in `software/` directory are missing or incomplete.

**Solution**:
1. **Check if software/ directory exists**:
   ```bash
   ls -lh /path/to/MfsFlow/software/
   # Should contain: STAR, featureCounts, samtools, pigz, seqkit
   ```

2. **Re-download the toolkit**:
   ```bash
   git clone https://github.com/lishuangshuang0616/MfsFlow.git
   # or update existing clone
   cd MfsFlow
   git pull origin main
   ```

3. **Check permissions**:
   ```bash
   chmod +x /path/to/MfsFlow/software/STAR \
     /path/to/MfsFlow/software/featureCounts \
     /path/to/MfsFlow/software/samtools \
     /path/to/MfsFlow/software/pigz \
     /path/to/MfsFlow/software/seqkit
   ```

   Recent versions repair a missing execute bit automatically. If the Python
   package directory is read-only, MfsFlow copies the affected tool to
   `~/.cache/mfsflow/tools/` and runs that copy. Set `MFSFLOW_TOOL_CACHE` to a
   writable directory when the default home cache is unavailable.

4. **Use system-installed tools** (fallback):
   Install tools via Conda or system package manager (see [Installation Guide](installation.md))

### Issue: Python dependencies installation fails

**Error message**:
```
pip install -r requirements.txt
...
ERROR: Failed building wheel for pysam
```

**Solution**:
```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get install zlib1g-dev libcurl4-openssl-dev libssl-dev

# Retry installation
pip install -r requirements.txt

# If still failing, try conda
conda install -c bioconda pysam
```

### Issue: pysam installation fails on macOS

**Error message**:
```
fatal error: 'htslib/hts.h' file not found
```

**Solution**:
```bash
# Install htslib via homebrew
brew install htslib

# Set environment variables
export CFLAGS="-I$(brew --prefix htslib)/include"
export LDFLAGS="-L$(brew --prefix htslib)/lib"

# Retry installation
pip install pysam
```

### Issue: MfsFlow command not found after installation

**Error message**:
```
mfsflow: command not found
```

**Solution**:
```bash
# Check if installed
pip show mfsflow

# If installed but not in PATH, reinstall with pip
pip install -e .

# Or run directly
python3 run_analysis_pipeline.py ...
```

## Configuration Issues

### Issue: Barcode detection fails

**Error message**:
```
ERROR: None of the expected barcodes were observed in the data.
```

**Cause**: Strict barcode mode (`--plate`, `--manual`, or `--expectBarcode`) but barcodes don't match.

**Solution**:
1. **Check barcode file**:
   ```bash
   # View expected barcodes
   cat XPRESS_PROCESSING/config/expect_id_barcode.tsv
   ```

2. **Check observed barcodes**:
   ```bash
   # Extract barcodes from FASTQ
   zcat /path/to/fastqs/Sample_R2.fq.gz | \
       awk 'NR%4==2 {print substr($0, length($0)-19)}' | \
       sort | uniq -c | sort -rn | head -20
   ```

3. **Use discovery mode** (just omit `--plate`/`--manual`/`--expectBarcode`):
   ```bash
   mfsflow --fastqs ... --genomeDir ... --sample ...
   ```

### Issue: Reference genome not found

**Error message**:
```
ERROR: STAR index not found at /path/to/reference/star
```

**Solution**:
1. **Check directory structure**:
   ```bash
   ls -lh /path/to/reference/
   # Should contain:
   # - star/ (STAR index directory)
   # - genes/genes.gtf (GTF file)
   ```

2. **Build STAR index** (see [Installation Guide](installation.md#installing-reference-genome))

3. **Specify correct path**:
   ```bash
   mfsflow ... --genomeDir /correct/path/to/reference
   ```

### Issue: Input FASTQ files not found

**Error message**:
```
ERROR: No FASTQ files found in /path/to/fastqs
```

**Solution**:
1. **Check directory exists**:
   ```bash
   ls -lh /path/to/fastqs/
   ```

2. **Check file naming**:
   - Files should end with `.fq.gz` or `.fastq.gz`
   - Should have both R1 and R2 files

3. **Specify correct path**:
   ```bash
   mfsflow ... --fastqs /correct/path/to/fastqs
   ```

## Runtime Errors

### Issue: STAR alignment fails

**Error message**:
```
EXITING because of FATAL ERROR: unknown parameter in input: ...
```

**Solution**:
1. **Check STAR version**:
   ```bash
   STAR --version
   # Should be 2.7.10b or later
   ```

2. **Check additional STAR params**:
   - Remove unsupported parameters from config

3. **Check reference index**:
   - Rebuild STAR index with correct parameters

### Issue: Out of memory (OOM) error

**Error message**:
```
Killed
```
or
```
ERROR: Cannot allocate memory
```

**Solution**:
1. **Reduce thread count**:
   ```bash
   mfsflow ... --threads 10  # Reduce from 20
   ```

2. **Use streaming mode** (already default):
   ```yaml
   performance_opts:
     stream_bc_correction: true
   ```

### Issue: FeatureCounts fails

**Error message**:
```
ERROR: failed to open the annotation file
```

**Solution**:
1. **Check GTF file exists**:
   ```bash
   ls -lh /path/to/reference/genes/genes.gtf
   ```

2. **Check GTF file format**:
   ```bash
   head -20 /path/to/reference/genes/genes.gtf
   # Should have 9 columns, tab-separated
   ```

3. **Rebuild reference** (if GTF is corrupted)

### Issue: Pipeline hangs or freezes

**Possible causes**:
- Insufficient disk space
- Network filesystem latency
- Deadlock in parallel processing

**Solution**:
1. **Check disk space**:
   ```bash
   df -h /path/to/output
   ```

2. **Check if process is running**:
   ```bash
   ps aux | grep mfsflow
   ```

3. **Check for errors in log**:
   ```bash
   tail -f XPRESS_PROCESSING/logs/pipeline.log
   ```

4. **Restart from last successful stage**:
   ```bash
   mfsflow ... --stage Counting  # Resume from Counting
   ```

## Output Issues

### Issue: Output directory already exists

**Error message**:
```
ERROR: Output directory already exists: /path/to/output
```

**Solution**:
1. **Remove existing directory** (if safe):
   ```bash
   rm -rf /path/to/output
   ```

2. **Specify different output directory**:
   ```bash
   mfsflow ... --outdir /path/to/output_v2
   ```

3. **Use timestamp in output directory**:
   ```bash
   mfsflow ... --outdir /path/to/output_$(date +%Y%m%d_%H%M%S)
   ```

### Issue: Missing output files

**Symptom**: Pipeline completes but some expected files are missing.

**Solution**:
1. **Check pipeline log for errors**:
   ```bash
   grep -i "error\|fail" XPRESS_PROCESSING/logs/pipeline.log
   ```

2. **Check which stage failed**:
   ```bash
   # Check timing file
   cat XPRESS_PROCESSING/logs/pipeline_timing.tsv
   ```

3. **Resume from failed stage**:
   ```bash
   mfsflow ... --stage <failed_stage>
   ```

### Issue: BAM file is corrupted

**Error message**:
```
[E::bgzf_read] Read block operation failed
```

**Solution**:
1. **Re-run pipeline** (intermediate files might be corrupted)

2. **Check disk integrity**:
   ```bash
   fsck /dev/sdX  # (Linux) or Disk Utility (macOS)
   ```

3. **Use different storage** (e.g., local SSD instead of network filesystem)

## Performance Issues

### Issue: Pipeline runs too slowly

**Possible causes**:
- Insufficient threads
- Slow storage (HDD instead of SSD)
- Network filesystem latency

**Solution**:
1. **Increase thread count**:
   ```bash
   mfsflow ... --threads 40
   ```

2. **Use fast storage**:
   - Copy input data to local SSD
   - Use `/dev/shm` for temporary files

3. **Disable unnecessary outputs**:
   ```yaml
   make_ub_bam: false  # Don't generate UB BAM
   make_h5ad: false     # Don't generate H5AD (generate later if needed)
   ```

### Issue: High memory usage

**Solution**:
1. **Reduce thread count** (see above)

2. **Use streaming mode**:
   ```yaml
   performance_opts:
     stream_bc_correction: true
   ```

3. **Process fewer cells**:
   - Use `barcodes.barcode_num` to limit number of cells

4. **Increase swap space** (temporary solution):
   ```bash
   sudo fallocate -l 64G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile
   ```

## Report Issues

### Issue: HTML report not generated

**Symptom**: Pipeline completes but no HTML report in `outs/`.

**Solution**:
1. **Check if report generation is enabled**:
   ```yaml
   make_stats: yes  # Should be "yes"
   ```

2. **Check for errors in log**:
   ```bash
   grep -i "report" XPRESS_PROCESSING/logs/pipeline.log
   ```

3. **Generate report manually**:
   ```bash
   python3 -c "import sys; sys.path.insert(0, 'src'); \
       import report; \
       report.generate_multi_report('SAMPLE', '/path/to/output', config)"
   ```

### Issue: Report shows empty values

**Symptom**: Report generates but tables/plots are empty or show zeros.

**Possible causes**:
- Pipeline didn't complete successfully
- Stats files are missing or corrupted
- Incorrect data in stats files

**Solution**:
1. **Check stats files**:
   ```bash
   ls -lh XPRESS_PROCESSING/stats/
   head -5 XPRESS_PROCESSING/stats/qc_stats.tsv
   ```

2. **Re-run stats generation**:
   - Delete stats files
   - Resume from Summarising stage

3. **Check input data quality**:
   - Low-quality data might produce empty stats

### Issue: Plots not displaying in report

**Symptom**: Report HTML loads but plots are missing or broken.

**Solution**:
1. **Check Plotly JS is loaded**:
   - Open browser developer tools (F12)
   - Check for JavaScript errors

2. **Try different browser**:
   - Chrome, Firefox, Safari

3. **Check if plots are generated**:
   ```bash
   ls -lh XPRESS_PROCESSING/stats/*.png
   ```

## Getting Help

If you've tried the above solutions and still have issues:

### 1. Prepare Bug Report

Include the following information:
- MfsFlow version (`mfsflow --version` or check `pyproject.toml`)
- Python version (`python3 --version`)
- Operating system (Linux/macOS, version)
- Complete error message
- Relevant log file excerpts
- Command used to run pipeline
- Sample size (number of cells, reads)

### 2. Open GitHub Issue

Go to: https://github.com/lishuangshuang0616/MfsFlow/issues

Include:
- Descriptive title
- Steps to reproduce
- Expected vs actual behavior
- All information from step 1

### 3. Contact Development Team

For urgent issues or private inquiries, contact the development team (contact information to be added).

## FAQ

### Q: Can I run MfsFlow on Windows?
**A**: Not directly. Use WSL2 (Windows Subsystem for Linux) or Docker.

### Q: How much disk space do I need?
**A**: Approximately 10-20x the size of your input FASTQ files (for intermediate files, BAM files, etc.).

### Q: Can I process multiple samples simultaneously?
**A**: Yes, but ensure sufficient system resources (RAM, CPU, disk). Run each sample in a separate output directory.

### Q: How do I update MfsFlow to the latest version?
**A**:
```bash
cd /path/to/MfsFlow
git pull origin main
pip install -e .
```

### Q: Can I use MfsFlow for single-end data?
**A**: The pipeline is designed for paired-end data. Single-end support is limited.

### Q: How do I cite MfsFlow?
**A**: (Citation information to be added once published)
