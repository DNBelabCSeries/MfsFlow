# Installation Guide

## System Requirements

**Minimum**: Linux (Ubuntu 18.04+/CentOS 7+) or macOS (10.14+), 4 cores, 16GB RAM, 50GB storage, Python 3.8+

**Recommended**: 20+ cores, 64+ GB RAM, 200+ GB SSD, Python 3.9+

---

## Quick Installation

```bash
# 1. Clone repository
git clone https://github.com/lishuangshuang0616/MfsFlow.git
cd MfsFlow

# 2. Install dependencies
pip install -r requirements.txt

# 3. Install MfsFlow
pip install .

# 4. Verify
mfsflow --help
```

**Note**: MfsFlow comes with bundled bioinformatics tools (STAR, Samtools, Pigz, Seqkit, FeatureCounts) in `mfsflow/software/` (inside the installed package). **No manual installation required**.

---

## Detailed Steps

### 1. Install Python Dependencies

Using virtual environment (recommended):

```bash
# Create virtual environment
python3 -m venv mfsflow_env
source mfsflow_env/bin/activate

# Install dependencies
pip install -r requirements.txt
```

**Dependencies**:
- PyYAML >=5.1 (config parsing)
- pysam >=0.21 (BAM manipulation)
- numpy >=1.20, pandas >=1.3, scipy >=1.7 (computation)
- matplotlib >=3.5 (plotting)
- anndata >=0.8, h5py >=3.1 (HDF5 support)

### 2. Install MfsFlow

**Development mode** (recommended, no reinstall needed after code changes):

```bash
pip install -e .
```

**Direct installation**:

```bash
pip install .
```

---

## Verify Installation

```bash
# Check CLI
mfsflow --help

# Run tests
python3 -m unittest discover -s tests
```

---

## Troubleshooting

### pysam installation fails

**Ubuntu/Debian**:
```bash
sudo apt-get install zlib1g-dev libcurl4-openssl-dev libssl-dev
pip install pysam
```

**macOS**:
```bash
brew install htslib
export CFLAGS="-I$(brew --prefix htslib)/include"
export LDFLAGS="-L$(brew --prefix htslib)/lib"
pip install pysam
```

### `mfsflow` command not found

```bash
pip show mfsflow  # Check if installed
pip install -e .   # Reinstall
```

### Update MfsFlow

```bash
cd /path/to/MfsFlow
git pull origin main
pip install -e .
```

---

## Next Steps

1. Read [User Guide](user_guide.md)
2. Prepare reference genome
3. Run test analysis

## Getting Help

- Check [Troubleshooting Guide](troubleshooting.md)
- Open an issue on GitHub: https://github.com/lishuangshuang0616/MfsFlow/issues
