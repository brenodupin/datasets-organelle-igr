# Organelle Intergenic Region Analysis

Data and analysis scripts accompanying the manuscript on organelle intergenic regions across eukaryotic lineages.

[![DOI](https://zenodo.org/badge/1154845769.svg)](https://zenodo.org/badge/latestdoi/1154845769)

## Repository Structure

```
.
├── code/                          # Jupyter notebooks for figure generation and analysis
│   ├── table1.ipynb               # Jupyter notebook with analysis needed for Table 1
│   ├── table1.tsv                 # Table 1 output TSV file
│   ├── QC_ANS_removed.tsv         # List of ANs removed after quality control
│   └── figure1.ipynb              # Jupyter notebook for generating Figure 1
│
├── scripts/                       # Processing scripts
│   ├── prepare_igr.sh             # Extract intergenic regions
│   ├── run_brms.sh                # Run Bayesian analysis
│   ├── brms.py                    # Python wrapper for brms analysis (called by run_brms.sh)
│   ├── summary_igs.py             # Python script to summarize intergenic region data (called by prepare_igr.sh)
│   └── create_brms.R              # R script to fit brms models (called by brms.py)
│
├── *_compressed.tar.zst           # Compressed data archives by taxonomic group
└── all_groups_lean.tar.zst        # Lean archive (IGR summary data only)
```

## Data Groups

- **Mitochondria**: fungi_mit, green_algae_mit, protists_mit, metazoans_mit, plants_mit
- **Plastid**: green_algae_plt, plants_plt, protists_plt

## Requirements

```bash
# Install tigre
pip install "tigre[all]"

# Python dependencies
pip install pandas numpy matplotlib polars ete4
```

```r
# R packages
install.packages(c("brms", "ape"))
```

## Reproducing the Analysis

### Full Pipeline (from raw data)

**1. Decompress data archives**

To decompress specific groups:
```bash
tar -xf fungi_mit_compressed.tar.zst
tar -xf green_algae_mit_compressed.tar.zst
```

To decompress all groups:
```bash
for file in *_compressed.tar.zst; do tar -xf "$file"; done
```

For older versions of `tar` (<1.31), use:
```bash
# Decompress all_groups_lean.tar.zst
zstd -dc all_groups_lean.tar.zst | tar -xf -

# Decompress all groups with zstd support
for file in *_compressed.tar.zst; do zstd -dc "$file" | tar -xf -; done
```

**2. Extract intergenic regions**

```bash
./scripts/prepare_igr.sh 2>&1 | tee prepare_igr.log
```

This extracts intergenic regions using tigre and generates summary statistics.

**3. Run Bayesian analysis**

```bash
./scripts/run_brms.sh 2>&1 | tee run_brms.log
```

This performs the brms statistical analysis on the extracted intergenic regions.

### Using the Lean Archive

The `all_groups_lean.tar.zst` archive contains pre-computed summary files, allowing you to skip steps 1-2:

```bash
# Decompress lean archive
tar -xf all_groups_lean.tar.zst

# Run Bayesian analysis only
./scripts/run_brms.sh 2>&1 | tee run_brms.log
```

#### Nohup Usage

The brms analysis can take several hours. We do recommend running it with `nohup` to avoid interruptions:

```bash
nohup ./scripts/run_brms.sh > run_brms_nohup.log 2>&1 &
```

### Processing Specific Groups

Both scripts accept group names as arguments. If no arguments are provided, all groups are processed.

**Process all groups (default):**
```bash
./scripts/prepare_igr.sh
./scripts/run_brms.sh
```

**Process specific groups:**
```bash
# Process only fungi and plants
./scripts/prepare_igr.sh fungi_mit plants_mit plants_plt
./scripts/run_brms.sh fungi_mit plants_mit plants_plt
```

**Available groups:**
- Mitochondrial: `fungi_mit`, `green_algae_mit`, `protists_mit`, `metazoans_mit`, `plants_mit`
- Plastid: `green_algae_plt`, `plants_plt`, `protists_plt`

### Figure 1 and Table 1 Generation

The `code/` directory contains Jupyter notebooks for generating Figure 1 and Table 1.

Those analysis scripts and notebooks are designed to be run with all the data groups available, but you can modify them to focus on specific groups if desired, by changing the data loading steps to filter for the groups of interest.