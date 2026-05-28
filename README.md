# Organelle Intergenic Region Analysis

Data and analysis scripts accompanying the manuscript on organelle intergenic regions across eukaryotic lineages.

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18614283-blue)](https://doi.org/10.5281/zenodo.18614283)

## Repository Structure

```
.
├── code/                           # Jupyter notebooks for figure generation and analysis
│   ├── figure1.ipynb               # Jupyter notebook with analysis needed for Figure 1
│   ├── figure2.ipynb               # Jupyter notebook with analysis needed for Figure 2
│   ├── figure3.ipynb               # Jupyter notebook with analysis needed for Figure 3
│   ├── table1.ipynb                # Jupyter notebook with analysis needed for Table 1
│   ├── QC_ANS_removed.tsv          # List of ANs removed after quality control
│   └── supplemental_material.ipynb # Jupyter notebook with analysis needed for Supplemental Material
│
├── scripts/                        # Processing scripts
│   ├── brms_2bin.R                 # R code for Bayesian analysis of intergenic region data (2-bin model)
│   ├── brms_3bin.R                 # R code for Bayesian analysis of intergenic region data (3-bin model, default)
│   ├── brms.py                     # Python wrapper for brms analysis (called by run_brms.sh)
│   ├── prepare_igr.sh              # Bash script to extract intergenic regions
│   ├── run_brms.sh                 # Bash script to run Bayesian analysis
│   └── summary_igs.py              # Python script to summarize intergenic region data (called by prepare_igr.sh)
│
├── *_compressed.tar.zst            # Compressed data archives by taxonomic group
└── all_groups_lean.tar.zst         # Lean archive (IGR summary data only)
```

## Data Groups

- **Mitochondria**: `fungi_mit`, `green_algae_mit`, `protists_mit`, `metazoans_mit`, `plants_mit`
- **Plastid**: `green_algae_plt`, `plants_plt`, `protists_plt`

## Requirements

The pipeline has three stages, each with its own dependencies. Install whatever applies to the steps you intend to run.
Python ([3.12 or higher](https://www.python.org/downloads)) is required for all stages, with [R](https://www.r-project.org/) (4.5.2) only needed for the Bayesian analysis stage.

**1. IGR extraction** (`prepare_igr.sh`):

Python packages:
```bash
pip install "tigre[all]"
```

**2. Bayesian analysis** (`run_brms.sh`):

Python packages:
```bash
pip install polars ete4
```
R packages:
```r
install.packages(c("brms", "ape", "posterior"))
```

**3. Figures and supplemental notebooks** (`code/*.ipynb`):

```bash
pip install polars matplotlib scipy seaborn
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

**2. Extract intergenic regions**

```bash
./scripts/prepare_igr.sh --verbose --overwrite 2>&1 | tee prepare_igr.log
```

This extracts intergenic regions using tigre and generates summary statistics.

**3. Run Bayesian analysis**

```bash
./scripts/run_brms.sh --verbose --overwrite 2>&1 | tee run_brms.log
```

This performs the brms statistical analysis on the extracted intergenic regions.

### Using the Lean Archive

The `all_groups_lean.tar.zst` archive contains pre-computed summary files, allowing you to skip steps 1-2:

```bash
# Decompress lean archive
tar -xf all_groups_lean.tar.zst

# Run Bayesian analysis only
./scripts/run_brms.sh --verbose --overwrite 2>&1 | tee run_brms.log
```

#### Nohup Usage

The brms analysis can take several hours. We do recommend running it with `nohup` to avoid interruptions:

```bash
nohup ./scripts/run_brms.sh --verbose --overwrite > run_brms_nohup.log 2>&1 &
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
### 2-bin Model

By default, the analysis runs with the 3-bin model . To run with the 2-bin model instead, add the `--2bin` flag to the `run_brms.sh` script:

```bash
./scripts/run_brms.sh --2bin --verbose --overwrite 2>&1 | tee run_brms_2bin.log
```
### NCBI Topology Tree

The `brms` pipeline uses the NCBI Topology Tree provided at `<group>/brms_*bin/tree.nwk`. The script can generate new trees if needed with the `--new-tree` flag in `run_brms.sh`. However, due the everchanging nature of NCBI Taxonomy, the new tree may need some tweaking around the taxon ids. If an AN needs a new taxon id, change it in the `<group>.tsv` file, column `ncbi_taxid`. The script will then use the new taxon id to generate the tree.

Contact us if you need help with this step!


### Figures Generation

The `code/` directory contains Jupyter notebooks for generating Figures 1, 2 and 3, as well as the Table 1 and the Supplemental Material.

Those analysis scripts and notebooks are designed to be run with at least the `all_groups_lean.tar.zst` decompressed and it's run on all groups, but you can modify them to focus on specific groups if desired, by changing the `groups` variable. Some tweaking may be needed.