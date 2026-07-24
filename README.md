# Organelle Intergenic Region Analysis

Data and analysis scripts accompanying the manuscript on organelle intergenic regions across eukaryotic lineages.

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18614283-blue)](https://doi.org/10.5281/zenodo.18614283)

## Repository Structure

```
.
├── code/                              # Jupyter notebooks for figure and table generation
│   ├── figure1.ipynb
│   ├── figure2.ipynb
│   ├── figure3.ipynb
│   ├── QC_ANS_removed.tsv              # List of ANs removed after quality control
│   ├── supplementary_figure1.ipynb
│   ├── supplementary_table1.ipynb
│   ├── supplementary_table2.ipynb
│   ├── supplementary_table3.ipynb
│   ├── supplementary_table4.ipynb
│   └── supplementary_table5.ipynb
│
├── scripts/                           # Processing scripts
│   ├── 01_prepare_igr.sh              # Bash script to extract intergenic regions (calls summary_igs.py)
│   ├── 02_prepare_brms.ipynb          # Notebook that prepares brms input (filtered_3bin.tsv + tree.nwk) per group
│   ├── 03_run_brms.sh                 # Bash script to run the Bayesian (brms) analyses
│   ├── brms_contrast.R                # R code for the divergent vs convergent posterior contrast
│   ├── brms_divergence.R              # R code for MCMC diagnostics (divergences, treedepth, Rhat, ESS)
│   ├── brms_polarity.R                # R code for the Polarity model
│   ├── brms_sigma_icc.R               # R code for the sigma + ICC companion analysis
│   ├── brms_type_length.R             # R code for the Type/Length models (with k-fold comparison)
│   ├── compile_all_txts.sh            # Bash script to gather all per-group brms result txts into one file (optional, after 03)
│   ├── summary_igs.py                 # Python script to summarize intergenic region data (called by 01_prepare_igr.sh)
│   └── taxonomy_tree.ipynb            # Notebook documenting how the NCBI topology trees were built (tree.nwk / all_tree.nwk)
│
├── all_tree.nwk                      # Phylogenetic tree spanning all taxonomic groups combined
│
├── <group>_compressed.tar.zst        # Full per-group archive (8 total), everything except logs and fitted models
├── all_groups_lean.tar.zst           # All groups, essentials only: no gff3, no fasta, no fitted models
└── all_groups_models.tar.zst         # Fitted brms models (.rds) for all groups
```

## Data Groups

- **Mitochondria**: `fungi_mit`, `green_algae_mit`, `protists_mit`, `metazoans_mit`, `plants_mit`
- **Plastid**: `green_algae_plt`, `plants_plt`, `protists_plt`

## Requirements

The pipeline has four stages, each with its own dependencies. Install whatever applies to the steps you intend to run.
Python ([3.12+](https://www.python.org/downloads)) is required for all stages, with [R](https://www.r-project.org/) (4.5.2) only needed for the Bayesian analysis stage.

**1. IGR extraction** (`01_prepare_igr.sh`):

Python packages:
```bash
pip install "tigre[all]"
```

**2. Bayesian analysis** (`02_prepare_brms.ipynb` + `03_run_brms.sh`):

Python packages:
```bash
pip install polars
```
R packages:
```r
install.packages(c("brms", "ape", "posterior", "loo", "callr", "matrixStats"))
```

**3. Rebuilding the NCBI trees** (`taxonomy_tree.ipynb`, not normally needed):

Only relevant if you want to reproduce how the shipped trees were generated. Use the provided trees for the analysis itself, see [NCBI Topology Tree](#ncbi-topology-tree).
```bash
pip install polars ete4 biopython python-dotenv
```
This notebook queries NCBI, so it expects `ENTREZ_EMAIL` and `ENTREZ_API_KEY` in a `.env` file at the project root.

**4. Figures and supplementary notebooks** (`code/*.ipynb`):

```bash
pip install polars matplotlib scipy seaborn
```

## Reproducing the Analysis

### Full Pipeline (from raw data)

**1. Decompress data archives**

Each taxonomic group has its own archive containing the whole group folder (gff3, fasta, tsv, gdict, summary and brms results; see [Data Archives](#data-archives)):
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
./scripts/01_prepare_igr.sh 2>&1 | tee prepare_igr.log
```

This extracts intergenic regions using `tigre` and generates summary statistics.

**3. Prepare brms input**

Run `scripts/02_prepare_brms.ipynb`, setting `brms_run_name` to `brms_polarity` and running it, then set it to `brms_type_length` and run it again. Each run builds a `<group>/<brms_run_name>/` folder containing `filtered_3bin.tsv` and a copy of `tree.nwk` for every group.

**4. Run Bayesian analysis**

```bash
./scripts/03_run_brms.sh 2>&1 | tee run_brms.log
```

Runtimes vary a lot between groups, from several hours to a few days for the larger ones (`plants_plt` and `metazoans_mit`). We do recommend running it with `nohup` to avoid interruptions:

```bash
nohup ./scripts/03_run_brms.sh > run_brms_nohup.log 2>&1 &
```

For each group, this runs five stages in order:

| Stage | Script | Reads from | Outputs |
|---|---|---|---|
| Polarity | `brms_polarity.R` | `brms_polarity/` | `brms_polarity.txt`, `brms_polarity_row.tsv` |
| Type_Length | `brms_type_length.R` | `brms_type_length/` | `brms_type_length.txt`, `brms_type_length_row.tsv` |
| Sigma_ICC | `brms_sigma_icc.R` | `brms_type_length/` | `brms_sigma_icc.txt`, `brms_sigma_icc_row.tsv` |
| Contrast | `brms_contrast.R` | `brms_type_length/` | `brms_contrast.txt`, `brms_contrast_row.tsv` |
| Divergence | `brms_divergence.R` | both folders | `brms_divergence_row.tsv` in each folder |

The Polarity and Type_Length stages fit the models and cache each fit as `brms_model_cache_*.rds` next to the results. Sigma_ICC and Contrast do not refit anything: they read those cached fits, so the corresponding stage must have run first. Sigma_ICC reports the sigma ratios and the taxon-level intraclass correlation coefficient (ICC), and Contrast reports the divergent-minus-convergent posterior contrast. Divergence is a diagnostics pass over every cached fit, reporting divergent transitions, treedepth hits, Rhat and ESS.

If any stage fails, the script keeps going and prints a summary of failed runs at the end, exiting with a non-zero status.

**5. Combine the result summaries (optional)**

Once `03_run_brms.sh` has finished, this optional step gathers every per-group human-readable summary into a single file for easier reading:

```bash
./scripts/compile_all_txts.sh
```

It writes `all_groups_brms_results.txt` at the project root, with one section per group covering the Polarity, Type_Length, Sigma_ICC and Contrast summaries. Any summary that has not been produced yet is reported as missing rather than causing a failure. The Divergence stage is not included, since it only writes a `_row.tsv`. Like the other scripts, it accepts group names as arguments and defaults to all groups.

### Data Archives

Three kinds of archive are provided, pick whichever matches what you want to do:

- **`<group>_compressed.tar.zst`** (one per taxonomic group): the complete group folder: gff3 and fasta files, `<group>.tsv`, `<group>.gdict`, `<group>_summary_igr.tsv`, `tree.nwk`, and the contents of `brms_polarity/` and `brms_type_length/`. Log files and fitted models (`.rds`) are excluded. **Use these to replicate the analysis in full**, from step 1.
- **`all_groups_lean.tar.zst`**: the same essentials for every group at once: `<group>.tsv`, `<group>.gdict`, `<group>_summary_igr.tsv`, `tree.nwk`, and everything in the two brms folders apart from the fitted models, plus the combined `all_tree.nwk`. **Use it to inspect the results**, or to rerun the analysis straight from step 4 (`03_run_brms.sh`). This is the minimum you need to run the figure and table generation notebooks in `code/`.
- **`all_groups_models.tar.zst`**: the fitted brms models (`.rds`) only, which none of the other archives contain. **Use it to inspect the models.**

`all_groups_models.tar.zst` is only available in [Zenodo](https://doi.org/10.5281/zenodo.18614283).

### Processing Specific Groups

Both scripts accept group names as arguments. If no arguments are provided, all groups are processed.

**Process all groups (default):**
```bash
./scripts/01_prepare_igr.sh
./scripts/03_run_brms.sh
```

**Process specific groups:**
```bash
# Process only fungi and plants
./scripts/01_prepare_igr.sh fungi_mit plants_mit plants_plt
./scripts/03_run_brms.sh fungi_mit plants_mit plants_plt
```

### NCBI Topology Tree

The `brms` pipeline uses the NCBI Topology Tree provided at `<group>/tree.nwk` (copied into `<group>/<brms_run_name>/tree.nwk` by `02_prepare_brms.ipynb`, then pruned to the taxa present in the filtered data).

**We recommend using the trees provided with the data.** They were built on **July 7th 2026, 14h**, and every model in the manuscript was fitted against them. `scripts/taxonomy_tree.ipynb` is included to document how they were generated, not as a step you need to run: it queries NCBI Taxonomy through `ete4` and Entrez, validates that every taxon id in `<group>.tsv` is present in the resulting topology, and writes the per-group `tree.nwk` files along with the combined `all_tree.nwk`.

Because NCBI Taxonomy changes continuously, rebuilding the trees against a newer release can produce topologies that no longer match the taxon ids in the shipped `<group>.tsv` files, which will cause incompatibilities downstream. If you do rebuild them and an AN needs a different taxon id, change it in the `<group>.tsv` file, column `ncbi_taxid`, before re-running `02_prepare_brms.ipynb`. Alternatively, you can drop the AN from the `<group>.tsv` file entirely, which removes it from the analysis.

Note that any change to a group's `<group>.tsv` requires re-running the pipeline from the start (`01_prepare_igr.sh`), because `<group>_summary_igr.tsv` is derived from it and has to be rebuilt before `02_prepare_brms.ipynb` can pick the change up.

Contact us if you need help with this step!


### Figures and Tables Generation

The `code/` directory contains the notebooks that produce Figures 1-3, Supplementary Figure 1 and Supplementary Tables 1-5. Each writes its output next to itself in `code/` (`figure1.png`, `supplementary_table1.tsv`, and so on).