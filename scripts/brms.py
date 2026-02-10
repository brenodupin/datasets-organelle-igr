#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pipeline to fit BRMS model to IGS data."""

import subprocess
import sys
import traceback
from pathlib import Path
import argparse
import logging

import polars as pl
from ete4 import NCBITaxa, Tree  # type: ignore[import-not-found]

R_CMD = (
    "Rscript {script_path} {tree_path} {igs_summary} {out_name} {rds_name} {table_name}"
)


def get_first_descendants(
    ncbi_taxa: NCBITaxa,
    taxids: set[int],
) -> dict[int, int]:
    """Get the first descendant for each taxid in the set.

    Args:
        ncbi_taxa: NCBITaxa instance.
        taxids: Set of taxids to find descendants for.

    Returns:
        Dictionary mapping original taxid to first descendant taxid.

    """
    query = f"""
    SELECT parent, MIN(taxid) as first_descendant
    FROM species
    WHERE parent IN ({','.join(map(str, taxids))})
    GROUP BY parent
    """
    return dict(ncbi_taxa.db.execute(query).fetchall())


def validade_tree(
    tree: Tree,
    taxa: set[str],
    log: logging.Logger,
) -> tuple[set[str], set[str]]:
    """Validate that the tree contains exactly the given taxa.

    Args:
        tree: ETE4 Tree instance.
        taxa: Set of taxids as strings.
        log: Logger instance.

    Returns:
        A tuple containing two sets:
            - Missing taxa in the tree.
            - Extra taxa in the tree.

    """
    tree_taxa = set(leaf.name for leaf in tree.leaves())
    missing_in_tree = taxa - tree_taxa
    extra_in_tree = tree_taxa - taxa
    if missing_in_tree:
        log.warning(f"Taxa missing in tree: {', '.join(missing_in_tree)}")
    if extra_in_tree:
        log.error(f"Extra taxa in tree: {', '.join(extra_in_tree)}")
    return missing_in_tree, extra_in_tree


def run_single(
    tsv_in: Path,
    igs_in: Path,
    out_folder: Path,
    an_column: str,
    log: logging.Logger,
    update_taxa: bool = False,
) -> None:
    """Run the BRMS model fitting.

    Args:
        tsv_in: Path to TSV file with required columns.
        igs_in: Path to IGS summary TSV file.
        out_folder: Output folder path.
        an_column: Name of the accession number column in IGS file.
        log: Logger instance.
        remove_not_in_tree: If True, remove entries not in the tree instead of failing.
        update_taxa: If True, update the NCBI taxonomy database before processing.

    """
    tsv = pl.read_csv(tsv_in, separator="\t")
    required_columns = {"ncbi_taxid", "Species", "ncbi_name", "Genome_length"}
    if not required_columns.issubset(tsv.columns):
        log.error(
            f"TSV file {tsv_in} is missing required columns. "
            "Make sure to run iga populate before this."
        )
        sys.exit(1)
    if "GI" in tsv.columns:
        tsv = tsv.drop("GI")

    ncbi = NCBITaxa()
    if update_taxa:
        log.info("Updating NCBI taxonomy database...")
        ncbi.update_taxonomy_database()

    out_folder = Path(out_folder).resolve()
    out_folder.mkdir(parents=True, exist_ok=True)
    tree_path = out_folder / "tree.nwk"

    log.info("Creating ultrametric tree for BRMS...")
    tree = ncbi.get_topology(tsv["ncbi_taxid"])
    tree.to_ultrametric()

    log.info("Tree created, now validating it...")
    taxa_in_tsv = set(tsv["ncbi_taxid"].cast(pl.String))
    missing_in_tree, extra_in_tree = validade_tree(tree, taxa_in_tsv, log)

    if not missing_in_tree and not extra_in_tree:
        log.info("Tree validation successful: all taxa match.")
        tsv = tsv.with_columns(pl.col("ncbi_taxid").alias("taxon_tree"))

    # first descendant
    if missing_in_tree:
        log.info("Replacing missing taxa with first descendant...")

        missing_int: set[int] = {int(x) for x in missing_in_tree}

        descendants = get_first_descendants(ncbi, missing_int)

        tsv = tsv.with_columns(
            pl.col("ncbi_taxid")
            .replace_strict(descendants, default=pl.col("ncbi_taxid"))
            .alias("taxon_tree")
        )

        tree = ncbi.get_topology(tsv["taxon_tree"])
        tree.to_ultrametric()

        log.info("Re-validation of the tree after replacement...")
        taxa_in_tsv2 = set(tsv["taxon_tree"].cast(pl.String))
        missing_in_tree2, _ = validade_tree(tree, taxa_in_tsv2, log)

        # first descendant from first descendant
        if missing_in_tree2:
            log.warning(
                "Some taxa are still missing in the tree after replacement, going "
                "for second round."
            )
            log.debug(f"Missing taxa: {', '.join(missing_in_tree2)}")

            missing_int2: set[int] = {int(x) for x in missing_in_tree2}

            descendants2 = get_first_descendants(ncbi, missing_int2)

            tsv = tsv.with_columns(
                pl.col("taxon_tree")
                .replace_strict(descendants2, default=pl.col("taxon_tree"))
                .alias("taxon_tree")
            )

            tree = ncbi.get_topology(tsv["taxon_tree"])
            tree.to_ultrametric()

            log.info("Re-re-validation of the tree after replacement...")
            taxa_in_tsv3 = set(tsv["taxon_tree"].cast(pl.String))
            missing_in_tree3, _ = validade_tree(tree, taxa_in_tsv3, log)

            if missing_in_tree3:
                log.error(
                    "Some taxa are still missing in the tree after second replacement, "
                    "we give up."
                )
                log.debug(f"Missing taxa: {', '.join(missing_in_tree3)}")
                sys.exit(1)
            else:
                log.info(
                    "All missing taxa replaced successfully by their first descendant, "
                    "first descendant."
                )

        else:
            log.info(
                "All missing taxa replaced successfully by their first descendant."
            )

    tree.write(outfile=tree_path, parser=1)

    log.info("Preparing IGS data for BRMS fitting...")
    igs = pl.scan_csv(igs_in, separator="\t")
    igs = igs.join(tsv.lazy(), on=an_column, how="left")

    igs = igs.with_columns(
        [
            pl.col("Polarity")
            .replace({"++": "same", "--": "same", "+-": "opposite", "-+": "opposite"})
            .alias("polarity_bin"),
            pl.col("Length").log10().alias("log10_length"),
        ]
    )
    igs_path = out_folder / "filtered.tsv"
    igs = igs.sort([an_column, "ID"], descending=[False, False])
    igs.sink_csv(igs_path, separator="\t")

    log_path = out_folder / "brms_fit.log"
    r_script = (Path(__file__).parent / "create_brms.R").resolve()

    cmd_str = R_CMD.format(
        script_path=r_script,
        tree_path=tree_path,
        igs_summary=igs_path,
        out_name=out_folder / "brms_result.txt",
        rds_name=out_folder / "brms_model.rds",
        table_name=out_folder / "brms_results_row.tsv",
    )
    log.debug(f"Running command: {cmd_str}")
    try:
        log.info(f"Fitting BRMS model, log at {log_path}")
        with open(log_path, "w") as log_file:
            subprocess.run(
                cmd_str.split(),
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=True,  # Raises CalledProcessError if returncode != 0
            )
        log.info("BRMS model fitting completed successfully.")

    except subprocess.CalledProcessError as e:
        err_file = out_folder / "brms_err.log"
        log.error(f"Error {e.returncode} in brms fit, check {log_path} for details.")
        err_file.write_text(traceback.format_exc())
        sys.exit(1)


def simple_log_setup(name: str, verbose: bool) -> logging.Logger:
    """Setup a simple logger."""
    log = logging.getLogger(name)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    log.propagate = False

    if not log.handlers:
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG if verbose else logging.INFO)
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        ch.setFormatter(formatter)
        log.addHandler(ch)

    return log


def main() -> None:
    """Main function to run BRMS fitting."""
    parser = argparse.ArgumentParser(description="Fit BRMS model to IGS data.")
    parser.add_argument(
        "--tsv",
        type=Path,
        required=True,
        help="Path to TSV file with required columns.",
    )
    parser.add_argument(
        "--igs",
        type=Path,
        required=True,
        help="Path to IGS summary TSV file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=False,
        default=Path.cwd() / "brms_output",
        help="Path to output folder. Defaults to 'brms_output' in the current directory.",
    )
    parser.add_argument(
        "--an-column",
        type=str,
        required=False,
        default="AN",
        help="Name of the accession number column in TSV/IGS file."
        " Defaults to 'AN'.",
    )
    parser.add_argument(
        "--update-taxa",
        action="store_true",
        required=False,
        help="Force update the NCBI taxonomy database before processing.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        required=False,
        help="Enable verbose logging output.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        required=False,
        help="Overwrite existing output folder if it exists.",
    )
    args = parser.parse_args()

    # Setup simple logger, only to console
    log = simple_log_setup("brms", args.verbose)

    # check input files
    args.tsv = Path(args.tsv).resolve()
    args.igs = Path(args.igs).resolve()
    args.output = Path(args.output).resolve()

    if not args.tsv.is_file():
        log.error(f"TSV input file {args.tsv} does not exist or is not a file.")
        sys.exit(1)

    if not args.igs.is_file():
        log.error(f"IGS input file {args.igs} does not exist or is not a file.")
        sys.exit(1)

    if args.output.exists():
        if args.overwrite:
            log.warning(
                f"Output folder {args.output} exists, overwriting as requested."
            )
        else:
            log.error(
                f"Output folder {args.output} already exists. "
                "Use --overwrite to overwrite."
            )
            sys.exit(1)

    run_single(
        tsv_in=args.tsv,
        igs_in=args.igs,
        out_folder=args.output,
        an_column=args.an_column,
        log=log,
        update_taxa=args.update_taxa,
    )


if __name__ == "__main__":
    main()
