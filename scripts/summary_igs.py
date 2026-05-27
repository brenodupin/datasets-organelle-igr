#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compile Summaries of intergenic runs."""

import argparse
import concurrent.futures as cf
import os
from pathlib import Path
import tigre

os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl

_RS_strand = r"\|([+-])\|"

MAX_CPU = os.cpu_count() or 1

output_cols = [
    "ID",
    "AN",
    "UP",
    "DOWN",
    "Pair",
    "Polarity",
    "Length",
    "Merged",
    "source_up",
    "source_dw",
]


def summary_igs(
    an: str,
    gff_path: Path,
) -> pl.DataFrame:
    """Compile Summary of IGS outputs."""
    df = tigre.load_gff3(
        gff_path,
        usecols=("type", "start", "end", "attributes"),
        return_polars=True,
    )

    df = df.with_columns(
        pl.col("attributes").str.extract(tigre.gff3_utils._RS_ID, 1).alias("ID"),
        (pl.col("end") - pl.col("start") + 1).alias("Length"),
        pl.lit(an).alias("AN"),
        (pl.col("type") != "intergenic_region").cast(pl.Int8).alias("Merged"),
        pl.col("attributes")
        .str.extract(tigre.clean._RS_name_up, 1)
        .alias("UP")
        .fill_null("RS"),
        pl.col("attributes")
        .str.extract(tigre.clean._RS_name_dw, 1)
        .alias("DOWN")
        .fill_null("RE"),
        pl.col("attributes")
        .str.extract(tigre.clean._RS_source_up, 1)
        .alias("source_up"),
        pl.col("attributes")
        .str.extract(tigre.clean._RS_source_dw, 1)
        .alias("source_dw"),
    )

    df = df.with_columns(
        (pl.col("UP") + "-" + pl.col("DOWN")).alias("Pair"),
        pl.col("source_up")
        .str.extract(_RS_strand, 1)
        .fill_null("+")
        .alias("up-strand"),
        pl.col("source_dw")
        .str.extract(_RS_strand, 1)
        .fill_null("+")
        .alias("dw-strand"),
    )

    df = df.with_columns((pl.col("up-strand") + pl.col("dw-strand")).alias("Polarity"))
    df = df.select(output_cols)
    return df


def igs_multiple(
    log: tigre.GDTLogger,
    tsv_path: Path,
    output_file: Path,
    workers: int,
    gff_in_ext: str,
    gff_in_suffix: str,
    an_column: str = "AN",
) -> None:
    """Compile Summaries of tool outputs."""
    tigre.gff3_utils._ensure_spawn(log)
    tsv = pl.read_csv(tsv_path, separator="\t")

    gff_in_builder = tigre.PathBuilder(gff_in_ext).use_folder_builder(
        tsv_path.parent, gff_in_suffix
    )

    tigre.check_files(
        log,
        tsv,
        gff_in_builder,
        an_column,
        should_exist=True,
    )

    log.info(f"Starting processing {tsv.shape[0]} ANs with {workers} workers...")
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        tasks = [
            executor.submit(
                summary_igs,
                an,
                gff_in_builder.build(an),
            )
            for an in tsv[an_column]
        ]

        log.info("Tasks submitted, waiting for completion...")
        # just wait all tasks complete first
        all_frames = []
        for future in cf.as_completed(tasks):
            all_frames.append(future.result())

    df = pl.concat(all_frames)
    df = df.sort(["AN"], maintain_order=True)
    df.write_csv(output_file, separator="\t")

    log.info("All processed.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compile Summaries of intergenic runs."
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        required=True,
        help="Path to TSV file with required columns.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=False,
        default=Path.cwd() / "summary_igs_intergenic.tsv",
        help="Path to output TSV file. Default is igs_summary.tsv",
    )
    parser.add_argument(
        "--workers",
        type=int,
        required=False,
        default=MAX_CPU,
        help="Number of worker processes to use. Default is number of "
        f"CPU cores ({MAX_CPU}).",
    )
    parser.add_argument(
        "--gff-in-ext",
        type=str,
        required=False,
        default=".gff3",
        help="Extension for GFF input files. Default is .gff3",
    )
    parser.add_argument(
        "--gff-in-suffix",
        type=str,
        required=False,
        default="_intergenic",
        help="Suffix for GFF input files. Default is _intergenic",
    )
    parser.add_argument(
        "--an-column",
        type=str,
        required=False,
        default="AN",
        help="Column name for AN in the TSV file.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file.",
    )
    args = parser.parse_args()

    args.tsv = Path(args.tsv).resolve()
    args.output = Path(args.output).resolve()

    if args.output.exists() and not args.overwrite:
        raise FileExistsError(
            f"Output file {args.output} already exists. Use --overwrite to overwrite."
        )

    # create simple logger, no file, just console
    log = tigre.create_logger(
        print_to_console=True,
        console_level="INFO",
        save_to_file=False,
    )

    igs_multiple(
        log,
        args.tsv,
        args.output,
        args.workers,
        args.gff_in_ext,
        args.gff_in_suffix,
        args.an_column,
    )


if __name__ == "__main__":
    main()
