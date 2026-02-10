#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compile Summaries of intergenic runs."""

import argparse
import concurrent.futures as cf
import os
import re
from pathlib import Path
from typing import Any, Hashable

import pandas as pd
import tigre

_RE_strand = re.compile(r"\|([+-])\|")
_RE_ID = re.compile(r"ID=([^;]+)")

MAX_CPU = os.cpu_count() or 1


def summary_igs(
    an: str,
    gff_path: Path,
) -> list[dict[Hashable, Any]]:
    """Compile Summary of IGS outputs."""
    df = pd.read_csv(
        gff_path,
        sep="\t",
        names=tigre.GFF3_COLUMNS,
        usecols=["type", "start", "end", "attributes"],
        comment="#",
    )
    df["ID"] = df["attributes"].str.extract(_RE_ID, expand=False)  # type: ignore[call-overload]
    df["Length"] = df["end"] - df["start"] + 1
    df["AN"] = an
    df["Merged"] = (df["type"] != "intergenic_region").astype(int)

    df["UP"] = df["attributes"].str.extract(tigre.clean._RE_name_up, expand=False)  # type: ignore[call-overload]
    df["UP"] = df["UP"].fillna("RS")
    df["DOWN"] = df["attributes"].str.extract(tigre.clean._RE_name_dw, expand=False)  # type: ignore[call-overload]
    df["DOWN"] = df["DOWN"].fillna("RE")

    df["source_up"] = df["attributes"].str.extract(tigre.clean._RE_source_up, expand=False)  # type: ignore[call-overload]
    df["source_dw"] = df["attributes"].str.extract(tigre.clean._RE_source_dw, expand=False)  # type: ignore[call-overload]

    df["Pair"] = df["UP"] + "-" + df["DOWN"]

    df["up-strand"] = df["source_up"].str.extract(_RE_strand, expand=False).fillna("+")  # type: ignore[call-overload]
    df["dw-strand"] = df["source_dw"].str.extract(_RE_strand, expand=False).fillna("+")  # type: ignore[call-overload]

    df["Polarity"] = df["up-strand"] + df["dw-strand"]

    # order columns by: AN UP DOWN Pair Polarity Length M source_up source_dw
    df = df[
        [
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
    ]

    return df.to_dict(orient="records")


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
    tsv = pd.read_csv(tsv_path, sep="\t")

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
        all_records = []
        for future in cf.as_completed(tasks):
            all_records.extend(future.result())

    df = pd.DataFrame(all_records)
    df = df.sort_values(by=["AN", "ID"]).reset_index(drop=True)
    df.to_csv(output_file, sep="\t", index=False)

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
