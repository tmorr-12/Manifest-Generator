#!/usr/bin/env python3

# TODO: output file containing all duplicates if duplicates are found

import argparse
import logging
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
from pandas._libs.missing import NAType

FASTQ_EXT = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

R1_PATTERN = re.compile(r"(?:_R?1)(?:_[\w]+)?$")
R2_PATTERN = re.compile(r"(?:_R?2)(?:_[\w]+)?$")

TIMESTAMP = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


class Fastq:
    def __init__(
        self,
        path: str | Path,
        ID: str | None = None,
        R1: Path | None = None,
        R2: Path | None = None,
        long_fastq: Path | None = None,
        genome_size: int | NAType = pd.NA,
    ) -> None:
        self.path = Path(path)
        self.ID = ID
        self.R1 = R1
        self.R2 = R2
        self.long_fastq = long_fastq
        self.genome_size = genome_size

    def _strip_extensions(self) -> str:
        """Remove extensions from file name"""
        name = self.path.name
        for ext in sorted(FASTQ_EXT, key=len, reverse=True):  # longest first, e.g. .fastq.gz before .gz
            if name.endswith(ext):
                return name[: -len(ext)]
        return name

    def _get_id(self) -> None:
        """Get read ID"""
        name = self._strip_extensions()
        name = R1_PATTERN.sub("", name)
        name = R2_PATTERN.sub("", name)
        self.ID = name

    def _detect_read_type(self) -> str:
        """Assign short reads to either R1 or R2, else assign long read"""
        name = self._strip_extensions()
        if R1_PATTERN.search(name):
            return "R1"
        if R2_PATTERN.search(name):
            return "R2"
        return "long_fastq"

    def estimate_genome_size(self, threads: int) -> int | None:
        """For use in long/hybrid manifests"""
        try:
            result = subprocess.run(
                ["lrge", "-t", f"{threads}", f"{self.path}"], check=True, text=True, capture_output=True
            )
            lines = result.stdout.strip().splitlines()
            return int(lines[-1])
        except subprocess.CalledProcessError as e:
            logging.warning(f"lrge process failed with exit status {e.returncode}: {e.cmd}")

    def assign_properties(self, mode) -> None:
        """Fill out class properties"""
        if mode != "hybrid":
            self._get_id()

        read_type = self._detect_read_type()
        if read_type == "R1":
            self.R1 = self.path
        elif read_type == "R2":
            self.R2 = self.path
        elif read_type == "long_fastq":
            self.long_fastq = self.path

    def to_dict(self) -> dict:
        """Create a dictionary of class properties for addition to DataFrame"""
        return {
            "ID": self.ID,
            "R1": self.R1,
            "R2": self.R2,
            "long_fastq": self.long_fastq,
            "genome_size": self.genome_size,
        }


class ManifestParser:
    def __init__(self, df: pd.DataFrame, mode: str, method: str) -> None:
        self.df = df
        self.mode = mode
        self.method = method

    def _detect_duplicates(self) -> tuple[bool, int]:
        """Helper function to determine whether manifest contains duplicates"""
        id_list = self.df["ID"].to_list()
        id_set = set(id_list)

        n_list = len(id_list)
        n_set = len(id_set)

        multipliers = {"long": 1, "short": 2, "hybrid": 3}
        expected = multipliers[self.mode] * n_set
        n_dup = max(0, n_list - expected)

        if n_list != expected and n_dup > 0:
            return (True, n_dup)
        else:
            return (False, 0)

    def _no_duplicates(self) -> pd.DataFrame:
        """Handle case where no duplicates are found"""
        manifest = self.df.groupby("ID", as_index=False).first()  # Collapse DataFrame
        if self.mode != "short":
            manifest = manifest.dropna(subset=["long_fastq"])
        else:
            manifest = manifest.dropna()
        return manifest

    def _pair_reads(self) -> pd.DataFrame:
        """Pair short and hybrid reads"""
        cols = {"R1": "R1", "R2": "R2"}
        if self.mode == "hybrid":
            cols["long_fastq"] = "long_fastq"

        dfs = [self.df[["ID", col]].sort_values(col).dropna() for col in cols]

        merged = dfs[0]
        for df in dfs[1:]:
            merged = merged.merge(df, on="ID", how="outer")
        return merged

    def _duplicate_handling_workflow(self, n_dup: int) -> pd.DataFrame | None:
        """Subworkflow for duplicate handling"""
        if self.method == "error":
            logging.error(f"Detected {n_dup} unacceptable duplicates, aborting...")
            sys.exit(1)

        df = self.df
        if self.mode != "long":
            paired_df = self._pair_reads()
            df = paired_df

        if self.method == "warn":
            logging.warning(f"Manifest contains {n_dup} duplicates")
            return df
        if self.method == "drop":
            logging.warning(f"Dropping {n_dup} duplicates from manifest")
            return df.groupby("ID", as_index=False).first()

    def handle_duplicates(self) -> pd.DataFrame | None:
        """Handle duplicates based on user specified method"""
        incidence = self._detect_duplicates()
        presence = incidence[0]
        n_dup = incidence[1]

        if presence:
            return self._duplicate_handling_workflow(n_dup)
        else:
            return self._no_duplicates()

    def publish_duplicates(self) -> None:
        """Publish csv containing duplicates"""
        pass


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        force=True,
    )
    logging.info("Logging initialized.")


class CustomFormatter(argparse.RawDescriptionHelpFormatter):
    def __init__(self, prog) -> None:
        super().__init__(prog, max_help_position=100, width=150)


def parse_args():
    banner = (
        "\033[1;96m"
        + r"""
 __  __    __    _  _  ____  ____  ____  ___  ____     ___  ____  _  _  ____  ____    __   ____  _____  ____
(  \/  )  /__\  ( \( )(_  _)( ___)( ___)/ __)(_  _)   / __)( ___)( \( )( ___)(  _ \  /__\ (_  _)(  _  )(  _ \
 )    (  /(__)\  )  (  _)(_  )__)  )__) \__ \  )(    ( (_-. )__)  )  (  )__)  )   / /(__)\  )(   )(_)(  )   /
(_/\/\_)(__)(__)(_)\_)(____)(__)  (____)(___/ (__)    \___/(____)(_)\_)(____)(_)\_)(__)(__)(__) (_____)(_)\_)
"""
        + "\033[0m"
        + "\n\n\033[1;95mWellcome to manifest generator, the swiss army knife on manifest generation!\n\033[0m"
    )

    parser = argparse.ArgumentParser(description=banner, formatter_class=CustomFormatter)

    usage = parser.add_argument_group("usage")

    usage.add_argument(
        "--from_dir",
        nargs="+",
        help="Create a manifest from all fastq's in a given directory",
    )
    usage.add_argument(
        "--from_dir_recursive",
        nargs="+",
        help="Create a manifest from all fastq's in a given parent directory (requires --max_depth flag)",
    )
    usage.add_argument(
        "--from_paths",
        type=Path,
        help="Create a manifest from a list of read paths",
    )
    usage.add_argument(
        "--from_paths_id",
        type=Path,
        help="Create a manifest from a csv with ID and read paths (required for hybrid manifests)",
    )

    flags = parser.add_argument_group("flags")

    flags.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=["short", "long", "hybrid"],
        default="short",
        help="Build a short, long or hybrid manifest (default = short)",
    )
    flags.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Path to output directory (default = '.')",
    )
    flags.add_argument(
        "-n",
        "--name",
        type=str,
        default=f"{TIMESTAMP}_manifest.csv",
        help="Now \033[4mYOU\033[0m can name your very own manifest! (default = <TIMESTAMP>_manifest.csv)",
    )
    flags.add_argument(
        "-d",
        "--max_depth",
        type=int,
        default=None,
        help="Set the maximum tree depth for --from_dir_recursive",
    )
    flags.add_argument(
        "--genome_size",
        action="store_true",
        help="Enable genome size estimation for long read",
    )
    flags.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads used for genome size estimation (default = 1)",
    )
    flags.add_argument(
        "--duplicate_handling",
        type=str,
        choices=["drop", "warn", "error"],
        default="error",
        help="Specify method with which to handle duplicate sample IDs in the manifest (default = error)",
    )
    flags.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="enable verbose logging",
    )
    return parser.parse_args()


def validate_args(args) -> argparse.Namespace:
    """
    Helper function to validate command line arguments

    Args:
        args: User supplied CLI arguments

    Returns:
        args

    Raises:
        Will run sys.exit(2) upon failed argument validation
    """
    if not any([args.from_dir, args.from_dir_recursive, args.from_paths, args.from_paths_id]):
        logging.error("one of --from_dirs, --from_dir_recursive, --from_paths or --from_paths_id is required")
        sys.exit(2)

    if args.from_dir_recursive and args.max_depth is None:
        logging.error("--max_depth must be supplied with --from_dir_recursive")
        sys.exit(2)

    for arg in (args.from_dir, args.from_dir_recursive):
        if arg is None:
            continue
        if isinstance(arg, str):
            arg = [arg]
        for dir in arg:
            if not Path(dir).is_dir():
                logging.error(f"{dir} is not a valid directory")
                sys.exit(2)

    if not args.name.endswith(".csv"):
        args.name += ".csv"

    if args.mode == "hybrid" and not args.from_paths_id:
        logging.error("Hybrid manifests can only be built from --from_paths_id")
        sys.exit(1)

    if args.from_paths_id and args.mode != "hybrid":
        logging.error("'--from_paths_id' is only for hybrid manifests, please try again with '--from_paths'")
        sys.exit(1)

    if args.genome_size:
        try:
            subprocess.run(["lrge", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception as e:
            logging.error(f"Could not run lrge, please load lrge and try again -> {e}")
            sys.exit(1)

    return args


def get_headers(mode: str) -> list[str]:
    """
    Get headers for dataframe

    Args:
        mode: Short/Long/Hybrid (args.mode)

    Returns:
        List columns headers for dataframe
    """
    if mode == "short":
        return ["ID", "R1", "R2"]
    if mode == "long":
        return ["ID", "long_fastq", "genome_size"]
    if mode == "hybrid":
        return ["ID", "R1", "R2", "long_fastq", "genome_size"]


def walk(root: Path, max_depth: int) -> list[Path]:
    """
    Search director(y/ies) for fastq files

    Args:
        root: Path to root directory
        max_depth: Max depth for recursive file search

    Returns:
        fastqs: List of fastq's
    """
    base_depth = str(root).count(os.sep)
    fastqs = []

    for path, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if d != "work"]  # skip nextflow work directories
        depth = path.count(os.sep) - base_depth

        if depth >= max_depth:
            dirs[:] = []

        for f in files:
            if f.endswith(FASTQ_EXT):
                fastqs.append(Path(path) / f)
    return fastqs


def collect_reads(dir_list: list[Path], max_depth: int) -> list[Path]:
    """
    Wrapper for walk function

    Args:
        dir_list: List of directories to search
        max_depth: Max depth for recursive dir search

    Returns:
        reads: List of Fastq's
    """
    reads = []
    for dir in dir_list:
        fastqs = walk(dir, max_depth)
        reads += fastqs
    return reads


def main():
    """
    Step 0: Parse args and setup logging

    Step 1: Collect a list of fastq file paths from all input flags

    Step 2: Initialise DataFrame and loop over main_list:
        - Extract ID and assign fastq read type (R1, R2, long_fastq)
        - Insert information into DataFrame

    Step 3: Check dataframe for duplicates and missing values
    """
    args = parse_args()

    setup_logging(args.verbose)

    validate_args(args)
    args.outdir.mkdir(parents=True, exist_ok=True)

    """
    Step 1
    """
    logging.info("Step 1: Gather fastq's")

    read_list = []

    if args.from_dir:  # add fastqs to read_list
        read_list += collect_reads(args.from_dir, 0)

    if args.from_dir_recursive:  # add fastqs to read_list
        read_list += collect_reads(args.from_dir_recursive, args.max_depth)

    if args.from_paths:  # add fastqs to read_list
        with open(args.from_paths) as f:
            read_list += [line.strip() for line in f]

    if args.from_paths_id:
        id_paths = pd.read_csv(args.from_paths_id)
        for _, row in id_paths.iterrows():
            read_list.append((row["ID"], row["path"]))

    """
    Step 2
    """
    logging.info("Step 2: Assign read properties")

    rows = []
    for read in read_list:
        if args.mode != "hybrid":
            read = Fastq(path=read)
        else:
            read = Fastq(ID=read[0], path=read[1])

        read.assign_properties(args.mode)
        if args.genome_size:
            logging.info(f"Estimating genome size for {read.ID}")
            size = read.estimate_genome_size(args.threads)
            read.genome_size = size

        rows.append(read.to_dict())

    headers = get_headers(args.mode)
    df = pd.DataFrame(rows, columns=headers).sort_values("ID")

    """
    Step 3
    """
    logging.info("Step 3: Process DataFrame")

    df = ManifestParser(df, args.mode, args.duplicate_handling).handle_duplicates()

    if df.empty:
        logging.error("Final manifest contains no data, exiting...")
        sys.exit(1)

    df = df.dropna(subset=[df.columns[1]])

    logging.info(f"Creating {args.mode} reads manifest with {len(df)} samples")

    outfile = args.outdir / f"{args.name}"

    logging.info(f"Writing manifest to {outfile}")

    df.to_csv(outfile, index=False)

    logging.info("Done! Thank you for using manifest generator")


if __name__ == "__main__":
    main()
