from __future__ import annotations

from pathlib import Path
import gzip
import shutil
import subprocess
from collections import Counter
import contextlib
import functools
import json
import threading
import time

import requests
import numpy as np
import pandas as pd
import psutil

import config


def save_file_from_url(source_url: str) -> Path:
    """
    Save a file from the given source to the given destination.

    If file is found to be a known zipped format, unzip it.

    Returns a Path to the resulting file.
    
    TODO: Make this stream files and generally be more efficient. Also note that this might not work for some URLs because HTTP headers are a nightmare.
    """
    response = requests.get(source_url)

    try:
        # TODO: This currently assumes that filename is always the second entry in the header; consider searching for "filename" instead
        header = response.headers["content-disposition"]
    except KeyError:
        # default to generating file name manually from the source url
        file_name = source_url.split("/")[-1]
    else:
        # This worked for a long time but needed to fix it due to NCBI responses
        # file_name = next(iter(re.findall("filename=(.+)", header)))
        file_name = header.split(";")[1].split("=")[1].strip('"')
    
    destination_file = config.raw_data_dir / file_name
    
    with destination_file.open("wb") as fp:
        fp.write(response.content)
    
    match destination_file.suffix:
        case ".gz":
            destination_file = gunzip(destination_file)
    
    return destination_file

def gunzip(gz_file: Path) -> Path:
    """
    Unzip the given file, deleting the original.

    Return a Path to the resulting file.
    """
    unzipped_file = gz_file.with_suffix("")
    with gzip.open(gz_file, "rb") as gz_fp, unzipped_file.open("wb") as unzipped_fp:
        shutil.copyfileobj(gz_fp, unzipped_fp)
    
    gz_file.unlink()

    return unzipped_file


def lift_over(
    old_file: Path,
    file_format: str,
    chain_file: Path,
    new_assembly_desc: str,
) -> Path:
    """
    Run the UCSC liftOver tool on the given inputs.

    Needs to be in bed, gff, or bedGraph format. Need to specify format manually because some downloaded files have non-standard names.

    Return a Path to the resulting lifted-over file.
    """    
    cmd_args = [config.user_config["executables"]["liftOver_exe"]]

    match file_format:
        case "bed":
            # BED files are often malformed; only guarantee that the first 6 columns are correct
            cmd_args.append("-bedPlus=6")
        case "gff":
            cmd_args.append("-gff")
        case "bedGraph":
            pass
        case _:
            raise ValueError("Format must be one of the strings 'bed', 'gff', or 'bedGraph'")
    
    new_file = config.processed_data_dir / old_file.with_suffix(f".{new_assembly_desc}.{file_format}").name
    unmapped_file = new_file.with_suffix(f".unmapped.{file_format}")
    cmd_args.extend([old_file, chain_file, new_file, unmapped_file])
    
    # Removed command printing to prevent printing of local paths in notebooks
    # print(f"liftOver command: {' '.join(str(x) for x in cmd_args)}")
    subprocess.run(cmd_args, check=True)

    print(f"lifted file: {new_file.relative_to(config.working_dir)}")

    with old_file.open("rb") as old_fp:
        old_count = sum(1 for _ in old_fp)
    with new_file.open("rb") as new_fp:
        new_count = sum(1 for _ in new_fp)
    print(f"original: {old_count}; lifted: {new_count}; percent: {new_count / old_count}")

    with unmapped_file.open("r") as unmapped_fp:
        unmapped_counts = Counter(line for idx, line in enumerate(unmapped_fp) if idx % 2 == 0)
    print(f"total unmapped: {unmapped_counts.total()}")
    for k, v in unmapped_counts.items():
        print(k, v)

    return new_file


def lift_over_bigWig(
    old_file: Path,
    chain_file: Path,
    chrom_size_file: Path,
    new_assembly_desc: str,
) -> Path:
    """
    Run the UCSC liftOver tool on a bedgraph file. Separate from the lift_over method because 
    this is a more involved process with multiple steps.

    Based on instructions from https://bioinfocore.com/blogs/liftover-bigwig-files/.
    """
    intermediate_bedgraph = old_file.with_suffix(".bedGraph")
    # !{ucsc_utils_dir}/bigWigToBedGraph {input_bw} {intermediate_bedgraph}
    subprocess.run(
        [config.user_config["executables"]["bigWigToBedGraph_exe"], old_file, intermediate_bedgraph],
        check=True
    )
    lifted_bedgraph = lift_over(intermediate_bedgraph, "bedGraph", chain_file, new_assembly_desc)

    lifted_sorted_bedgraph = lifted_bedgraph.with_suffix(".sorted.bedGraph")
    # !LC_COLLATE=C sort -k1,1 -k2,2n {lifted_bedgraph} > {lifted_sorted_bedgraph}
    with lifted_sorted_bedgraph.open("w") as fp:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", lifted_bedgraph],
            stdout=fp, encoding="utf-8"
        )
    
    output_bw = lifted_bedgraph.with_suffix(f".bigWig")
    # !{ucsc_utils_dir}/bedGraphToBigWig {lifted_sorted_bedgraph} {target_reference_id}.chrom.sizes {output_bw}
    subprocess.run(
        [config.user_config["executables"]["bedGraphToBigWig_exe"], lifted_sorted_bedgraph, chrom_size_file, output_bw],
        check=True
    )
    
    intermediate_bedgraph.unlink()
    lifted_bedgraph.unlink()
    lifted_sorted_bedgraph.unlink()


def load_bed_file(
    bed_file: Path,
    extra_columns: list = None
) -> pd.DataFrame:
    """
    Loads a bed file. Expected to be BED6+ format.
    """
    if extra_columns == None:
        extra_columns = []
    columns = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
    ] + extra_columns

    # Verify that the file matches the given set of columns
    with bed_file.open("r") as fp:
        ncols_first_line = len(next(fp).strip().split("\t"))
        if ncols_first_line != len(columns):
            raise ValueError("Mismatch between expected and observed number of columns; check that you are loading the correct type of BED file")

    return pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        index_col=False,
        names=columns
    )

def load_encode_narrowpeak_bed(bed_file: Path) -> pd.DataFrame:
    """
    Loads an ENCODE narrowPeak file, a bespoke BED6+4 format.
    """
    table = load_bed_file(
        bed_file,
        extra_columns=[
            "signalValue",
            "pValue",
            "qValue",
            "peak"
        ]
    )
    # TODO: I think it would be nice to do this replacement, but then writing bed files gets gnarly. Consider coming back to this.
    # table[["name", "strand"]] = table[["name", "strand"]].replace(".", None)
    # table[["pValue", "qValue", "peak"]] = table[["pValue", "qValue", "peak"]].replace(-1, np.nan)
    return table

def load_broadpeak_bed(bed_file: Path) -> pd.DataFrame:
    """
    Loads a broadPeak file, a bespoke BED6+3 format.
    """
    # TODO: This is the same as narrowPeak, just missing the last column. Maybe it can be consolidated.
    table = load_bed_file(
        bed_file,
        extra_columns=[
            "signalValue",
            "pValue",
            "qValue"
        ]
    )
    # TODO: I think it would be nice to do this replacement, but then writing bed files gets gnarly. Consider coming back to this.
    # table[["name", "strand"]] = table[["name", "strand"]].replace(".", None)
    # table[["pValue", "qValue", "peak"]] = table[["pValue", "qValue", "peak"]].replace(-1, np.nan)
    return table

def write_bed_file(table: pd.DataFrame, out_path: Path) -> None:
    """
    Write the bed file contained in the given DataFrame to the specified path.
    """
    table.to_csv(out_path, index=None, header=None, sep="\t")

def split_bed_by_quantiles(
    table: pd.DataFrame,
    column: str,
    q: int,
    dataset_id: str,
    avoid_duplicates: bool=False
) -> None:
    """
    Split a given bed into the desired number of quantiles, then write the resulting bed subsets to files.

    Makes a directory based on the dataset id in the processed data directory. Files are named after their quantile index.

    Args:
    * avoid_duplicates: if True, uses pd.DataFrame.rank to ensure there are no duplicate boundaries
    """
    output_dir = config.processed_data_dir / dataset_id
    output_dir.mkdir(exist_ok=True)

    labels = [str(i + 1) for i in range(q)]

    if avoid_duplicates:
        s = table[column].rank(method="first")
    else:
        s = table[column]
    quantile_assignments = pd.qcut(s, q=q, labels=labels)

    # In this case, observed True or False should not make a difference, but specifying True to squash FutureWarnings    
    for qidx, qgroup in table.groupby(quantile_assignments, observed=True):
        write_bed_file(qgroup, output_dir / f"q{qidx}.bed")

def split_bed_by_2_column_rank(
    table: pd.DataFrame,
    col1: str,
    col2: str,
    dataset_id: str
) -> None:
    """
    Split a given bed based on the 2-quantile ranking of both columns, then write the resulting bed subsets to files.
    
    In other words, for 2 columns col1 and col2, outputs the set of low/low, low/high, high/low, and high/high values, where high and low are above and below 50% for each column respectively.
    """
    output_dir = config.processed_data_dir / dataset_id
    output_dir.mkdir(exist_ok=True)

    cols = [col1, col2]
    labels = ["low", "high"]

    quantile_assignments = pd.concat((pd.qcut(table[col], q=2, labels=labels) for col in cols), axis=1)
    # In this case, observed True or False should not make a difference, but specifying True to squash FutureWarnings
    quantile_groups = quantile_assignments.groupby(cols, observed=True)

    print(quantile_groups.size())

    for (col1_rank, col2_rank), group in quantile_groups:
        write_bed_file(table.loc[group.index], output_dir / f"{col1}_{col2}.{col1_rank}_{col2_rank}.bed")

def apply_bed_offset(
    table: pd.DataFrame,
    offset: int,
    regions_5to3prime: bool = True
) -> pd.DataFrame:
    """
    Applies the specified offset in bp to the regions in the bed table.

    If regions_5to3prime is True, preserves strandedness.

    Examples:
        * offset=100, regions_5to3prime=True --> offset 100bp downstream
        * offset=-100, regions_5to3prime=True --> offset 100bp upstream
        * offset=100, regions_5to3prime=False --> offset 100bp in a positive direction
        * offset=-100, regions_5to3prime=False --> offset 100bp in a negative direction
    
    NOTE: This does not check for boundary conditions. For current needs this is unnecessary, but take note.
    """
    if regions_5to3prime:
        def offset_apply_fun(row):
            if row["strand"] == "+":
                return row[["chromStart", "chromEnd"]] + offset
            else:
                return row[["chromStart", "chromEnd"]] - offset
    else:
        def offset_apply_fun(row):
            return row[["chromStart", "chromEnd"]] + offset

    new_starts_ends = table.apply(offset_apply_fun, axis=1)
    return table.assign(
        chromStart=new_starts_ends["chromStart"],
        chromEnd=new_starts_ends["chromEnd"],
    )

def bed_peak_centers(table: pd.DataFrame) -> pd.Series:
    """
    Calculate the center of the peaks from the given bed file table.

    Current implementation is simple: start + ((end - start) / 2).
    This will result in floats, which may not be appropriate for many usecases, but should be fine for now.
    """
    centers = table["chromStart"] + ((table["chromEnd"] - table["chromStart"]) / 2)
    return pd.DataFrame({"chrom": table["chrom"], "peakCenter": centers})

def score_peaks_by_distance(
    source_peak_table: pd.DataFrame,
    query_peak_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Returns a minimal narrowPeak format bed table, replacing signalValue with this scoring
    """
    source_peak_centers = bed_peak_centers(source_peak_table)
    query_peak_centers = bed_peak_centers(query_peak_table)

    # TODO: This is done in a very dumb way. Not sure what the better way is...
    chrom_distance_vectors = []
    for chrom, chrom_source_centers in source_peak_centers.groupby("chrom"):
        chrom_distance_vectors.append(chrom_source_centers["peakCenter"].apply(
            lambda center: (query_peak_centers.query(f"chrom=='{chrom}'")["peakCenter"] - center).abs().min()
        ))
    source_peak_distances = pd.concat(chrom_distance_vectors)

    return source_peak_table.assign(signalValue=source_peak_distances)

"""Lightweight wall-time + peak-memory tracking for reproducibility reporting.

Measures the whole process tree (this process + all descendants), so it
correctly captures memory used by subprocesses -- including multi-threaded
compiled binaries (Rust/C) invoked via subprocess, and multi-process Python
workers.

Requires psutil:  pip install psutil

Usage
-----
# As a context manager:
    with track("pileup"):
        parse_bam.pileup(...)

# As a decorator (label defaults to the function name):
    @track(log_path="perf.jsonl")
    def run_pileup(...):
        parse_bam.pileup(...)

# Reading the numbers back out:
    with track("extract", quiet=True) as t:
        parse_bam.extract(...)
    print(t.elapsed, t.peak_gb)
"""

class track(contextlib.ContextDecorator):
    """Track wall-clock time and peak resident memory of the process tree.

    Parameters
    ----------
    label : str, optional
        Name for this measurement. When used as a decorator and left as None,
        defaults to the wrapped function's name.
    interval : float
        Polling interval in seconds. Smaller catches briefer spikes at a small
        CPU cost. 0.05 is fine for minute-scale jobs.
    log_path : str | Path, optional
        If given, append one JSON line per measurement to this file.
    quiet : bool
        If True, don't print a summary line on exit.
    """

    def __init__(self, label=None, *, interval=0.05, log_path=None, quiet=False):
        self.label = label
        self.interval = interval
        self.log_path = Path(log_path) if log_path else None
        self.quiet = quiet

    # Let the decorator form pick up the function name as the default label.
    def __call__(self, func):
        if self.label is None:
            self.label = func.__name__

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with self._fresh() as m:
                result = func(*args, **kwargs)
            # expose last measurement on the wrapper for convenience
            wrapper.last = m
            return result

        return wrapper

    # Each `with` / call gets its own state so repeated use is independent.
    def _fresh(self):
        return track(
            label=self.label,
            interval=self.interval,
            log_path=self.log_path,
            quiet=self.quiet,
        )

    def _tree_rss(self):
        total = self._proc.memory_info().rss
        for child in self._proc.children(recursive=True):
            try:
                total += child.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
        return total

    def _poll(self):
        while not self._stop.is_set():
            self.peak_bytes = max(self.peak_bytes, self._tree_rss())
            self._stop.wait(self.interval)

    def __enter__(self):
        self._proc = psutil.Process()
        self.peak_bytes = self._tree_rss()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._poll, daemon=True)
        self._thread.start()
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.elapsed = time.perf_counter() - self._t0
        self._stop.set()
        self._thread.join()
        # final sample, in case the peak occurred right before exit
        self.peak_bytes = max(self.peak_bytes, self._tree_rss())
        self.peak_gb = self.peak_bytes / 1e9        # decimal GB
        self.peak_gib = self.peak_bytes / 2 ** 30   # binary GiB

        if not self.quiet:
            print(f"[perf] {self.label}: {self.elapsed:.1f}s  "
                  f"peak {self.peak_gb:.2f} GB")

        if self.log_path:
            rec = {
                "label": self.label,
                "seconds": round(self.elapsed, 3),
                "peak_gb": round(self.peak_gb, 4),
                "peak_gib": round(self.peak_gib, 4),
                "timestamp": time.time(),
            }
            with open(self.log_path, "a") as f:
                f.write(json.dumps(rec) + "\n")

        return False  # don't suppress exceptions
