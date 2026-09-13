#!/usr/bin/env python3
"""Extract framewise-displacement summaries from fMRIPrep confounds files.

Walks a directory for ``*_desc-confounds_timeseries.tsv`` files, reads the
``framewise_displacement`` column out of each one, and writes a long-format
table with one row per subject x run:

    ID,run,day,n_vols,mean_FD,max_FD

The layout mirrors ``tSNR_qc.csv`` so that the two QC tables can be joined on
``ID``/``run``/``day`` and analysed with the same code. Unlike tSNR there is no
echo dimension: fMRIPrep estimates motion once per BOLD series, so FD is a
property of the run, not of the individual echoes.

Usage
-----
    python3 extract_fd.py /path/to/fmriprep -o fd_results.csv

    # also emit the wide layout used by tSNR_results.txt (ID, r1, r2, ... r6)
    python3 extract_fd.py /path/to/fmriprep -o fd_results.csv --wide fd_wide.txt

    # add columns counting volumes above an FD threshold, in mm
    python3 extract_fd.py . --fd-threshold 0.5
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

FD_COLUMN = "framewise_displacement"
DEFAULT_PATTERN = "*rest*_desc-confounds_timeseries.tsv"

# Run -> acquisition day, from eeg_fmri_ext_SNR_1_v2.m:63 (rest_day = [1 2 2 2 2 3]).
# Only used to cross-check the ses-<label> entity in the filename; the filename
# wins, so a dataset with a different layout still parses correctly.
EXPECTED_RUN_DAY = {1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 3}


# --------------------------------------------------------------------------- #
# filename parsing
# --------------------------------------------------------------------------- #

def parse_entities(path: Path) -> dict:
    """Split a BIDS-style filename into its key-value entities."""
    entities = {}
    for part in path.name.split("_"):
        if "-" in part:
            key, _, value = part.partition("-")
            entities[key] = value
    return entities


def parse_subject(entities: dict, strip_letters: bool) -> str:
    """Return the subject label, e.g. ``sub-10``.

    fMRIPrep filenames in this dataset use ``sub-S10`` while the tSNR tables use
    ``sub-10``, so by default any letters are dropped from the label to keep the
    two joinable. ``--keep-subject-label`` disables that.
    """
    label = entities.get("sub")
    if label is None:
        return ""
    if not strip_letters:
        return f"sub-{label}"
    digits = re.sub(r"\D", "", label)
    return f"sub-{digits}" if digits else f"sub-{label}"


def parse_run(entities: dict):
    """Run number, from the ``run`` entity or from digits in the task label.

    This dataset encodes the run in the task name (``task-rest1`` ... ``rest6``)
    rather than in a ``run`` entity, so both are supported.
    """
    if "run" in entities:
        digits = re.sub(r"\D", "", entities["run"])
        if digits:
            return int(digits)
    task = entities.get("task", "")
    match = re.search(r"(\d+)$", task)
    if match:
        return int(match.group(1))
    return None


def parse_day(entities: dict):
    """Acquisition day, from digits in the ``ses`` entity (``ses-Day2`` -> 2)."""
    match = re.search(r"(\d+)", entities.get("ses", ""))
    return int(match.group(1)) if match else None


# --------------------------------------------------------------------------- #
# FD extraction
# --------------------------------------------------------------------------- #

def read_fd(path: Path):
    """Return the FD column as a list of floats, dropping non-numeric entries.

    fMRIPrep writes ``n/a`` for the first volume, which has no preceding volume
    to be displaced from. That frame carries no information and is excluded from
    both the mean and the max, so ``n_vols`` is the number of usable FD values,
    one fewer than the number of volumes in the run.
    """
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("file is empty")
        if FD_COLUMN not in header:
            raise ValueError(f"no {FD_COLUMN!r} column")
        index = header.index(FD_COLUMN)

        values = []
        for row in reader:
            if index >= len(row):
                continue
            try:
                value = float(row[index])
            except ValueError:
                continue  # 'n/a', '', 'NaN'
            if not math.isnan(value):
                values.append(value)
    return values


def summarise(path: Path, strip_letters: bool, threshold):
    """One output record for one confounds file, or None if FD is unusable."""
    entities = parse_entities(path)
    fd = read_fd(path)
    if not fd:
        raise ValueError(f"no usable {FD_COLUMN} values")

    run = parse_run(entities)
    day = parse_day(entities)

    record = {
        "ID": parse_subject(entities, strip_letters),
        "run": run,
        "day": day,
        "n_vols": len(fd),
        "mean_FD": round(sum(fd) / len(fd), 6),
        "max_FD": round(max(fd), 6),
    }
    if threshold is not None:
        above = sum(1 for value in fd if value > threshold)
        record[f"n_FD_gt_{threshold:g}"] = above
        record[f"pct_FD_gt_{threshold:g}"] = round(100.0 * above / len(fd), 3)

    # Echo is kept only if fMRIPrep actually wrote per-echo confounds; see the
    # duplicate check in main().
    record["_echo"] = entities.get("echo")
    record["_file"] = path.name
    record["_ses"] = entities.get("ses")
    return record


# --------------------------------------------------------------------------- #
# output
# --------------------------------------------------------------------------- #

def sort_key(record: dict):
    digits = re.sub(r"\D", "", record["ID"])
    return (int(digits) if digits else 0, record["ID"],
            record["run"] if record["run"] is not None else 0)


def write_long(records, out_path: Path, columns):
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def write_wide(records, out_path: Path, measure: str):
    """Wide layout matching tSNR_results.txt: one row per subject, one column
    per run (``r1`` ... ``r6``), missing cells as ``NaN``."""
    runs = sorted({r["run"] for r in records if r["run"] is not None})
    subjects = sorted({r["ID"] for r in records},
                      key=lambda s: (int(re.sub(r"\D", "", s) or 0), s))
    table = {(r["ID"], r["run"]): r[measure] for r in records}

    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ID"] + [f"r{run}" for run in runs])
        for subject in subjects:
            writer.writerow(
                [subject] + [table.get((subject, run), "NaN") for run in runs]
            )


# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Extract mean and maximal framewise displacement from "
                    "fMRIPrep confounds files into a long-format table.")
    parser.add_argument("directory", nargs="?", default=".", type=Path,
                        help="directory to search (default: current)")
    parser.add_argument("-o", "--out", default="fd_results.csv", type=Path,
                        help="long-format output CSV (default: fd_results.csv)")
    parser.add_argument("--wide", type=Path, default=None,
                        help="also write this wide-format file, one column per "
                             "run, as in tSNR_results.txt")
    parser.add_argument("--wide-measure", default="mean_FD",
                        choices=["mean_FD", "max_FD"],
                        help="which measure the wide file holds "
                             "(default: mean_FD)")
    parser.add_argument("--pattern", default=DEFAULT_PATTERN,
                        help=f"glob for confounds files (default: {DEFAULT_PATTERN})")
    parser.add_argument("--no-recursive", action="store_true",
                        help="do not descend into subdirectories")
    parser.add_argument("--fd-threshold", type=float, default=None,
                        metavar="MM",
                        help="also count volumes with FD above this value in mm "
                             "(0.5 is the usual scrubbing threshold)")
    parser.add_argument("--keep-subject-label", action="store_true",
                        help="keep the subject label verbatim (sub-S10) instead "
                             "of stripping letters to match the tSNR tables")
    args = parser.parse_args(argv)

    if not args.directory.is_dir():
        parser.error(f"not a directory: {args.directory}")

    files = sorted(args.directory.glob(args.pattern) if args.no_recursive
                   else args.directory.rglob(args.pattern))
    if not files:
        print(f"No files matching {args.pattern!r} under {args.directory}",
              file=sys.stderr)
        return 1

    records, failures = [], 0
    for path in files:
        try:
            records.append(
                summarise(path, not args.keep_subject_label, args.fd_threshold))
        except (ValueError, OSError) as error:
            print(f"SKIP {path.name}: {error}", file=sys.stderr)
            failures += 1

    if not records:
        print("No usable FD data found.", file=sys.stderr)
        return 1

    # Warn about anything that would make the table ambiguous.
    for record in records:
        if record["run"] is None:
            print(f"WARNING no run number in {record['_file']}", file=sys.stderr)
        if record["day"] is None:
            print(f"WARNING no session/day in {record['_file']}", file=sys.stderr)
        elif EXPECTED_RUN_DAY.get(record["run"]) not in (None, record["day"]):
            print(f"WARNING {record['_file']}: run {record['run']} is labelled "
                  f"day {record['day']}, expected day "
                  f"{EXPECTED_RUN_DAY[record['run']]}", file=sys.stderr)

    # Per-echo confounds files would give several identical FD series per run.
    # Collapse them, but only after checking that they really are identical --
    # if they differ, motion was estimated per echo and that is worth knowing.
    seen = {}
    collapsed = []
    for record in records:
        key = (record["ID"], record["run"])
        if key in seen:
            first = seen[key]
            if abs(first["mean_FD"] - record["mean_FD"]) > 1e-6:
                print(f"WARNING {record['_file']}: duplicate {key} with a "
                      f"different mean FD ({first['mean_FD']} vs "
                      f"{record['mean_FD']}); keeping the first",
                      file=sys.stderr)
            else:
                print(f"NOTE {record['_file']}: duplicate {key} "
                      f"(per-echo confounds); keeping the first", file=sys.stderr)
            continue
        seen[key] = record
        collapsed.append(record)

    collapsed.sort(key=sort_key)

    columns = [key for key in collapsed[0] if not key.startswith("_")]
    write_long(collapsed, args.out, columns)
    print(f"wrote {args.out}  ({len(collapsed)} rows, "
          f"{len({r['ID'] for r in collapsed})} subjects, "
          f"{len({r['run'] for r in collapsed})} runs)")

    if args.wide is not None:
        write_wide(collapsed, args.wide, args.wide_measure)
        print(f"wrote {args.wide}  ({args.wide_measure}, wide layout)")

    if failures:
        print(f"{failures} file(s) skipped", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
