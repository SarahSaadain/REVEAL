#!/usr/bin/env python3
"""
Compare per-position indel statistics across samples and flag changes.

Input is two or more .indelstats.tsv files produced by so2indelstats.py.
For SNP comparison see compare_snpstats.py.

Flags
-----
  INDEL_GAIN       Indel (same position + type + length) present in some
                   samples but absent in others.
  INDEL_FREQ_SHIFT max(indel_freq) - min(indel_freq) >= --freq-shift.

Multiple flags are pipe-separated. Rows are sorted: flagged first, then by
freq_range descending, then seqid / pos.

Absent-sample convention
------------------------
When a sample has no indel at a given (position, type, length) combination its
numeric columns (cov, count, freq) are written as 0 in the output.

Join key
--------
(seqid, pos, type, indel_length) — different-length events at the same position
are biologically distinct and compared separately.

Usage
-----
    REVEAL indelcompare --indelstats s1.indelstats.tsv s2.indelstats.tsv
    REVEAL indelcompare --indelstats *.indelstats.tsv --freq-shift 0.8 --flagged-only
    REVEAL indelcompare --indelstats *.indelstats.tsv --outfile indel_comparison.tsv

Output defaults to comparison.indels.tsv.

Authors
-------
    Robert Kofler
    Sarah Saadain
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

PER_SAMPLE = ["cov", "indel_count", "indel_freq"]


# ── loading ───────────────────────────────────────────────────────────────────

def load_indelstats(files: list[Path]) -> pd.DataFrame:
    frames = []
    for fp in files:
        if not fp.is_file():
            log.error("File not found: %s", fp)
            sys.exit(1)
        log.info("Loading %s", fp.name)
        frames.append(pd.read_csv(fp, sep="\t"))

    combined = pd.concat(frames, ignore_index=True)
    combined["indel_length"] = combined["indel_length"].astype("Int64")

    dupes = combined.duplicated(subset=["seqid", "pos", "type", "indel_length", "sampleid"])
    if dupes.any():
        log.warning(
            "%d duplicate (seqid, pos, type, indel_length, sampleid) pairs — keeping first occurrence. "
            "Re-run so2indelstats to fix the source files:\n%s",
            int(dupes.sum()),
            combined[dupes][["seqid", "pos", "type", "indel_length", "sampleid"]].to_string(index=False),
        )
        combined = combined[~dupes].reset_index(drop=True)
    return combined


# ── comparison ────────────────────────────────────────────────────────────────

def compare_indels(long_df: pd.DataFrame, samples: list[str], freq_shift: float) -> pd.DataFrame:
    head       = ["seqid", "pos", "type", "indel_length", "flag", "flagged", "freq_range"]
    per_sample = [f"{m}__{s}" for m in PER_SAMPLE for s in samples]

    if long_df.empty:
        return pd.DataFrame(columns=head + per_sample)

    # One row per (seqid, pos, type, indel_length); columns: {metric}__{sampleid}
    wide = long_df.pivot(
        index=["seqid", "pos", "type", "indel_length"],
        columns="sampleid",
        values=PER_SAMPLE,
    )
    wide.columns = [f"{m}__{s}" for m, s in wide.columns]
    wide = wide.reset_index()

    # ── flags ──
    freq_cols = [f"indel_freq__{s}" for s in samples if f"indel_freq__{s}" in wide.columns]
    n_present = wide[freq_cols].notna().sum(axis=1)
    # treat absent sample as freq=0
    freq_filled = wide[freq_cols].fillna(0)
    freq_range  = freq_filled.max(axis=1) - freq_filled.min(axis=1)

    flags = pd.DataFrame(index=wide.index)
    flags["INDEL_GAIN"]       = (n_present > 0) & (n_present < len(samples))
    flags["INDEL_FREQ_SHIFT"] = freq_range >= freq_shift

    wide["freq_range"] = freq_range.round(4)
    wide["flag"]       = flags.apply(lambda r: "|".join(k for k, v in r.items() if v), axis=1)
    wide["flagged"]    = wide["flag"].str.len() > 0

    # fill absent-sample numeric columns with 0
    fill_cols = [f"{m}__{s}" for m in PER_SAMPLE for s in samples if f"{m}__{s}" in wide.columns]
    wide[fill_cols] = wide[fill_cols].fillna(0)

    # ── column order ──
    present_per_sample = [c for c in per_sample if c in wide.columns]
    rest  = [c for c in wide.columns if c not in head + present_per_sample]
    wide  = wide[head + present_per_sample + rest]

    return wide.sort_values(
        ["flagged", "freq_range", "seqid", "pos"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


# ── summary printer ───────────────────────────────────────────────────────────

def _print_summary(df: pd.DataFrame, outfile: Path) -> None:
    if df.empty:
        return
    n_flagged = int(df["flagged"].sum())
    print(f"\nIndels: {n_flagged} / {len(df)} events flagged  →  {outfile}")
    all_flags = df["flag"].str.split("|").explode()
    all_flags = all_flags[all_flags.str.len() > 0]
    if not all_flags.empty:
        print(all_flags.value_counts().to_string())


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Compare per-position indel stats across samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--indelstats", nargs="+", required=True, type=Path, metavar="FILE",
        help="Two or more .indelstats.tsv files produced by so2indelstats.py.",
    )
    parser.add_argument(
        "--outfile", "-o", type=Path, default=None, metavar="FILE",
        help="Output TSV (default: comparison.indels.tsv).",
    )
    parser.add_argument(
        "--flagged-only", action="store_true",
        help="Only output rows that have at least one flag.",
    )
    parser.add_argument(
        "--freq-shift", type=float, default=0.8, metavar="FLOAT",
        dest="freq_shift",
        help=(
            "Frequency shift threshold for INDEL_FREQ_SHIFT. "
            "Flag when max(indel_freq) - min(indel_freq) >= this value. "
            "Default: 0.8 (80%% swing)."
        ),
    )

    args = parser.parse_args()

    if len(args.indelstats) < 2:
        parser.error("--indelstats requires at least two files to compare.")

    outfile = args.outfile or Path("comparison.indels.tsv")

    long_df = load_indelstats(args.indelstats)
    samples = sorted(long_df["sampleid"].unique().tolist())
    log.info("Samples (%d): %s", len(samples), ", ".join(samples))

    wide = compare_indels(long_df, samples, args.freq_shift)

    if args.flagged_only and "flagged" in wide.columns:
        n_before = len(wide)
        wide = wide[wide["flagged"]]
        log.info("--flagged-only: %d / %d retained", len(wide), n_before)

    wide.to_csv(outfile, sep="\t", index=False)
    log.info("Indel comparison → %s  (%d events)", outfile, len(wide))
    _print_summary(wide, outfile)


if __name__ == "__main__":
    main()
