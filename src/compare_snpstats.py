#!/usr/bin/env python3
"""
Compare per-position SNP statistics across samples and flag changes.

Input is two or more .snpstats.tsv files produced by so2snpstats.py.
For indel comparison see compare_indelstats.py.

Flags
-----
  SNP_GAIN         SNP position present in some samples but absent in others.
  SNP_FLIP         Major alt allele changed between samples (e.g. A→T).
  SNP_FREQ_SHIFT   max(major_alt_freq) - min(major_alt_freq) >= --freq-shift.

Multiple flags are pipe-separated. Rows are sorted: flagged first, then by
freq_range descending, then seqid / pos.

Absent-sample convention
------------------------
When a sample has no SNP at a given position its numeric columns (cov, counts,
frequencies) are written as 0 in the output. major_alt is left empty.
This means a SNP at 80% in sample A that is absent in sample B produces
freq_range = 0.8 and is flagged as SNP_FREQ_SHIFT.

Usage
-----
    REVEAL snpcompare --snpstats s1.snpstats.tsv s2.snpstats.tsv
    REVEAL snpcompare --snpstats *.snpstats.tsv --freq-shift 0.8 --flagged-only
    REVEAL snpcompare --snpstats *.snpstats.tsv --outfile snp_comparison.tsv

Output defaults to comparison.snps.tsv.

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

PER_SAMPLE = [
    "cov",
    "A", "T", "C", "G",
    "A_freq", "T_freq", "C_freq", "G_freq",
    "major_alt", "major_alt_freq",
]


# ── loading ───────────────────────────────────────────────────────────────────

def load_snpstats(files: list[Path]) -> pd.DataFrame:
    frames = []
    for fp in files:
        if not fp.is_file():
            log.error("File not found: %s", fp)
            sys.exit(1)
        log.info("Loading %s", fp.name)
        frames.append(pd.read_csv(fp, sep="\t"))

    combined = pd.concat(frames, ignore_index=True)

    dupes = combined.duplicated(subset=["seqid", "pos", "sampleid"])
    if dupes.any():
        log.warning(
            "%d duplicate (seqid, pos, sampleid) pairs — keeping first occurrence. "
            "Re-run so2snpstats to fix the source files:\n%s",
            int(dupes.sum()),
            combined[dupes][["seqid", "pos", "sampleid"]].to_string(index=False),
        )
        combined = combined[~dupes].reset_index(drop=True)
    return combined


# ── comparison ────────────────────────────────────────────────────────────────

def compare_snps(long_df: pd.DataFrame, samples: list[str], freq_shift: float) -> pd.DataFrame:
    head       = ["seqid", "pos", "refc", "flag", "flagged", "freq_range"]
    per_sample = [f"{m}__{s}" for m in PER_SAMPLE for s in samples]

    if long_df.empty:
        return pd.DataFrame(columns=head + per_sample)

    # One row per (seqid, pos); columns: {metric}__{sampleid}
    wide = long_df.pivot(index=["seqid", "pos"], columns="sampleid", values=PER_SAMPLE)
    wide.columns = [f"{m}__{s}" for m, s in wide.columns]
    wide = wide.reset_index()

    refc = long_df.groupby(["seqid", "pos"])["refc"].first().reset_index()
    wide = wide.merge(refc, on=["seqid", "pos"], how="left")

    # ── flags ──
    major_alt_cols = [f"major_alt__{s}" for s in samples if f"major_alt__{s}" in wide.columns]
    freq_cols      = [f"major_alt_freq__{s}" for s in samples if f"major_alt_freq__{s}" in wide.columns]

    # n_present computed before fill so NaN correctly marks "absent in that sample"
    n_present = wide[major_alt_cols].notna().sum(axis=1)
    # treat absent sample as freq=0 so SNP_FREQ_SHIFT detects gain/loss correctly
    freq_filled = wide[freq_cols].fillna(0)
    freq_range  = freq_filled.max(axis=1) - freq_filled.min(axis=1)

    flags = pd.DataFrame(index=wide.index)
    flags["SNP_GAIN"]       = (n_present > 0) & (n_present < len(samples))
    flags["SNP_FLIP"]       = (n_present >= 2) & (wide[major_alt_cols].nunique(axis=1, dropna=True) > 1)
    flags["SNP_FREQ_SHIFT"] = freq_range >= freq_shift

    wide["freq_range"] = freq_range.round(4)
    wide["flag"]       = flags.apply(lambda r: "|".join(k for k, v in r.items() if v), axis=1)
    wide["flagged"]    = wide["flag"].str.len() > 0

    # fill absent-sample numeric columns with 0 (absent = 0 frequency/count)
    numeric_cols = [f"{m}__{s}" for m in PER_SAMPLE if m != "major_alt"
                    for s in samples if f"{m}__{s}" in wide.columns]
    wide[numeric_cols] = wide[numeric_cols].fillna(0)

    # ── column order ──
    major_cols = [f"{m}__{s}" for m in ("major_alt", "major_alt_freq") for s in samples]
    other_cols = [c for c in per_sample if c not in set(major_cols)]
    present_per_sample = [c for c in major_cols + other_cols if c in wide.columns]
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
    print(f"\nSNPs: {n_flagged} / {len(df)} positions flagged  →  {outfile}")
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
        description="Compare per-position SNP stats across samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--snpstats", nargs="+", required=True, type=Path, metavar="FILE",
        help="Two or more .snpstats.tsv files produced by so2snpstats.py.",
    )
    parser.add_argument(
        "--outfile", "-o", type=Path, default=None, metavar="FILE",
        help="Output TSV (default: comparison.snps.tsv).",
    )
    parser.add_argument(
        "--flagged-only", action="store_true",
        help="Only output rows that have at least one flag.",
    )
    parser.add_argument(
        "--freq-shift", type=float, default=0.8, metavar="FLOAT",
        dest="freq_shift",
        help=(
            "Frequency shift threshold for SNP_FREQ_SHIFT. "
            "Flag when max(major_alt_freq) - min(major_alt_freq) >= this value. "
            "Default: 0.8 (80%% swing)."
        ),
    )

    args = parser.parse_args()

    if len(args.snpstats) < 2:
        parser.error("--snpstats requires at least two files to compare.")

    outfile = args.outfile or Path("comparison.snps.tsv")

    long_df = load_snpstats(args.snpstats)
    samples = sorted(long_df["sampleid"].unique().tolist())
    log.info("Samples (%d): %s", len(samples), ", ".join(samples))

    wide = compare_snps(long_df, samples, args.freq_shift)

    if args.flagged_only and "flagged" in wide.columns:
        n_before = len(wide)
        wide = wide[wide["flagged"]]
        log.info("--flagged-only: %d / %d retained", len(wide), n_before)

    wide.to_csv(outfile, sep="\t", index=False)
    log.info("SNP comparison → %s  (%d positions)", outfile, len(wide))
    _print_summary(wide, outfile)


if __name__ == "__main__":
    main()
