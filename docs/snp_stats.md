# SNP Stats: `so2snpstats` and `snpcompare`

Per-position SNP statistics and cross-sample comparison.

For indel statistics see [indel_stats.md](indel_stats.md).

---

## Overview

| Script | Subcommand | Input | Output |
|---|---|---|---|
| `so2snpstats.py` | `REVEAL snpstats` | one `.so` file | one `.snpstats.tsv` per sample |
| `compare_snpstats.py` | `REVEAL snpcompare` | two or more `.snpstats.tsv` | one wide comparison `.tsv` |

Typical workflow:

```
# step 1: per-sample stats
REVEAL snpstats --so sample1.so --sample-id Dmel1933_SL07
REVEAL snpstats --so sample2.so --sample-id Dmel1940_SL12

# step 2: compare across samples
REVEAL snpcompare \
    --snpstats Dmel1933_SL07.snpstats.tsv Dmel1940_SL12.snpstats.tsv
```

---

## `so2snpstats.py` — per-position SNP stats

Streams through a `.so` file and emits one row per SNP position.

### Output columns

| Column | Description |
|---|---|
| `seqid` | Sequence name |
| `pos` | **1-based** position in the reference |
| `sampleid` | Value of `--sample-id` |
| `cov` | Coverage at this position |
| `refc` | Reference base |
| `A` `T` `C` `G` | Absolute read counts for each base |
| `A_freq` `T_freq` `C_freq` `G_freq` | `count / cov` for each base |
| `major_alt` | Non-reference base with the highest count |
| `major_alt_freq` | `count(major_alt) / cov` |

### Usage

```
REVEAL snpstats --so FILE --sample-id ID [--outfile FILE]
```

| Argument | Default | Description |
|---|---|---|
| `--so` | required | Input `.so` file |
| `--sample-id` | required | Sample identifier written into the output |
| `--outfile` / `-o` | `<sample_id>.snpstats.tsv` | Output TSV path |

If no SNPs are found, an empty TSV with the correct headers is still written.

---

## `compare_snpstats.py` — cross-sample SNP comparison

Loads two or more `.snpstats.tsv` files and joins positions across samples.
One row per `(seqid, pos)`.

### Absent-sample convention

When a sample has no SNP at a given position, its numeric columns (`cov`,
counts, frequencies) are written as **`0`**. `major_alt` is left empty.
This means a SNP at 80 % in sample A that is absent in sample B produces
`freq_range = 0.8` and is correctly flagged as `SNP_FREQ_SHIFT`.

### Flags

| Flag | Condition | Catches |
|---|---|---|
| `SNP_GAIN` | Position present in some samples but absent in others | New or lost SNP |
| `SNP_FLIP` | `major_alt` differs between any two samples that have the SNP | Allele switch (e.g. A→T) |
| `SNP_FREQ_SHIFT` | `max(major_alt_freq) − min(major_alt_freq) >= --freq-shift` | Large frequency change |

Multiple flags are pipe-separated in the `flag` column.

### Per-sample columns (wide format)

For each sample `S` the output contains:

`cov__S`, `A__S`, `T__S`, `C__S`, `G__S`,
`A_freq__S`, `T_freq__S`, `C_freq__S`, `G_freq__S`,
`major_alt__S`, `major_alt_freq__S`

### Usage

```
REVEAL snpcompare --snpstats FILE [FILE ...] [options]
```

| Argument | Default | Description |
|---|---|---|
| `--snpstats` | required (≥2) | `.snpstats.tsv` files from `so2snpstats.py` |
| `--outfile` / `-o` | `comparison.snps.tsv` | Output TSV path |
| `--freq-shift` | `0.8` | Frequency shift threshold for `SNP_FREQ_SHIFT` (flag when max−min ≥ this value) |
| `--flagged-only` | off | Only write rows that have at least one flag |

### Output sort order

Flagged rows first, then by `freq_range` descending, then `seqid` / `pos`.

---

## Notes

- SNPs in the `.so` file are already filtered during `.so` creation by minimum
  count (`--mc-snp`) and minimum frequency (`--mf-snp`) thresholds passed to
  `REVEAL bam2so`. Only positions passing those thresholds appear here.
- `major_alt_freq` close to 1.0 indicates a nearly fixed substitution;
  values around 0.5 suggest a heterozygous or population-segregating SNP.
- For indel analysis see [indel_stats.md](indel_stats.md).
- For copy-number and coverage analysis see [cov_stats.md](cov_stats.md).
