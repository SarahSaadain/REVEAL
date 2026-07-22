# Coverage Stats: `so2covstats` and `covcompare`

Per-sequence coverage statistics and cross-sample copy-number comparison.

---

## Overview

| Script | Subcommand | Input | Output |
|---|---|---|---|
| `so2covstats.py` | `REVEAL covstats` | one `.so` file | one `.covstats.tsv` per sample |
| `compare_covstats.py` | `REVEAL covcompare` | two or more `.covstats.tsv` | one wide comparison `.tsv` |

Typical workflow:

```
# step 1: per-sample stats
REVEAL covstats --so sample1.so --sample-id Dmel1933_SL07
REVEAL covstats --so sample2.so --sample-id Dmel1940_SL12

# step 2: compare across samples
REVEAL covcompare --stats Dmel1933_SL07.covstats.tsv Dmel1940_SL12.covstats.tsv \
    --cn-fc 1.0 --cn-abs 10 --outfile comparison.tsv
```

---

## `so2covstats.py` — per-sequence stats

Streams through a `.so` file and computes one summary row per sequence.
Memory usage is O(number of sequences), not O(number of positions).

### Output columns

#### Coverage

| Column | Description |
|---|---|
| `seqid` | Sequence name |
| `sampleid` | Value of `--sample-id` |
| `seq_len` | Number of positions in the sequence |
| `median_cov` | Median coverage across all positions. When data is SCG-normalised this is a copy-number proxy (1 = single-copy, 2 = duplicated, etc.) |
| `mad_cov` | Median absolute deviation (MAD) of per-position coverage. Defined as `median(|cov_i − median_cov|)`. Measures how spread-out coverage is around the median while being resistant to outliers (e.g. a single deep-coverage spike does not inflate it the way standard deviation would). A low MAD relative to the median means coverage is flat and uniform; a high MAD means coverage is patchy or uneven. |
| `cv_cov` | Coefficient of variation: `MAD / median_cov`. Scale-independent — a sequence with median 2 and MAD 0.4 has the same `cv_cov` (0.2) as one with median 50 and MAD 10, making it useful for comparing coverage evenness across sequences at very different copy numbers. Values close to 0 indicate flat, uniform coverage; values > 0.5 suggest substantial patchiness. Set to NaN when `median_cov = 0` (no coverage at all). |
| `max_cov` | Peak coverage; useful for spotting sharp spikes |
| `frac_low` | Fraction of positions with coverage < 0.1 — proxy for absent or deleted regions |

**Interpreting MAD and cv_cov together:**

| Scenario | `median_cov` | `mad_cov` | `cv_cov` | Likely biology |
|---|---|---|---|---|
| Single-copy, even coverage | ~1 | small (< 0.2) | < 0.2 | Clean single-copy element |
| Tandem repeat with uniform copies | ~N | small relative to median | < 0.2 | Stable repeat family at copy number N |
| Partially deleted / heterozygous | ~0.5 | moderate | 0.3–0.8 | Element partially present or hemizygous |
| Mosaic or chimeric insertion | any | high | > 0.5 | Uneven coverage from assembly artefacts or structural variation |
| Absent | ~0 | ~0 | NaN | Element not present in this sample |

> MAD is preferred over standard deviation here because TE coverage tracks frequently contain sharp spikes (e.g. from short repeated sub-sequences that attract multi-mapping reads). A single spike can double the standard deviation while barely moving the MAD.

#### SNP summary (aggregate, not per-position)

| Column | Description |
|---|---|
| `n_snps` | Total alt-allele observations across all positions |
| `snp_density` | `n_snps` per 100 bp (length-normalised) |
| `median_alt` | Median alt-allele count across SNP positions |

> For per-position SNP and indel detail, see [snp_stats.md](snp_stats.md).

### Usage

```
REVEAL covstats --so FILE --sample-id ID [--outfile FILE]
```

| Argument | Default | Description |
|---|---|---|
| `--so` | required | Input `.so` file |
| `--sample-id` | required | Sample identifier written into the output |
| `--outfile` / `-o` | `<sample_id>.covstats.tsv` | Output TSV path |

---

## `compare_covstats.py` — cross-sample comparison

Loads two or more `.covstats.tsv` files, pivots to wide format (one row per
sequence, one column block per sample), and flags sequences whose copy number shifts across samples.

### Flags

A sequence is flagged when **any** condition below is met. Multiple flags are pipe-separated in the `flag` column.

| Flag | Condition | Catches |
|---|---|---|
| `CN_FC` | `cn_log2fc >= --cn-fc` | Relative shifts at low copy number (e.g. 1 → 5 = log2FC 2.32) |
| `CN_ABS` | `cn_abs >= --cn-abs` | Large absolute shifts at high copy number (e.g. 50 → 70) |

### Cross-sample summary columns added

| Column | Formula |
|---|---|
| `cn_min` | Lowest `median_cov` across samples |
| `cn_max` | Highest `median_cov` across samples |
| `cn_abs` | `cn_max − cn_min` |
| `cn_log2fc` | `log2(cn_max / cn_min)` — NaN when `cn_min = 0` |

### Per-sample columns (wide format)

For each sample `S` and each metric `M` from the list below, the output contains a column `M__S`:

`median_cov`, `mad_cov`, `cv_cov`, `max_cov`, `frac_low`,
`n_snps`, `snp_density`, `median_alt`

### Usage

```
REVEAL covcompare --stats FILE [FILE ...] [options]
```

| Argument | Default | Description |
|---|---|---|
| `--stats` | required (≥2) | `.covstats.tsv` files from `so2covstats.py` |
| `--outfile` / `-o` | `comparison.tsv` | Output TSV path |
| `--cn-fc` | `2.0` | log2 fold-change threshold for `CN_FC` flag |
| `--cn-abs` | `10` | Absolute copy-number difference threshold for `CN_ABS` flag |
| `--flagged-only` | off | Only write sequences that have at least one flag |

### Output sort order

Flagged sequences first, then by `cn_abs` descending within each group.

---

## Notes

- `median_cov` is the recommended copy-number proxy when the `.so` file was produced from SCG-normalised coverage. See `REVEAL normalize --help`.
- `frac_low < 0.1` typically indicates a sequence present in the sample; `frac_low` close to 1.0 indicates absence or a large deletion.
- `n_snps` / `snp_density` / `median_alt` here are **aggregate** summaries.
  To compare per-position SNP allele frequencies and detect allele flips between samples, use `REVEAL snpstats` + `REVEAL snpcompare`
  (see [snp_stats.md](snp_stats.md)).
