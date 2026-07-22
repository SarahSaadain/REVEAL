# Indel Stats: `so2indelstats` and `indelcompare`

Per-position indel statistics and cross-sample comparison.

---

## Overview

| Script | Subcommand | Input | Output |
|---|---|---|---|
| `so2indelstats.py` | `REVEAL indelstats` | one `.so` file | one `.indelstats.tsv` per sample |
| `compare_indelstats.py` | `REVEAL indelcompare` | two or more `.indelstats.tsv` | one wide comparison `.tsv` |

Typical workflow:

```
# step 1: per-sample stats
REVEAL indelstats --so sample1.so --sample-id Dmel1933_SL07
REVEAL indelstats --so sample2.so --sample-id Dmel1940_SL12

# step 2: compare across samples
REVEAL indelcompare \
    --indelstats Dmel1933_SL07.indelstats.tsv Dmel1940_SL12.indelstats.tsv
```

---

## `so2indelstats.py` — per-position indel stats

Streams through a `.so` file and emits one row per indel event.
Both insertions (`ins`) and deletions (`del`) are included; use the `type`
column to filter.

### Output columns

| Column | Description |
|---|---|
| `seqid` | Sequence name |
| `pos` | **1-based** start position of the indel |
| `sampleid` | Value of `--sample-id` |
| `type` | `ins` or `del` |
| `cov` | Coverage at the position *immediately before* the indel start (matches the `toSeqEntry` convention in `modules.py`) |
| `indel_length` | Length of the insertion or deletion in bp |
| `indel_count` | Number of supporting reads |
| `indel_freq` | `indel_count / cov` |

> **Coverage note:** For both insertions and deletions, coverage is taken from
> `cov[pos − 1]` (0-based), not `cov[pos]`. This matches how indel frequency
> thresholds are applied during `.so` file creation in `SeqBuilder.toSeqEntry`.

### Usage

```
REVEAL indelstats --so FILE --sample-id ID [--outfile FILE]
```

| Argument | Default | Description |
|---|---|---|
| `--so` | required | Input `.so` file |
| `--sample-id` | required | Sample identifier written into the output |
| `--outfile` / `-o` | `<sample_id>.indelstats.tsv` | Output TSV path |

If no indels are found, an empty TSV with the correct headers is still written.

---

## `compare_indelstats.py` — cross-sample indel comparison

Loads two or more `.indelstats.tsv` files and joins events across samples.
One row per `(seqid, pos, type, indel_length)`.

### Join key

`(seqid, pos, type, indel_length)` — different-length events at the same
position are biologically distinct and compared separately. A 3 bp deletion
and a 5 bp deletion at the same position are two independent rows.

### Absent-sample convention

When a sample has no matching indel event, its numeric columns (`cov`,
`indel_count`, `indel_freq`) are written as **`0`** (not NaN). This means an
indel at 70 % frequency in sample A that is absent in sample B produces
`freq_range = 0.7` and is correctly flagged as `INDEL_FREQ_SHIFT`.

### Flags

| Flag | Condition | Catches |
|---|---|---|
| `INDEL_GAIN` | Indel present in some samples but absent in others | New or lost indel |
| `INDEL_FREQ_SHIFT` | `max(indel_freq) − min(indel_freq) >= --freq-shift` | Large frequency change |

Multiple flags are pipe-separated in the `flag` column.

### Per-sample columns (wide format)

For each sample `S` the output contains:

`cov__S`, `indel_count__S`, `indel_freq__S`

### Usage

```
REVEAL indelcompare --indelstats FILE [FILE ...] [options]
```

| Argument | Default | Description |
|---|---|---|
| `--indelstats` | required (≥2) | `.indelstats.tsv` files from `so2indelstats.py` |
| `--outfile` / `-o` | `comparison.indels.tsv` | Output TSV path |
| `--freq-shift` | `0.8` | Frequency shift threshold for `INDEL_FREQ_SHIFT` (flag when max−min ≥ this value) |
| `--flagged-only` | off | Only write rows that have at least one flag |

### Output sort order

Flagged rows first, then by `freq_range` descending, then `seqid` / `pos`.

---

## Notes

- Indels in the `.so` file are already filtered during `.so` creation by
  minimum count (`--mc-indel`) and minimum frequency (`--mf-indel`) thresholds
  passed to `REVEAL bam2so`. Only events passing those thresholds appear here.
- Each `(pos, type, length)` combination represents a distinct indel allele.
  High `indel_freq` (close to 1.0) suggests a fixed or nearly-fixed event;
  low `indel_freq` suggests a rare or heterozygous event.
- For SNP analysis see [snp_stats.md](snp_stats.md).
- For copy-number and coverage analysis see [cov_stats.md](cov_stats.md).
