# Coverage Binning for Long Sequences

## Problem

For long sequences (chromosomes, genes, large scaffolds), the `.plotable` file contains one row per base position for each of three coverage tracks (`cov`, `ambcov`, `mcov`). A 10 Mbp sequence produces ~30 million rows. Reading and rendering this in R is extremely slow.

## Solution

Coverage can be averaged into fixed-size windows (bins) before writing the `.plotable` file. A bin size of 100 reduces a 30M-row file to 300K rows with no meaningful loss of visual information.

## Usage

### During `.so` → `.plotable` conversion (recommended)

```bash
REVEAL so2plotable --so sample.so --bin-size 100 --outfile sample.plotable
```

The `--bin-size N` option averages coverage over windows of N positions. Each bin is written as a single row at the bin midpoint. SNPs, insertions, and deletions are not affected — they retain exact coordinates.

### During R plotting (for existing `.plotable` files)

```bash
Rscript visualize-plotable.R sample.plotable output.png --bin=100
```

The `--bin=N` flag bins coverage inside R after reading the file. This speeds up rendering but the full file is still read into memory first, so it is slower than pre-binning with `--bin-size`.

## Choosing a bin size

| Sequence length | Recommended `--bin-size` |
|-----------------|--------------------------|
| < 50 kbp        | 1 (no binning)           |
| 50 kbp – 1 Mbp  | 10 – 50                  |
| 1 – 10 Mbp      | 100                      |
| > 10 Mbp        | 500 – 1000               |

## What is binned vs. what is not

| Data type   | Binned? | Notes |
|-------------|---------|-------|
| `cov`       | Yes     | Mean coverage per bin, position = bin midpoint |
| `ambcov`    | Yes     | Mean ambiguous coverage per bin |
| `mcov`      | Yes     | Mean masked coverage per bin |
| SNPs        | No      | Exact base position retained |
| Insertions  | No      | Exact base position retained |
| Deletions   | No      | Start/end positions retained; arc heights (`startcov`, `endcov`) reflect original per-base values, which may differ slightly from the binned average at that position |
