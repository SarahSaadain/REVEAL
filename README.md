# REVEAL

A Python toolkit for converting BAM/SAM alignment files to sequence overview (SO) format, with support for variant detection (SNPs, indels), coverage normalization, and visualization-ready outputs.

## Overview

REVEAL processes short-read alignments to extract coverage, SNP, and indel information for each reference sequence. It is particularly useful for analyzing genomic regions of interest (e.g. TEs, genes, symbionts) and their coverage and variation.

If you want to use REVEAL as part of a automated pipeline, check out [pastForward](https://github.com/SarahSaadain/pastForward).

### Workflow

![REVEAL Pipeline](img/reveal_pipeline.svg)

## Features

- **BAM/SAM Processing**: Converts alignments to sequence overview (SO) format with pysam
- **Variant Detection**: Identifies SNPs and indels with configurable thresholds
- **Coverage Analysis**: Tracks read depth per position with ambiguous mapping support
- **Normalization**: Normalizes coverage using single-copy genes (SCGs) as reference
- **Copy Number Estimation**: Computes average copy number from coverage data
- **Visualization**: Generates tab-delimited output compatible with R/ggplot2
- **Batch Plotting**: Automates running the R visualization script across one or multiple sample folders

## Installation

### Requirements

- Python 3.7+
- pysam
- R + tidyverse (for visualization)
- conda (recommended)
- (Optional) samtools (for BAM/FASTA indexing)

### Install via conda (recommended)

Create or activate a conda environment, then run the installer:

```bash
conda create -n REVEAL
conda activate REVEAL
bash shell/install.sh
```

The installer will:
1. Install all dependencies (pysam, R, tidyverse) into the active environment via conda
2. Copy the scripts into `$CONDA_PREFIX/lib/reveal/`
3. Create the `REVEAL` entry point in `$CONDA_PREFIX/bin/` so the command is available anywhere while the environment is active

After installation:

```bash
REVEAL --help
```

### Manual install

If you prefer not to use the installer:

```bash
pip install pysam
# install R and tidyverse separately, then run scripts directly:
python src/reveal.py --help
```

## Usage

All tools are accessible through the central `REVEAL` command:

```bash
python src/reveal.py --help
```

```
usage: REVEAL [-h] <subcommand> ...

Subcommands:
  bam2so        Convert BAM/SAM → .so (sequence overview)
  normalize     Normalize .so coverage to single-copy genes
  estimate      Estimate per-entry coverage statistics
  so2plotable   Convert .so → R-plottable format
  plot          Render .plotable files to PNG via R
```

Each subcommand accepts `--help` for its full parameter list:

```bash
python src/reveal.py bam2so --help
```

---

### 1. Convert BAM/SAM to Sequence Overview Format

```bash
REVEAL bam2so \
  --infile alignments.bam \
  --fasta reference.fasta \
  --mapqth 5 \
  --mc-snp 5 \
  --mf-snp 0.1 \
  --mc-indel 3 \
  --mf-indel 0.01 \
  --outfile output.so
```

**Parameters:**
- `--infile`: Input BAM or SAM file (required)
- `--fasta`: Reference FASTA file (required)
- `--mapqth`: Mapping quality threshold (default: 5) — reads below this are counted as "ambiguous"
- `--mc-snp`: Minimum SNP count (default: 5)
- `--mf-snp`: Minimum SNP frequency (default: 0.1)
- `--mc-indel`: Minimum indel count (default: 3)
- `--mf-indel`: Minimum indel frequency (default: 0.01)
- `--outfile`: Output SO file; if omitted, prints to stdout
- `--log-level`: Logging verbosity — DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO)

---

### 2. Normalize Coverage by Single-Copy Genes

```bash
REVEAL normalize \
  --so input.so \
  --scg-end _scg \
  --end-distance 100 \
  --exclude-quantile 25 \
  --outfile normalized.so
```

**Parameters:**
- `--so`: Input SO file (required)
- `--scg-end`: Suffix used to identify single-copy genes (default: `_scg`)
- `--end-distance`: Number of bases to exclude from sequence ends during normalization (default: 100)
- `--exclude-quantile`: Exclude this percentage of extreme coverage values from both ends of the distribution (default: 25)
- `--outfile`: Output file; if omitted, prints to stdout
- `--log-level`: Logging verbosity (default: INFO)

---

### 3. Estimate Copy Number

```bash
REVEAL estimate \
  --so normalized.so \
  --end-distance 100 \
  --exclude-quantile 25 \
  --outfile copy_numbers.txt
```

Outputs tab-delimited format: `seqname <tab> length <tab> mean_coverage <tab> min_coverage <tab> max_coverage <tab> ...`

**Parameters:**
- `--so`: Input SO file (required)
- `--end-distance`: Distance from sequence ends to exclude (default: 100)
- `--exclude-quantile`: Quantile threshold for excluding extreme coverage values (default: 25)
- `--outfile`: Output file; if omitted, prints to stdout
- `--log-level`: Logging verbosity (default: INFO)

---

### 4. Convert to Plotable Format

```bash
# Write all sequences to a single file
REVEAL so2plotable \
  --so input.so \
  --sample-id year1933 \
  --seq-ids gypsy,act_scg \
  --outfile sequences.plotable

# Write each sequence to its own file in a directory
REVEAL so2plotable \
  --so input.so \
  --sample-id year1933 \
  --seq-ids ALL \
  --outdir myplotables/ \
  --prefix strain1_
```

Generates visualization-ready tab-delimited output with columns: `seqname`, `sampleid`, `feature` (cov / ambcov / mcov / snp / del / ins), `position`, `value`.

**Parameters:**
- `--so`: Input SO file (required)
- `--sample-id`: Label for this sample, embedded in every output row (default: `x`)
- `--seq-ids`: Comma-separated list of sequence IDs to include, or `ALL` (default: `ALL`)
- `--outfile`: Write all selected sequences into a single plotable file
- `--outdir`: Write one `.plotable` file per sequence into this directory (mutually exclusive with `--outfile`)
- `--prefix`: Filename prefix applied when using `--outdir` (only valid with `--outdir`)
- `--mask-bed`: BED file of regions to mask (0-based coordinates); masked positions are moved from `cov` to `mcov` and excluded from SNP/indel output
- `--mask-ymax`: Mask positions where coverage exceeds this integer value; excess coverage is capped and moved to `mcov`
- `--log-level`: Logging verbosity (default: INFO)

**Notes:**
- Both normalized and non-normalized SO files can be visualized.
- Plotable files from different samples and sequences can be freely concatenated, enabling joint visualization across samples (e.g. `cat sample1/*.plotable sample2/*.plotable > combined.plotable`).

---

### 5. Plot

`REVEAL plot` automates running `visualize-plotable.R` on one or more folders of `.plotable` files.

**Single folder — plot each file independently:**
```bash
REVEAL plot \
  --folder Dmel01_plottable \
  --outdir results/
```

**Multiple folders — merge same-named files across samples, then plot (enables facetting):**
```bash
REVEAL plot \
  --folders Dmel01_plottable Dmel02_plottable Dmel03_plottable \
  --outdir merged_results/
```

**Parameters:**
- `--folder`: Single sample folder; each `.plotable` file is plotted independently (mutually exclusive with `--folders`)
- `--folders`: One or more sample folders; `.plotable` files with matching names are merged before plotting (mutually exclusive with `--folder`)
- `--outdir` / `-o`: Output directory for the generated plots. Required when `--folders` is used; defaults to the source folder when `--folder` is used.
- `--log`: Use a logarithmic y-axis:
  - `--log` (no value): always use log scale
  - `--log 1000` (with a value): automatically switch to log scale if the maximum coverage exceeds the given threshold

**Examples:**
```bash
# Always use log scale
REVEAL plot --folder Dmel01_plottable --outdir results/ --log

# Auto-switch to log scale if max coverage exceeds 1000
REVEAL plot --folders Dmel01_plottable Dmel02_plottable --outdir merged_results/ --log 1000
```

---

### 6. Visualize directly with R

```bash
Rscript src/visualize-plotable.R sequences.plotable output.png
```

Supported output extensions: `.png`, `.pdf`, `.eps`, `.svg`.

Plotable files from different samples can be concatenated and passed directly to the R script to invoke facetting:
```bash
cat sample1.plotable sample2.plotable > combined.plotable
Rscript src/visualize-plotable.R combined.plotable output.png
```

An optional `--log` flag enables a logarithmic y-axis:
```bash
Rscript src/visualize-plotable.R sequences.plotable output.png --log
Rscript src/visualize-plotable.R sequences.plotable output.png --log=1000
```

---

## File Formats

### Sequence Overview (SO) Format

A tab-delimited format containing:
- Coverage per position
- Ambiguous coverage (low MAPQ reads)
- SNP calls with per-base allele counts
- Indel calls (insertions/deletions) with position, length, and count

### Plotable Format

Tab-delimited format for direct plotting in R/ggplot2:
```
seqname  sampleid  feature  position  value
TE_001   sample_1  cov      1         42.5
TE_001   sample_1  snp      100       A  T  3
```

Feature types:
- `cov`: per-position coverage (unmasked)
- `ambcov`: per-position coverage from low-MAPQ (ambiguously mapping) reads
- `mcov`: per-position masked coverage — positions masked via `--mask-bed` or `--mask-ymax`; visualized separately to indicate regions excluded from variant calling
- `snp`: SNP calls (columns: seqname, sampleid, snp, pos, ref_base, alt_base, count)
- `ins`: insertion calls (columns: seqname, sampleid, ins, pos, length, count)
- `del`: deletion calls (columns: seqname, sampleid, del, start, end, start_cov, end_cov, count)

---

## Examples

### Complete Pipeline

```bash
# 1. Convert BAM to SO
REVEAL bam2so \
  --infile reads.bam \
  --fasta te_library.fasta \
  --outfile raw.so

# 2. Normalize coverage
REVEAL normalize \
  --so raw.so \
  --scg-end _scg \
  --outfile normalized.so

# 3. Estimate copy numbers
REVEAL estimate \
  --so normalized.so \
  --outfile copy_numbers.txt

# 4. Format for plotting (one file per sequence)
REVEAL so2plotable \
  --so normalized.so \
  --sample-id my_sample \
  --seq-ids ALL \
  --outdir my_sample_plotables/

# 5. Plot all sequences in the folder
REVEAL plot \
  --folder my_sample_plotables/ \
  --outdir my_sample_plots/
```

### Adjusting Variant Thresholds

```bash
REVEAL bam2so \
  --infile reads.bam \
  --fasta te_library.fasta \
  --mc-snp 10 \
  --mf-snp 0.2 \
  --outfile output.so
```

### Comparing Multiple Samples

```bash
# Generate plotables per sample
for sample in Dmel01 Dmel02 Dmel03; do
  REVEAL so2plotable \
    --so ${sample}_normalized.so \
    --sample-id ${sample} \
    --seq-ids ALL \
    --outdir ${sample}_plotables/
done

# Merge and plot
REVEAL plot \
  --folders Dmel01_plotables Dmel02_plotables Dmel03_plotables \
  --outdir merged_plots/
```

---

## Core Modules

### `modules.py`

Contains shared utilities:
- **SeqEntry**: Data class for sequence overview records
- **SeqBuilder**: Accumulates reads for a sequence and extracts variants
- **SeqEntryReader**: Iterator over SO files (supports gzip)
- **Writer**: Output writer (file or stdout)
- **NormFactor**: Normalization using SCG coverage
- **SNP, Indel**: Variant data classes

### `bam2so.py`

Main conversion pipeline:
1. Loads reference FASTA
2. Iterates alignments (skips unmapped, secondary, supplementary)
3. Builds per-position coverage, SNPs, and indels
4. Writes SO format output; sequences with zero coverage are also written

### `normalize-so.py`

Normalizes coverage using the mean coverage of SCGs across the middle of each sequence (excluding ends and extreme quantiles).

### `estimate-so.py`

Computes per-sequence coverage statistics (useful for copy number estimation when coverage has been normalized to SCGs).

### `so2plotable.py`

Transforms SO entries into a tab-delimited format suitable for visualization. Supports writing all sequences to one file or splitting into per-sequence files in an output directory. Optionally masks coverage peaks via a BED file (`--mask-bed`) or a coverage ceiling (`--mask-ymax`); masked positions are output as a separate `mcov` track and excluded from SNP/indel reporting.

### `run_plotable.py`

Batch runner that automates calling `visualize-plotable.R` on folders of `.plotable` files. Supports single-sample mode (independent plots) and multi-sample mode (merge matching files for facetted plots).

### `REVEAL.py`

Central dispatcher. All subcommands (`bam2so`, `normalize`, `estimate`, `so2plotable`, `plot`) delegate to the corresponding script, so `REVEAL <subcommand> --help` always shows the original script's full help.

---

## Development

### Testing

Unit tests live in `tests/` and cover `modules.py`. Run them with:

```bash
pip install pytest numpy
pytest
```

---

## Tips & Best Practices

1. **FASTA indexing**: Create a `.fai` index for faster access:
   ```bash
   samtools faidx reference.fasta
   ```

2. **BAM indexing**: Ensure BAM files are sorted and indexed:
   ```bash
   samtools sort -o sorted.bam input.bam
   samtools index sorted.bam
   ```

3. **Threshold tuning**: Adjust `--mc-snp`, `--mf-snp`, etc. based on coverage depth and expected variation rates.

4. **Single-copy genes**: If normalizing, ensure your reference includes sequences with the `_scg` suffix (or your chosen `--scg-end`); coverage is normalized to the mean coverage of these sequences.

5. **Large files**: SO format can be space-intensive. Consider gzip compression for archival:
   ```bash
   gzip output.so
   # SO files ending in .gz are read automatically
   REVEAL so2plotable --so output.so.gz ...
   ```

6. **Multi-sample comparison**: Use `REVEAL plot --folders` to merge plotable files across samples automatically rather than concatenating files by hand.

---

## Troubleshooting

### "Reference 'X' not found in FASTA"
Ensure the BAM header references sequences present in your FASTA file:
```bash
samtools view -H alignments.bam | head
samtools faidx reference.fasta && head reference.fasta.fai
```

### Missing FASTA index
```bash
samtools faidx reference.fasta
```

### pysam import errors
```bash
pip install --upgrade pysam
```

### Normalization factor is zero
This usually means no sequences with the `--scg-end` suffix were found, or their coverage is effectively zero. Verify that your reference FASTA includes properly named SCG sequences and that reads mapped to them.

---

## Acknowledgements

REVEAL is based on [teplotter](https://github.com/RobertKofler/teplotter).


## Authors

[Robert Kofler](https://github.com/RobertKofler), [SarahSaadain](https://github.com/SarahSaadain)

## License

MIT License (MIT)

## Citation

If you use REVEAL in your research, please cite:
```
TODO
```
