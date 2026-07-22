#!/usr/bin/env python3
"""
Convert a BAM/SAM file to SO (sequence observation) format.

Iterates over aligned reads, skipping unmapped, secondary, and supplementary
alignments, and summarises per-position coverage, SNPs, and indels for every
reference sequence.  Reference sequences with no aligned reads are emitted as
empty entries so the output is complete across the whole genome.

Usage
-----
    REVEAL bam2so --infile reads.bam --fasta genome.fa --outfile output.so
    REVEAL bam2so --infile reads.bam --fasta genome.fa --mc-snp 10 --mf-snp 0.05

Authors
-------
    Robert Kofler
    Sarah Saadain
"""

import argparse
import logging
from pathlib import Path

import pysam

from modules import SeqBuilder, Writer, load_fasta
from version import __version__

log = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(
        description="Convert a BAM/SAM file to SO (sequence observation) format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--infile", required=True, type=Path, metavar="FILE",
        help="Input BAM or SAM file path.",
    )
    parser.add_argument(
        "--fasta", required=True, type=Path, metavar="FILE",
        help="The FASTA file to which reads were mapped.",
    )
    parser.add_argument(
        "--mapqth", type=int, default=5, metavar="N",
        help="Mapping quality threshold; reads below this are treated as ambiguous (default: 5).",
    )
    parser.add_argument(
        "--mc-snp", type=int, dest="mcsnp", default=5, metavar="N",
        help="Minimum count of SNPs (default: 5).",
    )
    parser.add_argument(
        "--mf-snp", type=float, dest="mfsnp", default=0.1, metavar="F",
        help="Minimum frequency of SNPs (default: 0.1).",
    )
    parser.add_argument(
        "--mc-indel", type=int, dest="mcindel", default=3, metavar="N",
        help="Minimum count of indels (default: 3).",
    )
    parser.add_argument(
        "--mf-indel", type=float, dest="mfindel", default=0.01, metavar="F",
        help="Minimum frequency of indels (default: 0.01).",
    )
    parser.add_argument(
        "--outfile", "-o", type=Path, metavar="FILE", default=None,
        help="Output file in SO format; if omitted, output is written to stdout.",
    )
    parser.add_argument(
        "--log-level", dest="loglevel", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()
    logging.getLogger().setLevel(args.loglevel)

    #if no output file is provided, don't write log to screen, otherwise it will mess up the output
    if args.outfile is None:
        logging.getLogger().setLevel("ERROR")

    if not args.infile.is_file():
        parser.error(f"File not found: {args.infile}")
    if not args.fasta.is_file():
        parser.error(f"File not found: {args.fasta}")

    writer = Writer(str(args.outfile) if args.outfile else None)
    reference_dict = load_fasta(str(args.fasta))

    builder = None
    seen_sequences = set()

    mode = 'rb' if str(args.infile).lower().endswith('.bam') else 'r'
    log.info("Processing file: %s with mode: %s", args.infile, mode)

    samfile = pysam.AlignmentFile(str(args.infile), mode)

    for read in samfile:
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue

        ref_name = read.reference_name
        pos = read.reference_start  # 0-based
        mapq = read.mapping_quality if read.mapping_quality is not None else 0
        cigar = read.cigarstring
        read_sequence = read.query_sequence.upper() if read.query_sequence is not None else ''

        if cigar is None or cigar == '*':
            continue

        if builder is None:
            builder = SeqBuilder(reference_dict[ref_name], ref_name, args.mapqth)

        if ref_name != builder.seqname:
            seq_entry = builder.toSeqEntry(
                mcsnp=args.mcsnp,
                mfsnp=args.mfsnp,
                mcindel=args.mcindel,
                mfindel=args.mfindel,
            )
            writer.write(str(seq_entry))
            seen_sequences.add(builder.seqname)
            builder = SeqBuilder(reference_dict[ref_name], ref_name, args.mapqth)

        builder.add_read(pos, cigar, mapq, read_sequence)

    samfile.close()

    if builder is not None:
        seq_entry = builder.toSeqEntry(args.mcsnp, args.mfsnp, args.mcindel, args.mfindel)
        writer.write(str(seq_entry))
        seen_sequences.add(seq_entry.seqname)

    for ref_name, ref_sequence in reference_dict.items():
        if ref_name not in seen_sequences:
            empty_builder = SeqBuilder(ref_sequence, ref_name, args.mapqth)
            empty_entry = empty_builder.toSeqEntry(args.mcsnp, args.mfsnp, args.mcindel, args.mfindel)
            writer.write(str(empty_entry))


if __name__ == "__main__":
    main()
