#!/usr/bin/env python3
"""
Normalize coverage in a .so file to single-copy gene (SCG) mean coverage.

Sequences whose names end with --scg-end are used to derive the normalization
factor. After normalization, coverage of 1.0 corresponds to single-copy depth.

Usage
-----
    REVEAL normalize --so sample.so --outfile sample.norm.so
    REVEAL normalize --so sample.so --scg-end _scg --outfile sample.norm.so

Authors
-------
    Robert Kofler
    Sarah Saadain
"""

import argparse
import logging
import sys
from pathlib import Path

from modules import SeqEntryReader, NormFactor, Writer
from version import __version__

log = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(
        description="Normalize coverage in a .so file to single-copy gene mean coverage.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--so", required=True, type=Path, metavar="FILE", dest="seqentry",
        help="A sequence overview (.so) file.",
    )
    parser.add_argument(
        "--scg-end", type=str, default="_scg", metavar="STR", dest="scgend",
        help="Suffix that identifies single-copy-gene entries (default: _scg).",
    )
    parser.add_argument(
        "--end-distance", type=int, default=100, metavar="N", dest="enddist",
        help="Distance from sequence ends excluded before computing the normalization factor (default: 100).",
    )
    parser.add_argument(
        "--exclude-quantile", type=int, default=25, metavar="N", dest="quantile",
        help="Exclude the most extreme coverage quantiles when normalizing (default: 25).",
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

    if not args.seqentry.is_file():
        parser.error(f"File not found: {args.seqentry}")

    writer = Writer(str(args.outfile) if args.outfile else None)
    normfactor = NormFactor.getNormalizationFactor(str(args.seqentry), args.scgend, args.enddist, args.quantile)

    if normfactor == 0:
        log.error(
            "Normalization factor is zero. This may indicate that no single copy genes were found "
            "or that the coverage is zero or very low. Please check your input data."
        )
        sys.exit(1)

    for se in SeqEntryReader(str(args.seqentry)):
        sen = se.normalize(normfactor)
        writer.write(str(sen))


if __name__ == "__main__":
    main()
