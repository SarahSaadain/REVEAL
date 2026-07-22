#!/usr/bin/env python3
"""
Estimate average coverage statistics for each entry in a .so file.

Average coverage corresponds to copy number when the input has been
normalized to single-copy gene depth (see REVEAL normalize).

Usage
-----
    REVEAL estimate --so sample.norm.so
    REVEAL estimate --so sample.norm.so --outfile sample.estimate.tsv

Authors
-------
    Robert Kofler
    Sarah Saadain
"""

import argparse
import logging
from pathlib import Path

from modules import SeqEntryReader, Writer, NormFactor
from version import __version__

log = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(
        description="Estimate average coverage statistics for each entry in a .so file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--so", required=True, type=Path, metavar="FILE", dest="seqentry",
        help="A sequence overview (.so) file.",
    )
    parser.add_argument(
        "--end-distance", type=int, default=100, metavar="N", dest="enddist",
        help="Distance from sequence ends excluded when estimating coverage (default: 100).",
    )
    parser.add_argument(
        "--exclude-quantile", type=int, default=25, metavar="N", dest="quantile",
        help="Exclude the most extreme coverage quantiles (default: 25).",
    )
    parser.add_argument(
        "--outfile", "-o", type=Path, metavar="FILE", default=None,
        help="Output file; if omitted, output is written to stdout.",
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

    for se in SeqEntryReader(str(args.seqentry)):
        cost = NormFactor.getCovStat(se, args.enddist, args.quantile)
        topr = [se.seqname, str(len(se.cov))]
        form = []
        for c in cost:
            if c is not None:
                form.append(f"{c:.2f}")
            else:
                form.append("na")
        topr.extend(form)
        writer.write("\t".join(topr))


if __name__ == "__main__":
    main()
