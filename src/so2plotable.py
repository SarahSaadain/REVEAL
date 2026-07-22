#!/usr/bin/env python3
"""
Convert a sequence overview (.so) file to R-plottable format.

Usage
-----
    REVEAL so2plotable --so sample.so --outfile sample.plotable
    REVEAL so2plotable --so sample.so --seq-ids chr1,chr2 --outfile sample.plotable
    REVEAL so2plotable --so sample.so --outdir plots/ --sample-id Dmel01

Authors
-------
    Robert Kofler
    Sarah Saadain
"""

import argparse
import logging
import re
from pathlib import Path

from modules import PlotableFormater, SeqEntryReader, Writer, load_bed
from version import __version__

log = logging.getLogger(__name__)

_PADTO = 9


def format_col(topr: list) -> str:
    if _PADTO == 0:
        return "\t".join(topr)
    tf = []
    for i in range(_PADTO):
        if i < len(topr):
            tf.append(str(topr[i]))
        else:
            tf.append("")
    return "\t".join(tf)


def sanitize_filename(filename: str, replacement: str = "_") -> str:
    filename = re.sub(r'[^\w\-.]', replacement, filename)
    filename = re.sub(rf'{re.escape(replacement)}+', replacement, filename)
    return filename.strip(replacement)


def parse_bin_size(s: str):
    """Parse --bin-size value into one of three modes:
    - int: fixed bin size for all sequences
    - ('target', int): bin_size = max(1, seq_len // N) per sequence
    - list of (max_len, bin_size): threshold rules sorted by max_len
    """
    s = s.strip()
    if re.fullmatch(r'\d+', s):
        return int(s)
    m = re.fullmatch(r'target:(\d+)', s, re.IGNORECASE)
    if m:
        n = int(m.group(1))
        if n < 1:
            raise argparse.ArgumentTypeError("target bin count must be >= 1")
        return ('target', n)
    rules = []
    for part in s.split(','):
        part = part.strip()
        if ':' not in part:
            raise argparse.ArgumentTypeError(f"invalid bin-size rule: {part!r} — expected 'len:binsize' or 'default:binsize'")
        key, val = part.split(':', 1)
        key, val = key.strip(), val.strip()
        try:
            bin_val = int(val)
        except ValueError:
            raise argparse.ArgumentTypeError(f"bin size must be an integer, got: {val!r}")
        if bin_val < 1:
            raise argparse.ArgumentTypeError(f"bin size must be >= 1, got: {bin_val}")
        if key.lower() == 'default':
            rules.append((float('inf'), bin_val))
        else:
            try:
                rules.append((int(key), bin_val))
            except ValueError:
                raise argparse.ArgumentTypeError(f"threshold must be an integer or 'default', got: {key!r}")
    if not rules:
        raise argparse.ArgumentTypeError("no valid bin-size rules found")
    if not any(max_len == float('inf') for max_len, _ in rules):
        raise argparse.ArgumentTypeError("threshold mode requires a 'default:N' rule for sequences longer than all thresholds")
    rules.sort(key=lambda x: x[0])
    return rules


def resolve_bin_size(parsed, seq_len: int) -> int:
    """Return the bin size for a sequence of the given length."""
    if isinstance(parsed, int):
        return parsed
    if isinstance(parsed, tuple):  # ('target', N)
        return max(1, seq_len // parsed[1])
    # threshold list of (max_len, bin_size)
    for max_len, bin_size in parsed:
        if seq_len <= max_len:
            return bin_size
    return 1  # unreachable if 'default' rule is present


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(
        description="Convert a sequence overview (.so) file to R-plottable format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--so", required=True, type=Path, metavar="FILE", dest="so",
        help="A sequence overview (.so) file.",
    )
    parser.add_argument(
        "--seq-ids", type=str, default="ALL", metavar="IDS", dest="seqids",
        help="Comma-separated sequence IDs to plot, or 'ALL' (default: ALL).",
    )
    parser.add_argument(
        "--sample-id", type=str, default="x", metavar="ID", dest="sampleid",
        help="Sample identifier written into the output (default: x).",
    )
    parser.add_argument(
        "--prefix", type=str, default="", metavar="STR", dest="prefix",
        help="Filename prefix for output files; only valid with --outdir.",
    )
    parser.add_argument(
        "--outdir", type=Path, metavar="DIR", default=None, dest="outputdir",
        help="Output directory; one .plotable file is written per FASTA entry.",
    )
    parser.add_argument(
        "--outfile", "-o", type=Path, metavar="FILE", default=None, dest="outfile",
        help="Output file in plotable format.",
    )
    parser.add_argument(
        "--mask-bed", type=Path, metavar="FILE", default=None, dest="maskbed",
        help="BED file for masking; regions in the file are set to zero (0-based).",
    )
    parser.add_argument(
        "--mask-ymax", type=int, metavar="N", default=None, dest="ymax",
        help="Mask positions with coverage exceeding this value.",
    )
    parser.add_argument(
        "--bin-size", type=str, default="1", metavar="SPEC", dest="binsize",
        help=(
            "Bin size specification. Three modes: "
            "(1) fixed integer: '--bin-size 100' applies the same bin size to every sequence; "
            "(2) target bins: '--bin-size target:1000' auto-computes bin_size = max(1, seq_len // 1000) per sequence; "
            "(3) length thresholds: '--bin-size 10000:1,100000:10,1000000:100,default:500' uses bin size 1 for sequences <= 10000 bp, "
            "10 for <= 100000 bp, and 100 for anything longer. "
            "Default: 1 (no binning)."
        ),
    )
    parser.add_argument(
        "--log-level", dest="loglevel", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()
    logging.getLogger().setLevel(args.loglevel)

    if args.outfile is not None and args.outputdir is not None:
        parser.error("provide either --outfile or --outdir, not both")
    if args.prefix != "" and args.outputdir is None:
        parser.error("--prefix only works when --outdir is provided")

    #if no output file is provided, don't write log to screen, otherwise it will mess up the output
    if args.outfile is None:
        logging.getLogger().setLevel("ERROR")

    if not args.so.is_file():
        parser.error(f"File not found: {args.so}")

    try:
        binspec = parse_bin_size(args.binsize)
    except argparse.ArgumentTypeError as e:
        parser.error(str(e))

    writer = Writer(str(args.outfile) if args.outfile else None)
    tomask = load_bed(str(args.maskbed) if args.maskbed else None)

    if "," in args.seqids:
        seqset = set(args.seqids.split(","))
    else:
        seqset = {args.seqids}

    is_print_all_requested = (args.seqids.upper() == "ALL")

    if args.outputdir is not None:
        args.outputdir.mkdir(parents=True, exist_ok=True)

    for se in SeqEntryReader(str(args.so)):
        if is_print_all_requested or se.seqname in seqset:
            binsize = resolve_bin_size(binspec, len(se.cov))
            log.debug("sequence %s: length=%d, bin_size=%d", se.seqname, len(se.cov), binsize)
            tmp = PlotableFormater.prepareForPrint(se, args.sampleid, tomask, args.ymax, binsize)
            tr = [format_col(i) for i in tmp]
            tp = "\n".join(tr)

            if args.outputdir is not None:
                filename = args.prefix + sanitize_filename(se.seqname) + ".plotable"
                with open(args.outputdir / filename, "w") as f:
                    f.write(tp)
            else:
                writer.write(tp)


if __name__ == "__main__":
    main()
