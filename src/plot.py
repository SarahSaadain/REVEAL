#!/usr/bin/env python3
"""
Run visualize-plotable.R on .plotable files from one or multiple sample folders.

Usage
-----

Single file mode:
    REVEAL plot --file sample.plotable
    REVEAL plot --file sample.plotable --outdir results/
    REVEAL plot --file sample.plotable --log
    REVEAL plot --file sample.plotable --ymax 500

Multi-file (merged/facet) mode — same seqid across individuals:
    REVEAL plot --files ind1/seq.plotable ind2/seq.plotable --outdir merged_results/
    REVEAL plot --files ind1/seq.plotable ind2/seq.plotable --outdir merged_results/ --log
    REVEAL plot --files ind1/seq.plotable ind2/seq.plotable --outdir merged_results/ --merged-dir kept_merges/
    REVEAL plot --files ind1/seq.plotable ind2/seq.plotable --outdir merged_results/ --ymax 500

Single folder mode:
    REVEAL plot --folder Dmel01_plottable
    REVEAL plot --folder Dmel01_plottable --output results/
    REVEAL plot --folder Dmel01_plottable --log
    REVEAL plot --folder Dmel01_plottable --log 1000
    REVEAL plot --folder Dmel01_plottable --ymax 500

Multi-folder (merged/facet) mode:
    REVEAL plot --folders Dmel01_plottable Dmel02_plottable --output merged_results/
    REVEAL plot --folders Dmel01_plottable Dmel02_plottable --output merged_results/ --log
    REVEAL plot --folders Dmel01_plottable Dmel02_plottable --output merged_results/ --merged-dir kept_merges/
    REVEAL plot --folders Dmel01_plottable Dmel02_plottable --output merged_results/ --ymax 500

Authors
-------
    Robert Kofler
    Sarah Saadain
"""

import argparse
import contextlib
import logging
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import tempfile

log = logging.getLogger(__name__)


def find_plotables(folder: Path) -> dict[str, Path]:
    """Return a dict of {filename: path} for all .plotable files in folder."""
    return {f.name: f for f in folder.glob("*.plotable")}


R_SCRIPT = Path(__file__).parent / "visualize-plotable.R"


def run_rscript(input_path: Path, output_path: Path, extra_args: list[str] | None = None, capture: bool = False):
    """Call the R script with input and output paths, and optional extra flags."""
    cmd = ["Rscript", str(R_SCRIPT), str(input_path), str(output_path)]
    if extra_args:
        cmd.extend(extra_args)
    log.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=capture)
    if result.returncode != 0:
        log.warning("Rscript exited with code %d for %s", result.returncode, input_path.name)
        if capture and result.stderr:
            log.warning("stderr: %s", result.stderr.decode(errors="replace").strip())


def merge_plotables(file_paths: list[Path], dest: Path):
    """Concatenate multiple .plotable files into dest."""
    with open(dest, "wb") as out_fh:
        for i, fp in enumerate(file_paths):
            with open(fp, "rb") as in_fh:
                content = in_fh.read()
                if i > 0 and content and not content.startswith(b"\n"):
                    out_fh.write(b"\n")
                out_fh.write(content)


def single_file_mode(file: Path, output: Path, extra_args: list[str] | None = None):
    """Plot a single .plotable file."""
    output.mkdir(parents=True, exist_ok=True)
    out_path = output / (file.stem + ".png")
    log.info("[%s]", file.name)
    run_rscript(file, out_path, extra_args)
    log.info("Done.")


def multi_file_mode(files: list[Path], output: Path, extra_args: list[str] | None = None,
                    merged_dir: Path | None = None):
    """Merge specific .plotable files (same seqid, one per individual) and plot the merged result."""
    output.mkdir(parents=True, exist_ok=True)

    name = files[0].name
    if any(f.name != name for f in files[1:]):
        log.warning("Input files have different names — using '%s' as output stem.", name)

    log.info("Merging %d file(s) for '%s'", len(files), name)

    if merged_dir is not None:
        merged_dir.mkdir(parents=True, exist_ok=True)
        log.info("Merged .plotable file will be saved to: %s", merged_dir)
        tmp_ctx = contextlib.nullcontext(merged_dir)
    else:
        tmp_ctx = tempfile.TemporaryDirectory()

    with tmp_ctx as tmp:
        merged_file = Path(tmp) / name
        merge_plotables(files, merged_file)
        out_path = output / (Path(name).stem + ".png")
        log.info("Plotting merged file to %s", out_path)
        run_rscript(merged_file, out_path, extra_args)

    log.info("Done.")


def single_folder_mode(folder: Path, output: Path, extra_args: list[str] | None = None, threads: int = 1):
    """Plot each .plotable file in folder independently."""
    output.mkdir(parents=True, exist_ok=True)
    plotables = find_plotables(folder)

    if not plotables:
        log.error("No .plotable files found in %s", folder)
        sys.exit(1)

    log.info("Found %d .plotable file(s) in %s using %d thread(s)", len(plotables), folder, threads)
    capture = threads > 1

    def _plot(item):
        name, src_path = item
        out_path = output / (src_path.stem + ".png")
        log.info("[%s]", name)
        run_rscript(src_path, out_path, extra_args, capture=capture)

    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(_plot, item): item[0] for item in sorted(plotables.items())}
        for fut in as_completed(futures):
            exc = fut.exception()
            if exc:
                log.error("Error plotting %s: %s", futures[fut], exc)

    log.info("Done.")


def multi_folder_mode(folders: list[Path], output: Path, extra_args: list[str] | None = None, threads: int = 1,
                      merged_dir: Path | None = None):
    """Merge same-named .plotable files across folders and plot each merged file."""
    output.mkdir(parents=True, exist_ok=True)

    all_names: dict[str, list[Path]] = {}
    for folder in folders:
        for name, path in find_plotables(folder).items():
            all_names.setdefault(name, []).append(path)

    if not all_names:
        log.error("No .plotable files found in any of the provided folders.")
        sys.exit(1)

    n_folders = len(folders)
    for name, paths in sorted(all_names.items()):
        if len(paths) < n_folders:
            log.warning("'%s' missing from %d folder(s) — will merge available copies only.",
                        name, n_folders - len(paths))

    log.info("Found %d unique .plotable file name(s) across %d folder(s) using %d thread(s).",
             len(all_names), n_folders, threads)
    capture = threads > 1

    if merged_dir is not None:
        merged_dir.mkdir(parents=True, exist_ok=True)
        log.info("Merged .plotable files will be saved to: %s", merged_dir)
        tmp_ctx = contextlib.nullcontext(merged_dir)
    else:
        tmp_ctx = tempfile.TemporaryDirectory()

    with tmp_ctx as tmp:
        tmp_path = Path(tmp)

        def _merge_and_plot(name_paths):
            name, paths = name_paths
            merged_file = tmp_path / name
            log.info("[%s] — merging %d file(s)", name, len(paths))
            merge_plotables(paths, merged_file)
            out_path = output / (Path(name).stem + ".png")
            log.info("[%s] — plotting merged file to %s", name, out_path)
            run_rscript(merged_file, out_path, extra_args, capture=capture)

        with ThreadPoolExecutor(max_workers=threads) as pool:
            futures = {pool.submit(_merge_and_plot, item): item[0] for item in sorted(all_names.items())}
            for fut in as_completed(futures):
                exc = fut.exception()
                if exc:
                    log.error("Error processing %s: %s", futures[fut], exc)

    log.info("Done.")


def build_extra_args(log_value: str | None, ymax: float | None) -> list[str]:
    """Convert Python CLI args to R script flag strings."""
    extra: list[str] = []
    if log_value is not None:
        extra.append("--log" if log_value == "" else f"--log={log_value}")
    if ymax is not None:
        extra.append(f"--ymax={ymax}")
    return extra


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(
        description="Run visualize-plotable.R on .plotable files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--file",
        type=Path,
        metavar="FILE",
        help="Single .plotable file to plot directly.",
    )
    mode.add_argument(
        "--files",
        nargs="+",
        type=Path,
        metavar="FILE",
        help="Multiple .plotable files for the same sequence (one per individual). Merged before plotting.",
    )
    mode.add_argument(
        "--folder",
        type=Path,
        metavar="DIR",
        help="Single sample folder. Each .plotable file is plotted independently.",
    )
    mode.add_argument(
        "--folders",
        nargs="+",
        type=Path,
        metavar="DIR",
        help="Multiple sample folders. Same-named .plotable files are merged before plotting.",
    )

    parser.add_argument(
        "--outdir", "-o",
        type=Path,
        metavar="DIR",
        default=None,
        help=(
            "Output directory for plots. "
            "Required when --folders is used. "
            "Defaults to the source folder when --folder is used."
        ),
    )

    parser.add_argument(
        "--merged-dir",
        type=Path,
        metavar="DIR",
        default=None,
        help=(
            "Directory where merged .plotable files are written (only used with --folders). "
            "If omitted, a temporary directory is used and deleted automatically after plotting."
        ),
    )

    parser.add_argument(
        "--threads", "-j",
        type=int,
        default=1,
        metavar="N",
        help="Number of parallel Rscript processes (default: 1).",
    )

    parser.add_argument(
        "--log",
        nargs="?",        # 0 or 1 argument: bare --log or --log N
        const="",         # value when --log is given with no argument
        default=None,     # value when --log is not given at all
        metavar="N",
        help=(
            "Use logarithmic y-axis. "
            "Without a value (--log): always use log scale. "
            "With a value (--log 1000): auto-switch to log if max coverage exceeds N."
        ),
    )

    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        metavar="N",
        help="Maximum value for the y-axis. If omitted, the axis scales dynamically.",
    )

    args = parser.parse_args()

    extra_args = build_extra_args(args.log, args.ymax)

    # Validate
    if args.files is not None and args.outdir is None:
        parser.error("--outdir is required when using --files")
    if args.folders is not None and args.outdir is None:
        parser.error("--outdir is required when using --folders")

    if args.file is not None:
        if not args.file.is_file():
            parser.error(f"File not found: {args.file}")
        output = args.outdir if args.outdir else args.file.parent
        single_file_mode(args.file, output, extra_args)

    elif args.files is not None:
        for f in args.files:
            if not f.is_file():
                parser.error(f"File not found: {f}")
        multi_file_mode(args.files, args.outdir, extra_args, merged_dir=args.merged_dir)

    elif args.folder is not None:
        if not args.folder.is_dir():
            parser.error(f"Folder not found: {args.folder}")
        output = args.outdir if args.outdir else args.folder
        single_folder_mode(args.folder, output, extra_args, threads=args.threads)

    else:
        for f in args.folders:
            if not f.is_dir():
                parser.error(f"Folder not found: {f}")
        multi_folder_mode(args.folders, args.outdir, extra_args, threads=args.threads,
                          merged_dir=args.merged_dir)


if __name__ == "__main__":
    main()