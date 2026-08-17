import gzip
from io import StringIO

import numpy as np
import pytest

from modules import (
    Indel,
    NormFactor,
    PlotableFormater,
    SeqBuilder,
    SeqEntry,
    SeqEntryReader,
    SNP,
    Writer,
    isssnp,
    load_bed,
    load_fasta,
)


def test_computeNormalization():
    ses = [SeqEntry("t", [10] * 10, [], [], []), SeqEntry("t", [2] * 10, [], [], [])]
    nf = NormFactor.computeNormFactorForSe(ses, 0, 0)
    assert nf == 6, "test1"

    # median-of-medians: SCG1 median=10, SCG2 median=2, median-of-[10,2]=(10+2)/2=6
    ses = [
        SeqEntry("t", [1, 10, 10, 10, 10, 10, 1], [], [], []),
        SeqEntry("t", [1, 2, 2, 2, 2, 2, 1], [], [], []),
    ]
    nf = NormFactor.computeNormFactorForSe(ses, 0, 0)
    assert nf == 6, "test2"


def test_computeNormalization_insufficient_coverage_reports_zero_percentage():
    # at very low depth, a gene can have >=50% zero-coverage positions purely by
    # chance (Poisson zero-inflation), which drives the median -- and thus the
    # overall factor -- to 0. The error should call this out as insufficient
    # coverage, with how much of the gene has zero depth, not a bare division-by-zero.
    se = SeqEntry("t", [0, 0, 0, 0, 1, 1, 1], [], [], [])
    assert np.median(se.cov) == 0

    with pytest.raises(Exception) as excinfo:
        NormFactor.computeNormFactorForSe([se], 0, 0)

    msg = str(excinfo.value)
    assert "insufficient coverage" in msg.lower()
    assert "57.1%" in msg  # 4 of 7 positions have zero depth


def test_computeNormalization_all_scgs_truly_uncovered_reports_100_percent():
    # if every single-copy gene genuinely has zero reads everywhere, the error
    # should reflect that plainly as 100% zero coverage.
    ses = [SeqEntry("t", [0] * 10, [], [], []), SeqEntry("t", [0] * 5, [], [], [])]

    with pytest.raises(Exception) as excinfo:
        NormFactor.computeNormFactorForSe(ses, 0, 0)

    msg = str(excinfo.value)
    assert "insufficient coverage" in msg.lower()
    assert "100.0%" in msg


def test_computeNormalization_end_distance_trims_before_median():
    # contig edges often carry distorted coverage (assembly/mapping artifacts), so
    # end-distance trimming must be applied before computing the per-gene median.
    se = SeqEntry("t", [0, 0, 1, 1, 1, 1, 0, 0], [], [], [])

    nf_no_trim = NormFactor.computeNormFactorForSe([se], 0, 0)
    assert nf_no_trim == 0.5  # untrimmed median of [0,0,1,1,1,1,0,0]

    nf_trimmed = NormFactor.computeNormFactorForSe([se], 2, 0)
    assert nf_trimmed == 1.0  # after trimming 2 positions off each end: [1,1,1,1]


def test_covstat():
    se = SeqEntry(
        "t",
        [0, 2, 2, 2, 2, 2, 1, 2, 3, 2, 2, 0],
        [99, 5, 5, 5, 6, 4, 5, 5, 5, 5, 5, 99],
        [],
        [],
    )
    cs = NormFactor.getCovStat(se, 1, 10)

    assert cs[0] == 2
    assert cs[1] == 1
    assert cs[2] == 3
    assert cs[3] == 5
    assert cs[4] == 4
    assert cs[5] == 6


def test_covstat_does_not_mutate_se_cov():
    # se.cov/ambcov are numpy arrays once loaded via SeqEntry.parse; slicing them
    # yields a view, so NormFactor._getCovTriplet must sort a copy, not the view,
    # or it silently reorders se.cov in place.
    cov = np.array([0, 2, 2, 2, 2, 2, 1, 2, 3, 2, 2, 0], dtype=np.float64)
    ambcov = np.array([99, 5, 5, 5, 6, 4, 5, 5, 5, 5, 5, 99], dtype=np.float64)
    cov_before = cov.copy()
    ambcov_before = ambcov.copy()
    se = SeqEntry("t", cov, ambcov, [], [])
    NormFactor.getCovStat(se, 1, 10)
    assert (se.cov == cov_before).all()
    assert (se.ambcov == ambcov_before).all()


def test_seqentryreader_need_ambcov_false_skips_ambcov():
    se_in = SeqEntry("chr1", [1.0, 2.0, 3.0], [9.0, 9.0, 9.0], [], [])
    so_text = str(se_in) + "\n"

    reader = SeqEntryReader(StringIO(so_text), need_ambcov=False)
    entries = list(reader)

    assert len(entries) == 1
    se = entries[0]
    assert list(se.cov) == [1.0, 2.0, 3.0]
    assert se.ambcov is None


def test_normalize():
    s = SNP("chr1", 1, "A", 5, 6, 7, 1)
    sn = s.normalize(2.0)
    assert sn.ref == "chr1"
    assert sn.pos == 1
    assert sn.refc == "A"
    assert sn.ac == 2.5
    assert sn.tc == 3
    assert sn.cc == 3.5
    assert sn.gc == 0.5

    id = Indel("chr2", "ins", 5, 2, 11)
    idn = id.normalize(2.0)
    assert idn.ref == "chr2"
    assert idn.type == "ins"
    assert idn.pos == 5
    assert idn.count == 5.5
    assert idn.length == 2

    deli = Indel("chr3", "del", 5, 2, 20)
    de = deli.normalize(5.0)
    assert de.ref == "chr3"
    assert de.type == "del"
    assert de.pos == 5
    assert de.count == 4
    assert de.length == 2

    id = Indel("chr2", "ins", 5, 2, 11)
    deli = Indel("chr3", "del", 5, 2, 20)
    s = SNP("chr1", 1, "A", 5, 6, 7, 1)
    se = SeqEntry("te1", [5, 6, 6, 4, 2], [2, 3, 4, 6, 1], [s], [id, deli])
    sn = se.normalize(2)
    assert sn.cov[0] == 2.5
    assert sn.cov[1] == 3
    assert sn.cov[4] == 1
    assert sn.ambcov[0] == 1
    assert sn.ambcov[1] == 1.5
    assert sn.ambcov[4] == 0.5
    assert sn.ambcov[3] == 3
    assert sn.snplist[0].ac == 2.5
    assert sn.indellist[0].count == 5.5
    assert sn.indellist[1].count == 10


def test_getSNP():
    sb = SeqBuilder("AAATTTCCCGGG", "hans", 5)
    sb.add_read(0, "3M", 5, "AAT")
    sb.add_read(0, "3M", 5, "AAT")
    sb.add_read(0, "3M", 5, "TAT")
    sb.add_read(0, "3M", 5, "TCT")
    se = sb.toSeqEntry(2, 0.1, 2, 0.1)

    assert len(se.snplist) == 2
    assert se.snplist[0].pos == 0
    assert se.snplist[0].ac == 2
    assert se.snplist[0].tc == 2
    assert se.snplist[1].pos == 2
    assert se.snplist[1].ac == 0
    assert se.snplist[1].tc == 4


def test_getInsertion():
    sb = SeqBuilder("AAATTTCCCGGG", "hans", 5)
    # 123456---789012
    # 012345---678901 0-based = (6,3) insertions
    # AAATTT---CCCGGG
    #    TTTAAACCC
    sb.add_read(3, "3M3I3M", 5, "TTTAAACCC")
    sb.add_read(3, "3M3I3M", 5, "TTTAAACCC")
    se = sb.toSeqEntry(2, 0.1, 2, 0.1)

    assert len(se.indellist) == 1
    assert se.indellist[0].pos == 6
    assert se.indellist[0].length == 3
    assert se.indellist[0].count == 2
    assert se.indellist[0].type == "ins"


def test_getDeletion():
    sb = SeqBuilder("AAATTTCCCGGG", "hans", 5)
    # 123456890123
    # 012345678901.  0-based = (6,3) deletion
    # AAATTTCCCGGG
    #    TTT---AAA
    sb.add_read(3, "3M3D3M", 5, "TTTAAA")
    sb.add_read(2, "4M3D3M", 5, "TTTTAAA")
    sb.add_read(3, "3M3D3M", 5, "TTTAAC")
    se = sb.toSeqEntry(2, 0.1, 2, 0.1)

    assert len(se.indellist) == 1
    assert se.indellist[0].pos == 6
    assert se.indellist[0].length == 3
    assert se.indellist[0].count == 3
    assert se.indellist[0].type == "del"


def test_Seq_Builder_add():
    # 012345678901
    # AAATTTCCCGGG
    # AAA
    # TTT
    sb = SeqBuilder("AAATTTCCCGGG", "hans", 5)
    sb.add_read(0, "3M", 4, "ACC")
    sb.add_read(0, "3M", 5, "TGG")

    assert sb.covar[0] == 2
    assert sb.ambcovar[0] == 1
    assert sb.covar[1] == 2
    assert sb.ambcovar[1] == 1
    assert sb.covar[2] == 2
    assert sb.ambcovar[2] == 1
    assert sb.covar[3] == 0
    assert sb.ambcovar[3] == 0
    assert sb.snpar[0]["A"] == 1
    assert sb.snpar[0]["T"] == 1

    # 123456---789012
    # 012345---678901
    # AAATTT---CCCGGG
    #    TTTAAACCC
    sb.add_read(3, "3=3I3X", 5, "TTTAAACCC")
    assert sb.covar[3] == 1
    assert sb.covar[4] == 1
    assert sb.covar[5] == 1
    assert sb.covar[6] == 1
    assert sb.covar[7] == 1
    assert sb.covar[8] == 1
    assert sb.covar[9] == 0
    assert sb.snpar[3]["T"] == 1
    assert sb.snpar[6]["A"] == 0
    assert sb.snpar[6]["C"] == 1
    assert sb.inscol[0] == (6, 3), f"got {sb.inscol[0]}"

    # 123456---789012
    # 012345---678901
    # AAATTTCCCGGG
    #    TTT---AAA
    sb.add_read(3, "3M3D3M", 5, "TTTAAA")
    assert sb.covar[3] == 2
    assert sb.covar[4] == 2
    assert sb.covar[5] == 2
    assert sb.covar[6] == 1
    assert sb.covar[7] == 1
    assert sb.covar[8] == 1
    assert sb.covar[9] == 1
    assert sb.covar[10] == 1
    assert sb.covar[11] == 1
    assert sb.delcol[0] == (6, 3), f"got {sb.delcol[0]}"

    sb.add_read(11, "3M", 5, "TTT")


def test_Seq_Builder_init():
    sb = SeqBuilder("AAATTTCCCGGG", "hans", 5)
    assert sb.seq == "AAATTTCCCGGG", "sequence"
    assert sb.seqname == "hans", "seqname"
    assert sb.minmapq == 5, "minmapq"
    assert len(sb.covar) == 12, "length of covar"
    assert len(sb.ambcovar) == 12, "length of ambcovar"
    assert len(sb.snpar) == 12, "length of snpar"
    assert len(sb.inscol) == 0, "length of inscol"
    assert len(sb.delcol) == 0, "length of delcol"


def test_fasta_loader():
    test_content = """>seq1 some description
ACGTACGT
GCTA
>seq2
NNNNNNNNNN
>seq3 empty sequence

>seq4
ATGCATGCATGC
"""

    result = load_fasta(StringIO(test_content))

    expected = {
        "seq1": "ACGTACGTGCTA",
        "seq2": "NNNNNNNNNN",
        "seq3": "",
        "seq4": "ATGCATGCATGC",
    }

    assert len(result) == 4, f"Expected 4 sequences, got {len(result)}"
    assert "seq3" in result, "Missing empty sequence entry"
    assert result["seq3"] == "", "Empty sequence should be empty string"
    assert result == expected, "Dictionary content doesn't match expected"


def test_convert_to_portable():
    se = SeqEntry("tr1", [], [], [], [])
    se.snplist.append(SNP("t", 100, "A", 2, 3, 4, 0))
    sl = PlotableFormater.prepareSNPForPrint(se, "tamtam", {})

    assert len(sl) == 2
    assert sl[0][3] == "101"  # conversion 100->101 R is 1-based
    assert sl[0][6] == "3"
    assert sl[1][3] == "101"  # conversion 100->101 R is 1-based
    assert sl[1][6] == "4"

    se = SeqEntry("tr1", [], [], [], [])
    se.indellist.append(Indel("t", "ins", 200, 3, 10))
    ins = PlotableFormater.prepareIndelForPrint(se, "tamtam", {})
    assert len(ins) == 1
    # conversion of 200 -> 200 (position is now one position before insertion;
    # instead of 1 position after insertion)
    assert ins[0][3] == "200"

    se = SeqEntry("tr1", [i for i in range(1000, 1400)], [], [], [])
    se.indellist.append(Indel("t", "del", 300, 10, 20))
    dele = PlotableFormater.prepareIndelForPrint(se, "tamtam", {})
    assert len(dele) == 1
    assert dele[0][3] == "300"  # conversion of 300 -> 300 (first coordinate before deletion, 1-based)
    assert dele[0][4] == "311"  # conversion of 310 -> 311 (first coordinate after deletion, 1-based)

    cov = PlotableFormater.prepareCoveragForPrint("hans", [20, 30], "sepp", "cov")
    assert len(cov) == 4
    assert cov[0][3] == "1"
    assert cov[3][3] == "2"


def test_filter_portable():
    se = SeqEntry("tr1", [], [], [], [])
    se.snplist.append(SNP("t", 11, "A", 2, 3, 0, 0))
    se.snplist.append(SNP("t", 12, "A", 2, 3, 0, 0))
    se.snplist.append(SNP("t", 13, "A", 2, 3, 0, 0))
    sl = PlotableFormater.prepareSNPForPrint(se, "tamtam", {12: True})

    assert len(sl) == 2
    assert sl[0][3] == "12"  # 11 +1 (remember conversion from 0-based to 1-based)
    assert sl[1][3] == "14"  # 13+1; hence 12+1 should be missing

    se = SeqEntry("tr1", [], [], [], [])
    se.indellist.append(Indel("t", "ins", 111, 3, 10))
    se.indellist.append(Indel("t", "ins", 112, 3, 10))
    se.indellist.append(Indel("t", "ins", 113, 3, 10))
    ins = PlotableFormater.prepareIndelForPrint(se, "tamtam", {112: True})
    assert len(ins) == 2
    assert ins[0][3] == "111"
    assert ins[1][3] == "113"


def test_load_bed_half_open(tmp_path):
    # BED is 0-based, half-open [start, end): "chr1 10 20" covers positions 10..19, not 20
    bed = tmp_path / "mask.bed"
    bed.write_text("chr1\t10\t20\n")

    result = load_bed(str(bed))

    assert 9 not in result["chr1"]
    assert 10 in result["chr1"]
    assert 19 in result["chr1"]
    assert 20 not in result["chr1"]


def test_load_bed_none_path_returns_empty():
    result = load_bed(None)
    assert result["anything"][0] is False


def test_load_bed_skips_comments_and_headers(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("# comment\ntrack name=x\nbrowser position chr1\nchr1\t0\t2\n")

    result = load_bed(str(bed))

    assert list(result["chr1"].keys()) == [0, 1]


def test_getInsertion_at_position_zero():
    # insertion right at the start of the contig: no reference base precedes it,
    # so coverage lookup must not wrap around to the last base via covar[-1]
    sb = SeqBuilder("AAATTTCCCGGG", "hans", 5)
    sb.add_read(0, "3I3M", 5, "GGGAAA")
    sb.add_read(0, "3I3M", 5, "GGGAAA")
    se = sb.toSeqEntry(2, 0.1, 2, 0.1)

    assert len(se.indellist) == 1
    assert se.indellist[0].pos == 0
    assert se.indellist[0].type == "ins"
    assert se.indellist[0].count == 2


def test_getDeletion_at_position_zero():
    # deletion right at the start of the contig: coverage lookup must use covar[0]
    # (real local depth from other reads), not wrap around to covar[-1] (the
    # unrelated last base of the contig, which is 0 here and would otherwise
    # cause this valid deletion to be silently dropped)
    sb = SeqBuilder("AAATTTCCCGGG", "hans", 5)
    sb.add_read(0, "1M", 5, "A")
    sb.add_read(0, "1M", 5, "A")
    sb.add_read(0, "1M", 5, "A")
    sb.add_read(0, "2D3M", 5, "TTT")
    sb.add_read(0, "2D3M", 5, "TTT")
    se = sb.toSeqEntry(2, 0.1, 2, 0.1)

    assert sb.covar[-1] == 0  # last base of contig is untouched by any read
    assert len(se.indellist) == 1
    assert se.indellist[0].pos == 0
    assert se.indellist[0].type == "del"
    assert se.indellist[0].count == 2


def test_prepareIndelForPrint_deletion_at_position_zero():
    # startpos==0 must read cov[0], not wrap around to the last element of cov
    se = SeqEntry("tr1", [i for i in range(1000, 1400)], [], [], [])
    se.indellist.append(Indel("t", "del", 0, 10, 20))

    dele = PlotableFormater.prepareIndelForPrint(se, "tamtam", {})

    assert len(dele) == 1
    assert dele[0][5] == "1000.0"  # startcov must be cov[0], not cov[-1] (1399.0)


def test_isssnp_zero_coverage_is_false():
    assert isssnp("A", {"A": 0, "T": 0, "C": 0, "G": 0}, 0, 1, 0.1) is False


def test_isssnp_true_when_count_and_freq_thresholds_met():
    hash = {"A": 0, "T": 5, "C": 0, "G": 0}
    assert isssnp("A", hash, 10, 2, 0.1) is True


def test_isssnp_false_below_freq_threshold():
    hash = {"A": 0, "T": 1, "C": 0, "G": 0}
    assert isssnp("A", hash, 100, 1, 0.5) is False


def test_isssnp_ref_base_never_counted():
    # even with high count/freq, the reference base itself is never called a SNP
    hash = {"A": 50, "T": 0, "C": 0, "G": 0}
    assert isssnp("A", hash, 50, 1, 0.1) is False


def test_writer_writes_to_file(tmp_path):
    p = tmp_path / "out.txt"
    w = Writer(str(p))
    w.write("hello")
    w.write("world")
    w.__exit__(None, None, None)

    assert p.read_text() == "hello\nworld\n"


def test_writer_writes_to_stdout_when_no_outfile(capsys):
    w = Writer(None)
    w.write("hello")

    assert capsys.readouterr().out == "hello\n"


def test_writer_usable_as_context_manager(tmp_path):
    p = tmp_path / "out.txt"
    with Writer(str(p)) as w:
        w.write("hello")

    assert p.read_text() == "hello\n"
    assert w.file_handle is None  # closed on exit


def test_getNormalizationFactor_reads_scgs_from_file(tmp_path):
    se1 = SeqEntry("contig1_scg", [10.0] * 10, [1.0] * 10, [], [])
    se2 = SeqEntry("contig2_scg", [2.0] * 10, [1.0] * 10, [], [])
    se3 = SeqEntry("contig3", [100.0] * 10, [1.0] * 10, [], [])  # not a SCG, must be ignored
    content = "\n".join(str(se) for se in (se1, se2, se3)) + "\n"
    p = tmp_path / "in.so"
    p.write_text(content)

    nf = NormFactor.getNormalizationFactor(str(p), "_scg", 0, 0)

    assert nf == 6  # mean of medians 10 and 2


def test_getNormalizationFactor_no_scg_suffix_match_raises(tmp_path):
    se1 = SeqEntry("contig1", [10.0] * 5, [1.0] * 5, [], [])
    p = tmp_path / "in.so"
    p.write_text(str(se1) + "\n")

    with pytest.raises(Exception):
        NormFactor.getNormalizationFactor(str(p), "_scg", 0, 0)


def test_seqentryreader_reads_gzip_file(tmp_path):
    se_in = SeqEntry("chr1", [1.0, 2.0], [3.0, 4.0], [], [])
    p = tmp_path / "in.so.gz"
    with gzip.open(p, "wt") as f:
        f.write(str(se_in) + "\n")

    entries = list(SeqEntryReader(str(p)))

    assert len(entries) == 1
    assert entries[0].seqname == "chr1"
    assert list(entries[0].cov) == [1.0, 2.0]


def test_applyMask_masks_ymax_and_localmask_without_mutating_inputs():
    cov = [5, 10, 3, 20]
    ambcov = [1, 2, 3, 4]
    localmask = {1: True}
    cov_before, ambcov_before, localmask_before = list(cov), list(ambcov), dict(localmask)

    new_cov, new_ambcov, mcov, mask = PlotableFormater.applyMask(cov, ambcov, localmask, ymax=15)

    assert new_cov == [5, 0, 3, 0]
    assert new_ambcov == [1, 0, 3, 0]
    assert mcov == [0, 10, 0, 15]
    assert mask == {1: True, 3: True}  # position 3 newly masked for exceeding ymax
    # inputs must be untouched
    assert cov == cov_before
    assert ambcov == ambcov_before
    assert localmask == localmask_before


def test_prepareForPrint_does_not_mutate_se_or_tomask():
    se = SeqEntry("chr1", [5.0, 25.0], [1.0, 2.0], [], [])
    tomask = {"chr1": {}}

    lines = PlotableFormater.prepareForPrint(se, "s1", tomask, ymax=15, bin_size=1)

    assert se.cov == [5.0, 25.0]
    assert se.ambcov == [1.0, 2.0]
    assert tomask["chr1"] == {}
    assert len(lines) > 0
