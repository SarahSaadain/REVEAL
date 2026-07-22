from io import StringIO

import numpy as np

from modules import (
    Indel,
    NormFactor,
    PlotableFormater,
    SeqBuilder,
    SeqEntry,
    SeqEntryReader,
    SNP,
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
