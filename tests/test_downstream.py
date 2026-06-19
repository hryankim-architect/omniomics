"""Offline unit tests for the RRBS + ChIP downstream helpers and the mm10 context geometry.
No network, no pipeline output — these guard the parsing/geometry logic that turns nf-core output
into the report tables."""
import gzip
import importlib.util
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import omniomics_mouse as mm
import run_meth_concordance as meth
import run_chip_redistribution as chip


# ---------- mm10 geometry ----------
def test_merge_and_subtract():
    assert mm.merge_intervals([(1, 5), (3, 8), (20, 25), (24, 30)]) == [(1, 8), (20, 30)]
    cgi = mm.merge_intervals([(1000, 2000)])
    flank = mm.merge_intervals([(500, 2500)])
    assert mm.subtract(flank, cgi) == [(500, 1000), (2000, 2500)]


def _ctx():
    return mm.Contexts({
        "promoter": {"chr1": [(0, 100)]},
        "cgi": {"chr1": [(1000, 2000)]},
        "cgi_shore": {"chr1": [(500, 1000), (2000, 2500)]},
    })


def test_classify_point_priority():
    c = _ctx()
    assert c.classify_point("chr1", 50) == "promoter"
    assert c.classify_point("chr1", 1500) == "cgi"
    assert c.classify_point("chr1", 2300) == "cgi_shore"
    assert c.classify_point("chr1", 9999) == "distal"
    assert c.classify_point("chrX", 50) == "distal"  # unknown chrom


def test_classify_interval_and_labelset():
    c = _ctx()
    assert c.classify_interval("chr1", 40, 60) == "promoter"
    assert c.classify_interval("chr1", 1400, 1600) == "cgi"
    assert c.classify_interval("chr1", 100, 200) == "distal"  # touches but no overlap
    assert {"promoter", "cgi"} <= c.label_set("chr1", 90, 1500)


# ---------- RRBS cov parser ----------
def test_parse_cov_line():
    assert meth.parse_cov_line("chr1\t100\t100\t75\t6\t2\n", min_cov=5) == ("chr1", 100, 75.0)
    assert meth.parse_cov_line("chr1\t100\t100\t50\t1\t1\n", min_cov=5) is None  # under-covered
    assert meth.parse_cov_line("garbage line") is None
    chrom, pos, pct = meth.parse_cov_line("chr2\t50\t50\t0\t0\t10\n", min_cov=5)
    assert (chrom, pos, pct) == ("chr2", 50, 0.0)


def test_sample_of_longest_match():
    samples = ["WT_RRBS_1", "WT_RRBS_2", "WT"]
    assert meth.sample_of("/x/WT_RRBS_1.bismark.cov.gz", samples) == "WT_RRBS_1"
    assert meth.sample_of("/x/unrelated.cov.gz", samples) is None


# ---------- ChIP narrowPeak parser ----------
def test_read_peaks_and_merge(tmp_path):
    np = tmp_path / "WT_Dnmt3a2_peaks.narrowPeak"
    np.write_text("chr1\t10\t50\tp1\t100\t.\t5\t9\t8\t20\n"
                  "chr1\t40\t90\tp2\t100\t.\t5\t9\t8\t20\n"
                  "bad\tline\n")
    pk = chip.read_peaks(str(np))
    assert pk == [("chr1", 10, 50), ("chr1", 40, 90)]
    m, st = chip.merged_by_chrom(pk)
    assert m["chr1"] == [(10, 90)]  # overlapping peaks merged
    assert mm.interval_overlaps(m["chr1"], st["chr1"], 80, 85)
    assert not mm.interval_overlaps(m["chr1"], st["chr1"], 200, 300)


# ---------- preflight ----------
def test_preflight_all_pass():
    spec = importlib.util.spec_from_file_location("preflight_check", REPO / "nextflow" / "preflight_check.py")
    pf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pf)
    checks = pf.run()
    failed = [name for name, ok, _ in checks if not ok]
    assert not failed, f"preflight failures: {failed}"
    assert len(checks) >= 15
