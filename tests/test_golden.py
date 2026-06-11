"""Golden-task regression tests (run in CI).

- test_gse57577_reproduction: downloads public GEO data, recomputes, asserts the verified numbers.
- test_recorded_brca_metrics: fast, data-free guard on the committed BRCA result CSVs (catches
  accidental edits to recorded headline numbers; the heavy BRCA pipeline runs on the swarm host).
"""
import pathlib
import pandas as pd
import pytest
from omniomics import golden_check

REPO = pathlib.Path(__file__).resolve().parents[1]


def test_gse57577_reproduction():
    ok, metrics, checks = golden_check.run(verbose=False)
    failed = [name for name, passed in checks if not passed]
    assert ok, f"golden failures: {failed} (metrics={metrics})"


def _read_metric(name, key):
    df = pd.read_csv(REPO / name)
    if "metric" in df.columns and "value" in df.columns:
        return float(dict(zip(df["metric"], df["value"]))[key])
    df = df.set_index(df.columns[0])
    col = "AUROC_mean" if "AUROC_mean" in df.columns else df.columns[-1]
    return float(df[col][key])


@pytest.mark.parametrize("fname,key,lo,hi", [
    ("brca_harmonization_metrics.csv", "PC1_cohort_var_after", 0.0,  0.05),
    ("brca_harmonization_metrics.csv", "xcohort_AUROC_after",  0.84, 0.90),
    ("meth_arm_metrics.csv",           "pct_negative",         70.0, 78.0),
    ("meth_arm_metrics.csv",           "AUROC_RNA",            0.88, 0.93),
    ("dmoi_enhancer_auroc.csv",        "DMOI enhancer (6: +disagree)", 0.90, 0.93),
])
def test_recorded_brca_metrics(fname, key, lo, hi):
    f = REPO / fname
    if not f.exists():
        pytest.skip(f"{fname} not committed")
    v = _read_metric(fname, key)
    assert lo <= v <= hi, f"{key}={v} outside [{lo}, {hi}]"
