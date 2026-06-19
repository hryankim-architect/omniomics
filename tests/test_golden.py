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
    ("dmoi_v2_auroc.csv",              "orig DMOI v1 (LR)",            0.905, 0.922),
    ("dmoi_v2_auroc.csv",              "DMOI v2 no-reliability (GBT)", 0.915, 0.930),
    ("dmoi_v2_auroc.csv",              "RNA poles only (LR)",          0.900, 0.915),
    ("dmoi_v2_auroc.csv",              "DMOI v2 gene-level (GBT)",     0.910, 0.925),
    ("dmoi_v2_auroc.csv",              "DMOI v2 (permuted-label null)", 0.44, 0.56),
])
def test_recorded_brca_metrics(fname, key, lo, hi):
    f = REPO / fname
    if not f.exists():
        pytest.skip(f"{fname} not committed")
    v = _read_metric(fname, key)
    assert lo <= v <= hi, f"{key}={v} outside [{lo}, {hi}]"


def test_dmoi_v2_gain_guard():
    """DMOI v2's nonlinear-interaction variant must beat v1 by a clear, significant margin -- the
    recorded headline of reports/DMOI_method_assessment.md (v1's linear 'disagreement' term is
    inert; the cross-omics signal needs interactions + a nonlinear learner). Skip-safe until
    run_dmoi_v2.py has produced the CSV."""
    f = REPO / "dmoi_v2_auroc.csv"
    if not f.exists():
        pytest.skip("dmoi_v2_auroc.csv not committed (run_dmoi_v2.py not run yet)")
    df = pd.read_csv(f).set_index("model")
    v1 = float(df.loc["orig DMOI v1 (LR)", "AUROC_mean"])
    v2 = float(df.loc["DMOI v2 no-reliability (GBT)", "AUROC_mean"])
    p = float(df.loc["DMOI v2 no-reliability (GBT)", "wilcoxon_p_vs_v1"])
    assert v2 - v1 > 0.005, f"v2 interaction gain too small: v2={v2:.4f} v1={v1:.4f}"
    assert p < 0.05, f"v2 gain not significant vs v1: p={p}"
    # permuted-label null must collapse to chance -> the gain is not a pipeline artefact
    if "DMOI v2 (permuted-label null)" in df.index:
        null = float(df.loc["DMOI v2 (permuted-label null)", "AUROC_mean"])
        assert 0.45 <= null <= 0.55, f"permuted-label null not ~0.5: {null}"


def test_auto_integrate_runner_guard():
    """run_auto_integrate.py must pick the dominant modality per task and never fall below it.
    Skip-safe until the runner has produced the CSV."""
    f = REPO / "auto_integrate_results.csv"
    if not f.exists():
        pytest.skip("auto_integrate_results.csv not committed (run_auto_integrate.py not run yet)")
    df = pd.read_csv(f).set_index("endpoint")
    lum = df.loc["LumA_vs_LumB"]
    assert str(lum["anchor"]) == "RNA", f"LumA/B anchor should be RNA, got {lum['anchor']}"
    assert 0.92 <= float(lum["auroc_combined"]) <= 0.96
    assert float(lum["auroc_combined"]) >= float(lum["auroc_anchor"]) - 0.005   # never below the anchor
    pc = df.loc["methylation_defined_posctrl"]
    assert str(pc["anchor"]).startswith("methylation"), f"pos-ctrl anchor should be methylation, got {pc['anchor']}"
    assert float(pc["auroc_combined"]) >= float(pc["auroc_anchor"]) - 0.005


def test_external_subtype_anchor_guard():
    """External validation on REAL TCGA-BRCA subtype labels: auto_integrate must route to the
    modality that defines each task -- methylation for the methylation-defined clusters, RNA for the
    expression-defined PAM50 calls -- and never fall below that leader. Skip-safe until the CSV exists."""
    f = REPO / "external_subtype_results.csv"
    if not f.exists():
        pytest.skip("external_subtype_results.csv not committed (dmoi_external_subtype.py not run yet)")
    df = pd.read_csv(f)
    meth = df[df["kind"] == "methylation-defined"]
    expr = df[df["kind"] == "expression-defined"]
    assert len(meth) >= 3 and (meth["anchor"] == "methylation").all(), \
        f"methylation-defined endpoints should anchor on methylation: {meth[['endpoint','anchor']].values.tolist()}"
    assert len(expr) >= 3 and (expr["anchor"] == "RNA").all(), \
        f"expression-defined endpoints should anchor on RNA: {expr[['endpoint','anchor']].values.tolist()}"
    assert (df["auroc_combined"] >= df["auroc_anchor"] - 0.01).all()   # never below the chosen anchor


def test_fusion_gain_guard():
    """The curated-biomarker fusion-gain result: on NORMAL-tissue age the Horvath clock methylation
    must beat RNA and be the chosen anchor, while a matched random-CpG set must not; on tumour RNA wins.
    Skip-safe until reports/dmoi_fusion_gain.py has produced the CSV."""
    f = REPO / "fusion_gain_results.csv"
    if not f.exists():
        pytest.skip("fusion_gain_results.csv not committed (dmoi_fusion_gain.py not run yet)")
    df = pd.read_csv(f).set_index("tissue")
    nt = df.loc["normal_tissue"]
    assert float(nt["auroc_clock_score"]) > float(nt["auroc_rna"]), "clock should beat RNA on normal-tissue age"
    assert float(nt["clock_minus_rna"]) > 0 and float(nt["clock_gt_rna_frac"]) >= 0.8
    assert str(nt["auto_anchor"]).startswith("clock"), "anchor should be the clock methylation"
    assert float(nt["auroc_clock_score"]) - float(nt["auroc_random"]) > 0.15, "curation must beat random CpGs"
    assert 0.42 <= float(nt["perm_null"]) <= 0.58, "permuted-label null should be ~0.5"
    if "tumor" in df.index:
        tm = df.loc["tumor"]
        assert float(tm["auroc_rna"]) > float(tm["auroc_clock_score"]) and str(tm["auto_anchor"]) == "RNA"


def test_immune_axis_guard():
    """Immune-axis fusion test (honest negative): no super-additive gain on either external label —
    histology is RNA-anchored with the gate adding ~nothing, node status is near-chance for all
    modalities. Skip-safe until reports/dmoi_immune_fusion.py has produced the CSV."""
    f = REPO / "immune_axis_results.csv"
    if not f.exists():
        pytest.skip("immune_axis_results.csv not committed (dmoi_immune_fusion.py not run yet)")
    df = pd.read_csv(f).set_index("endpoint")
    h = df.loc["histology_idc_ilc"]
    assert str(h["auto_anchor"]) == "RNA" and float(h["auto_delta"]) <= 0.01   # no fusion gain
    assert float(h["auroc_rna"]) >= float(h["auroc_meth_genomewide"])           # RNA leads
    nd = df.loc["node_status"]
    assert float(nd[["auroc_rna", "auroc_meth_genomewide", "auroc_immune_deconv", "auroc_immune_cpgs"]].max()) < 0.65
    assert float(nd["auto_delta"]) <= 0.01


def test_discordance_signal_guard():
    """Three 'is cross-omics disagreement signal or noise?' tests on gene-matched LumA/B: the linear
    difference is provably inert, the interaction does not clear the pairing-permutation null, and the
    layers are redundant (not synergistic). Skip-safe until dmoi_discordance_tests.py produces the CSV."""
    f = REPO / "discordance_test_results.csv"
    if not f.exists():
        pytest.skip("discordance_test_results.csv not committed (dmoi_discordance_tests.py not run yet)")
    df = pd.read_csv(f)
    v = {(t, m): float(val) for t, m, val in df[["test", "metric", "value"]].itertuples(index=False)}
    assert abs(v[("pairing_permutation", "linear_diff_increment")]) <= 0.005   # linear diff provably inert
    assert v[("pairing_permutation", "increment_int")] <= 0.02                 # interaction ~ does not clear null
    assert v[("synergy", "info_synergy_bits")] <= 0.05                         # redundant, not synergistic
    assert v[("synergy", "auroc_joint_nonlinear")] <= v[("synergy", "auroc_rna")] + 0.01  # joint not above RNA alone


def test_modern_de_concordance_guard():
    """If the nf-core/DESeq2 reanalysis has been run, its '2015 vs 2026' direction must hold;
    otherwise this is a no-op (heavy run happens on the Linux node)."""
    checks = golden_check.modern_de_concordance()
    if checks is None:
        pytest.skip("modern_de_concordance.csv not committed (nf-core run not done yet)")
    failed = [name for name, ok in checks if not ok]
    assert not failed, f"modern-DE direction failures: {failed}"


def test_modern_de_parser_logic(tmp_path):
    """Data-free unit test of the concordance parser/direction logic (runs in CI)."""
    csv = ("contrast,DESeq2_padj<0.05,our_n2_moderated,named_targets_sig(/6)\n"
           "WWD_vs_WT,1503,1888,5\nR_vs_WT,14,9,\nTKO_vs_WT,6,3,\n")
    (tmp_path / "modern_de_concordance.csv").write_text(csv)
    checks = golden_check.modern_de_concordance(repo=str(tmp_path))
    assert checks is not None and all(ok for _, ok in checks), checks
    # and the failing direction (WWD not dominant) must be caught
    bad = ("contrast,DESeq2_padj<0.05,our_n2_moderated,named_targets_sig(/6)\n"
           "WWD_vs_WT,5,1888,0\nR_vs_WT,400,9,\nTKO_vs_WT,6,3,\n")
    (tmp_path / "modern_de_concordance.csv").write_text(bad)
    checks = golden_check.modern_de_concordance(repo=str(tmp_path))
    assert not all(ok for _, ok in checks)
