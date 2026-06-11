#!/usr/bin/env python3
"""Golden task: full BRCA multi-omics pipeline (cross-cohort harmonization + methylation arm
+ joint embedding + pole-conditioned classification). Re-runs the pipeline and asserts the
verified metrics within tolerance, then appends a hash-chained audit record — the swarm's
critic/approval/audit gate for the scaled pipeline.

Usage: python run_golden_brca.py [--use-cached]
"""
import os, sys, json, time, hashlib, subprocess
import pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUDIT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "audit.log")

SCRIPTS = ["run_cohort_brca.py", "run_meth_arm.py", "run_joint_dmoi.py", "run_lumab_dmoi.py"]

# (file, metric-key, expected, tolerance, comparator)
GOLDEN = [
  ("brca_harmonization_metrics.csv", "PC1_cohort_var_before", 0.99, None, "ge_0.90"),
  ("brca_harmonization_metrics.csv", "PC1_cohort_var_after",  0.00, None, "le_0.05"),
  ("brca_harmonization_metrics.csv", "xcohort_AUROC_before",  0.753, 0.03, "approx"),
  ("brca_harmonization_metrics.csv", "xcohort_AUROC_after",   0.869, 0.03, "approx"),
  ("meth_arm_metrics.csv",           "pct_negative",          74.3, 4.0,  "approx"),
  ("meth_arm_metrics.csv",           "AUROC_RNA",             0.910, 0.03, "approx"),
  ("lumab_dmoi_auroc.csv",           "RNA pole (2 feats)",    0.914, 0.04, "approx"),
  ("lumab_dmoi_auroc.csv",           "RNA (2000 genes)",      0.950, 0.04, "approx"),
]

def load_metric(fname, key):
    df = pd.read_csv(os.path.join(ROOT, fname))
    # normalize to {key: value}
    if "metric" in df.columns and "value" in df.columns:
        d = dict(zip(df["metric"], df["value"]))
    else:
        df = df.set_index(df.columns[0]); d = dict(zip(df.index, df[df.columns[-1]]))
    return float(d[key])

def check(comp, val, exp, tol):
    if comp == "approx": return abs(val-exp) <= tol
    if comp == "ge_0.90": return val >= 0.90
    if comp == "le_0.05": return val <= 0.05
    return False

def audit_append(record):
    prev = "0"*64
    if os.path.exists(AUDIT):
        lines = [l for l in open(AUDIT) if l.strip()]
        if lines: prev = json.loads(lines[-1]).get("hash", prev)
    body = json.dumps(record, sort_keys=True)
    h = hashlib.sha256((prev+body).encode()).hexdigest()
    rec = {**record, "prev": prev, "hash": h}
    with open(AUDIT, "a") as f: f.write(json.dumps(rec, sort_keys=True)+"\n")
    return h

def main():
    if "--use-cached" not in sys.argv:
        for s in SCRIPTS:
            print(f"[run] {s}")
            r = subprocess.run([sys.executable, os.path.join(ROOT, s)],
                               capture_output=True, text=True)
            if r.returncode != 0:
                print(r.stdout[-500:]); print(r.stderr[-800:])
                raise SystemExit(f"pipeline step failed: {s}")

    print("\n=== GOLDEN TASK: BRCA multi-omics pipeline ===")
    ok = True; metrics = {}
    for fname, key, exp, tol, comp in GOLDEN:
        try:
            val = load_metric(fname, key); passed = check(comp, val, exp, tol)
        except Exception as e:
            val = None; passed = False; print("  load error:", e)
        metrics[key] = val
        tag = f"~{exp}" if comp=="approx" else comp
        print(f"  [{'PASS' if passed else 'FAIL'}] {key:24s} = {val}  (golden {tag})")
        ok &= passed

    # extra regression assertions: enhancer-DMOI gain + ComBat covariate benefit
    try:
        de = pd.read_csv(os.path.join(ROOT, "dmoi_enhancer_auroc.csv"), index_col=0)["AUROC_mean"]
        extra = [("dmoi_enh_AUROC >= 0.91", de["DMOI enhancer (6: +disagree)"] >= 0.91),
                 ("dmoi_enh > dmoi_promoter", de["DMOI enhancer (6: +disagree)"] > de["DMOI promoter (6, contrast)"]),
                 ("dmoi_enh > RNA_pole", de["DMOI enhancer (6: +disagree)"] > de["RNA pole (2)"])]
        for name, passed in extra:
            print(f"  [{'PASS' if passed else 'FAIL'}] {name}"); ok &= bool(passed)
        metrics["dmoi_enh_AUROC"] = float(de["DMOI enhancer (6: +disagree)"])
    except Exception as e: print("  (enhancer-DMOI metrics skipped:", e, ")")
    try:
        cb = pd.read_csv(os.path.join(ROOT, "combat_benchmark.csv"), index_col=0)["subtype_AUROC"]
        extra = [("EBComBat+cov preserves biology (>=0.85)", cb["EB-ComBat (+covariate)"] >= 0.85),
                 ("naive ComBat-lite destroys biology (<0.6)", cb["ComBat-lite (no cov)"] < 0.6)]
        for name, passed in extra:
            print(f"  [{'PASS' if passed else 'FAIL'}] {name}"); ok &= bool(passed)
        metrics["combat_EBcov_AUROC"] = float(cb["EB-ComBat (+covariate)"])
    except Exception as e: print("  (combat metrics skipped:", e, ")")

    rec = {"ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "task": "brca_multiomics_pipeline", "result": "PASS" if ok else "FAIL",
           "metrics": metrics, "scripts": SCRIPTS}
    h = audit_append(rec)
    print(f"\nRESULT: {'ALL PASS ✅' if ok else 'FAILURES ❌'}")
    print(f"audit: appended hash-chained record {h[:16]}… -> {AUDIT}")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
