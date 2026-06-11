#!/usr/bin/env python3
"""Golden-task verification: the engine must reproduce the verified GSE57577 numbers,
unattended, from the manifest. This is the regression gate for the whole roadmap."""
import os, sys, tarfile, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import geo, loaders, expression as ex

DATA = os.path.join(os.path.dirname(__file__), "data", "GSE57577")
os.makedirs(DATA, exist_ok=True)

def ensure_rnaseq():
    tar = geo.download("GSE57575", "GSE57575_RAW.tar", DATA)
    with tarfile.open(tar) as t: t.extractall(DATA)
    for gz in glob.glob(os.path.join(DATA, "*.gz")):
        out = gz[:-3]
        if not os.path.exists(out):
            import gzip, shutil
            with gzip.open(gz) as fi, open(out, "wb") as fo: shutil.copyfileobj(fi, fo)

def main():
    ensure_rnaseq()
    mat, names, loci = loaders.load_cufflinks_fpkm_dir(DATA)
    print(f"[engine] expression matrix {mat.shape[0]} genes x {mat.shape[1]} samples")
    # expression + ncRNA filter (required for a well-behaved eBayes prior)
    keep = ((mat >= 1).sum(axis=1) >= 2)
    keep &= ~pd.Series({i: ex.is_noise(names.get(i, i)) for i in mat.index}).reindex(mat.index).values
    Lall = ex.build_logmatrix(mat)
    L = Lall.loc[keep]
    print(f"[engine] genes tested after filter: {L.shape[0]}")

    checks = []
    def de(a, b):
        r = ex.paired_moderated_de(L, a, b, pairs=[(f"{a}_Set1", f"{b}_Set1"), (f"{a}_Set2", f"{b}_Set2")])
        r["gene"] = [names.get(i, i) for i in r.index]
        return r
    res = {c: de(c, "WT") for c in ["WWD", "R", "TKO"]}
    n = {c: int((res[c].FDR < 0.05).sum()) for c in res}
    print(f"[engine] FDR<0.05  WWD={n['WWD']}  R={n['R']}  TKO={n['TKO']}")

    # golden expectations (from the verified reproduction)
    checks.append(("expr matrix shape", mat.shape == (38227, 8)))
    checks.append(("WWD DE ~1888 (±60)", abs(n["WWD"] - 1888) <= 60))
    checks.append(("R DE small (<40)", n["R"] < 40))
    checks.append(("TKO DE small (<40)", n["TKO"] < 40))
    targets = ["Gata4","Dab2","Lama1","Col4a1","Col4a2","Enc1"]
    sig = res["WWD"][res["WWD"]["gene"].isin(targets) & (res["WWD"].FDR < 0.05)]["gene"].tolist()
    checks.append((f"named targets sig ({len(sig)}/6)", len(sig) >= 5))

    # ChIP ratio check
    try:
        chip = geo.download("GSE57574", "GSE57574_H3K4me3_density.txt.gz", DATA)
        d = pd.read_csv(chip, sep="\t")
        col = {("Dnmt3a2_WWD" in c): c for c in d.columns}
        wwd = d[[c for c in d.columns if "Dnmt3a2_WWD" in c][0]].mean()
        wt  = d[[c for c in d.columns if "Dnmt3a2_WT" in c][0]].mean()
        ratio = wwd / wt
        print(f"[engine] ChIP WWD/WT Dnmt3a2 ratio = {ratio:.2f}")
        checks.append(("ChIP WWD/WT ~1.64 (±0.15)", abs(ratio - 1.64) <= 0.15))
    except Exception as e:
        checks.append((f"ChIP check (skipped: {e})", False))

    print("\n=== GOLDEN TASK ===")
    ok = True
    for name, passed in checks:
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
        ok &= passed
    print(f"\nRESULT: {'ALL PASS ✅' if ok else 'FAILURES ❌'}")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
