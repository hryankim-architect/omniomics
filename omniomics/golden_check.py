"""GSE57577 golden task as an importable function — used by the CLI, run_golden.py, and pytest.
Reproduces the verified Noh et al. (2015) numbers from public GEO data, unattended."""
import os, glob, tarfile, gzip, shutil
import numpy as np, pandas as pd
from . import geo, loaders, expression as ex, config

def _default_data_dir():
    return os.environ.get("OMNIOMICS_GOLDEN_DATA",
                          os.path.join(config.repo_dir(), "data", "GSE57577"))

def _ensure_rnaseq(d):
    tar = geo.download("GSE57575", "GSE57575_RAW.tar", d)
    with tarfile.open(tar) as t:
        try: t.extractall(d, filter="data")     # py3.12+ safe extraction; future-proof for 3.14
        except TypeError: t.extractall(d)
    for gz in glob.glob(os.path.join(d, "*.gz")):
        out = gz[:-3]
        if not os.path.exists(out):
            with gzip.open(gz) as fi, open(out, "wb") as fo: shutil.copyfileobj(fi, fo)

def compute(data_dir=None) -> dict:
    """Download + compute the GSE57577 golden metrics. Returns a metrics dict."""
    d = data_dir or _default_data_dir(); os.makedirs(d, exist_ok=True)
    _ensure_rnaseq(d)
    mat, names, _ = loaders.load_cufflinks_fpkm_dir(d)
    keep = ((mat >= 1).sum(axis=1) >= 2)
    keep &= ~pd.Series({i: ex.is_noise(names.get(i, i)) for i in mat.index}).reindex(mat.index).values
    L = ex.build_logmatrix(mat).loc[keep]
    def de(a, b):
        r = ex.paired_moderated_de(L, a, b, pairs=[(f"{a}_Set1", f"{b}_Set1"), (f"{a}_Set2", f"{b}_Set2")])
        r["gene"] = [names.get(i, i) for i in r.index]; return r
    res = {c: de(c, "WT") for c in ["WWD", "R", "TKO"]}
    n = {c: int((res[c].FDR < 0.05).sum()) for c in res}
    targets = ["Gata4", "Dab2", "Lama1", "Col4a1", "Col4a2", "Enc1"]
    sig = res["WWD"][res["WWD"]["gene"].isin(targets) & (res["WWD"].FDR < 0.05)]["gene"].tolist()
    chip = geo.download("GSE57574", "GSE57574_H3K4me3_density.txt.gz", d)
    dd = pd.read_csv(chip, sep="\t")
    wwd = dd[[c for c in dd.columns if "Dnmt3a2_WWD" in c][0]].mean()
    wt  = dd[[c for c in dd.columns if "Dnmt3a2_WT"  in c][0]].mean()
    return {"matrix_shape": tuple(mat.shape), "n_filtered": int(L.shape[0]),
            "de": n, "named_targets_sig": len(sig), "chip_wwd_wt_ratio": float(wwd / wt)}

def evaluate(m) -> list:
    """Return [(check_name, passed), ...] for a metrics dict."""
    return [
        ("expr matrix shape 38227x8",   m["matrix_shape"] == (38227, 8)),
        ("WWD DE ~1888 (±60)",          abs(m["de"]["WWD"] - 1888) <= 60),
        ("R DE small (<40)",            m["de"]["R"] < 40),
        ("TKO DE small (<40)",          m["de"]["TKO"] < 40),
        ("named targets sig >=5/6",     m["named_targets_sig"] >= 5),
        ("ChIP WWD/WT ~1.64 (±0.15)",   abs(m["chip_wwd_wt_ratio"] - 1.64) <= 0.15),
    ]

def run(data_dir=None, verbose=True):
    m = compute(data_dir); checks = evaluate(m); ok = all(p for _, p in checks)
    if verbose:
        print(f"[golden] matrix {m['matrix_shape']} filtered {m['n_filtered']}; "
              f"DE WWD={m['de']['WWD']} R={m['de']['R']} TKO={m['de']['TKO']}; "
              f"targets {m['named_targets_sig']}/6; ChIP {m['chip_wwd_wt_ratio']:.2f}")
        for name, p in checks: print(f"  [{'PASS' if p else 'FAIL'}] {name}")
        print("RESULT:", "ALL PASS ✅" if ok else "FAILURES ❌")
    return ok, m, checks
