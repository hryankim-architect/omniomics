#!/usr/bin/env python3
"""Reproduce the Noh et al. (2015) graphical-abstract pattern with DMOI, from the full GSE57577 data.

Self-contained: downloads its own ChIP (GSE57574 H3K4me1/2/3 density), RRBS (GSE57576 methylcall),
and RNA-seq (GSE57575 FPKM) into ./data/gse57577/ and runs anywhere.

Run locally:
    pip install pandas numpy matplotlib        # one-time
    python3 run_gse57577_dmoi.py
Outputs: figures/gse57577_dmoi_pattern.png + gse57577_localization.csv / gse57577_ctx_methylation.csv
"""
import os, sys, glob, tarfile, gzip, shutil
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from omniomics import loaders, geo
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data", "gse57577")
FIG  = os.path.join(BASE, "figures"); os.makedirs(FIG, exist_ok=True)
VAR  = ["WT", "WWD", "R"]
CTX  = {"H3K4me1 (enhancer)":    "GSE57574_H3K4me1_density.txt.gz",
        "H3K4me2 (gene body)":   "GSE57574_H3K4me2_density.txt.gz",
        "H3K4me3 (CGI/promoter)":"GSE57574_H3K4me3_density.txt.gz"}
def _wrap(c): return c.replace(" (", "\n(")   # two-line labels for the figure only

def ensure_data():
    """Download all three assays from GEO if not already present (idempotent)."""
    cdir = os.path.join(DATA, "chip"); rdir = os.path.join(DATA, "rrbs"); ndir = os.path.join(DATA, "rna")
    for f in CTX.values():
        geo.download("GSE57574", f, cdir)
    for v in VAR:
        geo.download("GSE57576", f"GSE57576_methylcall.Sample_{v}.mincov10.txt.gz", rdir)
    if not glob.glob(os.path.join(ndir, "*.fpkm_tracking")):
        tar = geo.download("GSE57575", "GSE57575_RAW.tar", ndir)
        with tarfile.open(tar) as t: t.extractall(ndir)
        for gz in glob.glob(os.path.join(ndir, "*.gz")):
            if not os.path.exists(gz[:-3]):
                with gzip.open(gz) as fi, open(gz[:-3], "wb") as fo: shutil.copyfileobj(fi, fo)
    return cdir, rdir, ndir

print("[dmoi] ensuring GSE57577 data (ChIP + RRBS + RNA) ...")
CHIPDIR, RRBSDIR, RNADIR = ensure_data()

# ---- 1. Dnmt3a2 localization per context per variant (ChIP) ----
bind = {}; windows = {}
for ctx, fname in CTX.items():
    df = pd.read_csv(os.path.join(CHIPDIR, fname), sep="\t")     # pandas reads .gz directly
    ren = {}
    for c in df.columns:
        for v in VAR:
            if c.startswith(f"Dnmt3a2_{v} "): ren[c] = f"D_{v}"
    df = df.rename(columns=ren)
    bind[ctx] = {v: df[f"D_{v}"].mean() for v in VAR}
    df["chr"] = df["Chr"]; df["s"] = df["Start"]; df["e"] = df["End"]
    windows[ctx] = {c: sub[["s","e"]].sort_values("s").values for c, sub in df.groupby("chr")}
L = pd.DataFrame(bind).T[VAR]
print("\n[dmoi] Dnmt3a2 localization (mean normalized binding):"); print(L.round(3).to_string())
print("\n[dmoi] redistribution vs WT (variant/WT):"); print((L.div(L["WT"], axis=0)).round(2).to_string())

# ---- 2. methylation per context per variant (RRBS, vectorized overlap) ----
def ctx_meth(variant):
    d = pd.read_csv(os.path.join(RRBSDIR, f"GSE57576_methylcall.Sample_{variant}.mincov10.txt.gz"),
                    sep="\t", usecols=["chr","base","freqC"])
    out = {}
    for ctx, wins in windows.items():
        tot = 0.0; n = 0
        for chrom, sub in d.groupby("chr"):
            arr = wins.get(chrom)
            if arr is None: continue
            pos = sub["base"].values; fq = sub["freqC"].values; st = arr[:,0]; en = arr[:,1]
            idx = np.searchsorted(st, pos, side="right") - 1
            ok = (idx >= 0) & (idx < len(st)); hit = ok.copy(); hit[ok] = pos[ok] <= en[idx[ok]]
            tot += fq[hit].sum(); n += hit.sum()
        out[ctx] = tot / n if n else np.nan
    return out
M = pd.DataFrame({v: ctx_meth(v) for v in VAR})[VAR]
print("\n[dmoi] DNA methylation at each context (%):"); print(M.round(2).to_string())

# ---- 3. RNA differentiation-program score per variant ----
mat, names, _ = loaders.load_cufflinks_fpkm_dir(RNADIR)
gm = pd.DataFrame({v: mat[[c for c in mat.columns if c.startswith(v+"_")]].mean(axis=1) for v in VAR})
sym = pd.Series({i: str(names.get(i, i)) for i in gm.index})
DEV = ["Gata4","Gata6","Dab2","Lama1","Col4a1","Col4a2","Sox17","Foxa2","Sox7","Enc1","Hnf4a"]
dev_idx = sym[sym.isin(DEV)].index
diff_score = {v: np.log2(gm.loc[dev_idx, v] + 1).mean() for v in VAR}
print("\n[dmoi] differentiation-program score (mean log2 FPKM of endoderm/lineage genes):")
for v in VAR: print(f"   {v}: {diff_score[v]:.2f}")

# ---- 4. DMOI binding-vs-methylation disagreement ----
zL = (L.sub(L.mean(axis=1), axis=0)).div(L.std(axis=1) + 1e-9, axis=0)
zM = (M.sub(M.mean(axis=1), axis=0)).div(M.std(axis=1) + 1e-9, axis=0)
DIS = zL - zM
print("\n[dmoi] binding-vs-methylation DISAGREEMENT (z bind - z meth):"); print(DIS.round(2).to_string())
L.to_csv(os.path.join(BASE, "gse57577_localization.csv")); M.to_csv(os.path.join(BASE, "gse57577_ctx_methylation.csv"))

# ---- figures (assayed contexts only) ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"width_ratios":[1,1,0.5]})
ctxs = list(CTX)
Lb = L.T[ctxs]
im = ax[0].imshow(Lb.values, cmap="Greens", aspect="auto")
ax[0].set_xticks(range(3)); ax[0].set_xticklabels([_wrap(c) for c in ctxs], fontsize=8); ax[0].set_yticks(range(3)); ax[0].set_yticklabels(VAR)
for i in range(3):
    for j in range(3): ax[0].text(j,i,f"{Lb.values[i,j]:.2f}",ha="center",va="center",fontsize=9,color="white" if Lb.values[i,j]>0.7 else "black")
plt.colorbar(im, ax=ax[0], label="mean Dnmt3a2 binding")
ratios = (L["WWD"] / L["WT"]).reindex(ctxs)
ratio_str = " / ".join(f"{r:.2f}" for r in ratios)
ax[0].set_title("Dnmt3a2 genomic localization (reproduced)\nWWD/WT = " + ratio_str + "×")
DD = DIS.T[ctxs]
im2 = ax[1].imshow(DD.values, cmap="coolwarm", vmin=-1.5, vmax=1.5, aspect="auto")
ax[1].set_xticks(range(3)); ax[1].set_xticklabels([_wrap(c) for c in ctxs], fontsize=8); ax[1].set_yticks(range(3)); ax[1].set_yticklabels(VAR)
for i in range(3):
    for j in range(3): ax[1].text(j,i,f"{DD.values[i,j]:+.2f}",ha="center",va="center",fontsize=9)
plt.colorbar(im2, ax=ax[1], label="bind - methyl (z)")
ax[1].set_title("DMOI binding↔methylation disagreement\n(+ = binds but doesn't methylate)")
ax[2].bar(VAR, [diff_score[v] for v in VAR], color=["#4c9a2a","#c44","#4c9a2a"])
for i, v in enumerate(VAR): ax[2].text(i, diff_score[v]+0.03, f"{diff_score[v]:.2f}", ha="center", fontsize=8)
ax[2].set_ylabel("differentiation-program score\n(lineage-gene log2 FPKM)"); ax[2].set_title("ESC differentiation\n(WT ✓ / WWD ✗ / R ✓)")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "gse57577_dmoi_pattern.png"), dpi=140); plt.close()
print("\n[dmoi] saved figure -> figures/gse57577_dmoi_pattern.png  + matrices (CSV)")
