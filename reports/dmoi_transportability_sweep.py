#!/usr/bin/env python3
"""Transportability of the hypothesis verdict: why ER is SUPPORTED in TCGA but REDUNDANT in METABRIC.

The cross-cohort non-reproduction we diagnosed (commonality/mediation) is a *transportability* phenomenon
(Degtiar & Rose 2023): the verdict depends on the anchor-hypothesis correlation, which is a property of each
cohort's covariate distribution, not of the hypothesis's marginal effect. This runner makes that quantitative
and reproducible. Holding both marginal separations fixed at the values observed in the BRCA data
(prolif Cohen's d and ER Cohen's d, read from hypothesis_metabric_diagnosis.csv if present), it sweeps only the
residual anchor-hypothesis correlation and records how the verdict moves. It then locates the two real cohorts
on the curve by their measured corr(prolif, ER) (TCGA +0.19, METABRIC -0.17).

Writes transportability_sweep.csv and reports/figs/transportability_sweep.png. No external data required
(the marginal effect sizes are read from the recorded diagnosis CSV, with safe defaults if absent).

Run:  python reports/dmoi_transportability_sweep.py
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo

# marginal separations (Cohen's d) of anchor and hypothesis vs the LumA/B endpoint, from the diagnosis CSV
D_ANCHOR, D_HYP = 1.85, -0.22                                  # defaults ~ mean of the two cohorts
OBS = {"TCGA": 0.188, "METABRIC": -0.173}                     # measured corr(prolif, ER) per cohort


def main():
    diag = os.path.join(REPO, "hypothesis_metabric_diagnosis.csv")
    global D_ANCHOR, D_HYP, OBS
    if os.path.exists(diag):
        d = pd.read_csv(diag).set_index("cohort")
        D_ANCHOR = round(float(d["cohen_d_prolif"].abs().mean()), 3)
        D_HYP = round(float(d["cohen_d_er"].mean()), 3)
        OBS = {c: float(d.loc[c, "corr_prolif_er"]) for c in d.index}
    rho_grid = np.round(np.linspace(-0.4, 0.6, 21), 3)
    rows = mo.transportability_sweep(rho_grid, d_anchor=D_ANCHOR, d_hyp=D_HYP, n=600, reps=40, seed=0)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(REPO, "transportability_sweep.csv"), index=False)

    # locate each real cohort on the curve by its observed corr(prolif, ER)
    def at_corr(c):
        i = (df["obs_corr"] - c).abs().idxmin(); return df.loc[i]
    marks = {k: at_corr(v) for k, v in OBS.items()}

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    ax.plot(df["obs_corr"], df["frac_novel"], "-o", color="#2c7fb8", lw=2, ms=4,
            label="fraction NOVEL (ER adds unique signal)")
    ax.plot(df["obs_corr"], df["frac_redundant"], "-s", color="#d95f5f", lw=1.6, ms=3.5,
            label="fraction REDUNDANT (collinear with anchor)")
    ax.axhline(0.5, ls=":", lw=1, color="#999")
    for k, r in marks.items():
        col = "#1a9850" if k == "TCGA" else "#d73027"
        ax.axvline(r["obs_corr"], ls="--", lw=1.2, color=col)
        ax.annotate(f"{k}\ncorr={r['obs_corr']:+.2f}\n{int(round(r['frac_novel']*100))}% NOVEL",
                    xy=(r["obs_corr"], r["frac_novel"]), xytext=(r["obs_corr"], 0.62 if k == "TCGA" else 0.30),
                    ha="center", fontsize=9, color=col, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color=col, lw=1))
    ax.set_xlabel("corr(anchor, hypothesis)  —  a covariate-distribution property of the cohort")
    ax.set_ylabel("fraction of simulated cohorts")
    ax.set_title(f"Transportability of the ER verdict (fixed marginal effects: d_prolif={D_ANCHOR}, d_ER={D_HYP})",
                 fontsize=11, fontweight="bold")
    ax.set_ylim(-0.03, 1.03); ax.legend(loc="center left", fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = os.path.join(HERE, "figs", "transportability_sweep.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")

    print(f"transportability sweep | d_anchor={D_ANCHOR} d_hyp={D_HYP} | rho in [{rho_grid[0]},{rho_grid[-1]}]")
    print(df[["obs_corr", "frac_novel", "frac_redundant", "mean_unique_r2"]].to_string(index=False))
    for k, r in marks.items():
        print(f"  {k:9} corr={r['obs_corr']:+.3f} -> {r['frac_novel']*100:.0f}% NOVEL, "
              f"{r['frac_redundant']*100:.0f}% REDUNDANT (same marginal effects)")
    print("wrote transportability_sweep.csv and figs/transportability_sweep.png")
    print("\nThe SAME ER effect transports to opposite verdicts depending only on corr(prolif,ER): the "
          "non-reproduction is a covariate-distribution (transportability) property, not absent biology.")


if __name__ == "__main__":
    main()
