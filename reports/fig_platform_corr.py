#!/usr/bin/env python3
"""Why the LumA-vs-LumB (proliferation->ER) label tracks measurement technology: the SIGN of the
proliferation-ER correlation flips between RNA-seq and microarray platforms.

Reproducible from the recorded endpoint_panel.csv (no data access). For the LumA-vs-LumB endpoint it plots
corr(proliferation, ER) per cohort, coloured by platform family (RNA-seq vs microarray), with the resulting
commonality label. The two RNA-seq cohorts (TCGA RNA-seq, SCAN-B) have a POSITIVE correlation -> ER carries
suppression/unique variance -> NOVEL; the two microarrays (TCGA Agilent, METABRIC) have a NEGATIVE correlation
-> ER is collinear with proliferation in the LumB direction -> REDUNDANT. This is the mechanism behind the
transportability sweep (the verdict is governed by corr(anchor, hypothesis)).

Run:  python reports/fig_platform_corr.py
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
PLATFORM = {"TCGA_RNAseq": "RNA-seq", "SCAN-B": "RNA-seq", "TCGA_Agilent": "microarray", "METABRIC": "microarray"}
PCOL = {"RNA-seq": "#2c7fb8", "microarray": "#d95f5f"}


def main():
    d = pd.read_csv(os.path.join(REPO, "endpoint_panel.csv"))
    d = d[d["endpoint"] == "LumA_vs_LumB"].copy()
    d["platform"] = d["cohort"].map(PLATFORM)
    d = d.sort_values(["platform", "corr_anchor_hyp"])
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    x = np.arange(len(d))
    bars = ax.bar(x, d["corr_anchor_hyp"], color=[PCOL[p] for p in d["platform"]], width=0.62, alpha=0.9)
    ax.axhline(0, color="#333", lw=1)
    for xi, (c, lab, n) in enumerate(zip(d["corr_anchor_hyp"], d["collinearity_label"], d["n"])):
        va = "bottom" if c >= 0 else "top"
        ax.text(xi, c + (0.012 if c >= 0 else -0.012), f"{c:+.2f}", ha="center", va=va, fontsize=10, fontweight="bold")
        ax.text(xi, 0.005 if c < 0 else -0.02, lab, ha="center", va="bottom" if c < 0 else "top",
                fontsize=9, fontweight="bold", color="#2c7fb8" if lab == "NOVEL" else "#d95f5f")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n(n={n})" for c, n in zip(d["cohort"], d["n"])], fontsize=9)
    ax.set_ylabel("corr(proliferation, ER)  in Luminal A vs B")
    ax.set_title("Luminal A/B label tracks platform: corr(prolif, ER) sign flips RNA-seq vs microarray",
                 fontsize=11, fontweight="bold")
    lim = max(0.25, float(d["corr_anchor_hyp"].abs().max()) * 1.5); ax.set_ylim(-lim, lim)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=PCOL[p]) for p in ["RNA-seq", "microarray"]]
    ax.legend(handles, ["RNA-seq (→ NOVEL)", "microarray (→ REDUNDANT)"], loc="upper right", fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = os.path.join(HERE, "figs", "platform_corr.png")
    fig.savefig(out, dpi=150, bbox_inches="tight"); print("wrote", out)
    print(d[["cohort", "platform", "corr_anchor_hyp", "redundancy", "collinearity_label"]].to_string(index=False))


if __name__ == "__main__":
    main()
