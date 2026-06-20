#!/usr/bin/env python3
"""Render Figure 2 (hypothesis-as-anchor) from the recorded CSVs — fully reproducible, no data access needed.

Panels:
  A  hypothesis-as-anchor test on TCGA-BRCA LumA/B: signal each hypothesis adds beyond the proliferation
     anchor (hypothesis_anchor_results.csv), coloured by verdict.
  B  library screen of the 50 MSigDB Hallmarks (hypothesis_screen_results.csv): delta beyond the anchor,
     proliferation-type hallmarks (E2F/G2M/MYC) highlighted as EXPLAINED (add ~0).
  C  commonality/mediation re-characterization of the ER hypothesis across cohorts
     (hypothesis_metabric_diagnosis.csv): unique-vs-common R² stacked bars with the collinearity label —
     TCGA = NOVEL (ER variance unique beyond proliferation) vs METABRIC = REDUNDANT (shared/mediated).

Run:  python reports/fig_hypothesis_anchor.py
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
FIGS = os.path.join(HERE, "figs"); os.makedirs(FIGS, exist_ok=True)
CV = {"SUPPORTED": "#2c7fb8", "EXPLAINED_BY_TEXTBOOK": "#bdbdbd", "REFUTED": "#d95f5f"}


def main():
    a = pd.read_csv(os.path.join(REPO, "hypothesis_anchor_results.csv"))
    s = pd.read_csv(os.path.join(REPO, "hypothesis_screen_results.csv"))
    d = pd.read_csv(os.path.join(REPO, "hypothesis_metabric_diagnosis.csv"))
    tp_path = os.path.join(REPO, "transportability_sweep.csv")
    tp = pd.read_csv(tp_path) if os.path.exists(tp_path) else None

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2))
    axA, axB, axC, axD = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # --- Panel A: the three demonstration hypotheses --------------------------------------------------
    a = a.sort_values("delta_hyp_given_textbook")
    labels = a["hypothesis"].str.replace("_", " ").values
    axA.barh(labels, a["delta_hyp_given_textbook"], color=[CV.get(v, "#888") for v in a["verdict"]])
    for yi, (dv, vd) in enumerate(zip(a["delta_hyp_given_textbook"], a["verdict"])):
        axA.text(dv + 0.0008, yi, f"{dv:+.3f}", va="center", fontsize=9)
    axA.axvline(0.01, ls="--", lw=1, color="#555"); axA.text(0.0105, -0.45, "add-threshold", fontsize=7.5, color="#555")
    axA.set_xlabel("Δ signal added beyond the textbook anchor")
    axA.set_title("A  Hypothesis-as-anchor (TCGA LumA/B)", loc="left", fontweight="bold", fontsize=11)
    axA.set_xlim(-0.004, max(0.05, a["delta_hyp_given_textbook"].max() * 1.3))

    # --- Panel B: library screen ----------------------------------------------------------------------
    s = s.sort_values("delta_beyond_textbook", ascending=False)
    top = s.head(8).iloc[::-1]
    prolif_mask = top["hypothesis"].str.contains("E2F_TARGETS|G2M_CHECKPOINT|MYC_TARGETS", regex=True)
    names = top["hypothesis"].str.replace("HALLMARK_", "").str.replace("_", " ")
    axB.barh(names, top["delta_beyond_textbook"], color=[CV.get(v, "#888") for v in top["verdict"]])
    # show the proliferation hallmarks (add ~0) as a reference cluster
    pl = s[s["hypothesis"].str.contains("E2F_TARGETS|G2M_CHECKPOINT|MYC_TARGETS")]
    axB.axvline(0.0, lw=1, color="#333")
    axB.set_xlabel("Δ beyond the proliferation anchor")
    axB.set_title("B  Hallmark library screen (50 sets)", loc="left", fontweight="bold", fontsize=11)
    axB.annotate(f"proliferation hallmarks\n(E2F/G2M/MYC): Δ≈0 → EXPLAINED",
                 xy=(0.0, 0.4), xytext=(0.006, 0.8), fontsize=8, color="#444",
                 arrowprops=dict(arrowstyle="->", color="#888", lw=1))

    # --- Panel C: ER commonality TCGA vs METABRIC -----------------------------------------------------
    d = d.set_index("cohort").reindex(["TCGA", "METABRIC"])
    uniq = d["unique_er_r2"].clip(lower=0).values
    comm = d["common_r2"].clip(lower=0).values
    x = np.arange(len(d))
    axC.bar(x, uniq, color="#2c7fb8", label="unique to ER (beyond proliferation)")
    axC.bar(x, comm, bottom=uniq, color="#cfd8e3", label="shared with proliferation (common)")
    for xi, (lab, r, pm, u) in enumerate(zip(d["collinearity_label"], d["redundancy"],
                                             d["prop_mediated"], d["unique_er_r2"])):
        ytop = uniq[xi] + comm[xi]
        axC.text(xi, ytop + 0.004, lab, ha="center", fontweight="bold", fontsize=10,
                 color="#2c7fb8" if lab == "NOVEL" else "#d95f5f")
        note = f"unique R²={u:.3f}\northogonal to prolif." if lab == "NOVEL" else \
               f"redundancy {r:.2f}\n{pm*100:.0f}% mediated by prolif."
        axC.text(xi, -0.013, note, ha="center", fontsize=8, color="#444")
    axC.set_xticks(x); axC.set_xticklabels(d.index)
    axC.set_ylabel("R² explained in y (LumA vs LumB)")
    axC.set_ylim(-0.03, max(0.14, (uniq + comm).max() * 1.35))
    axC.set_title("C  ER hypothesis: NOVEL vs REDUNDANT across cohorts", loc="left", fontweight="bold", fontsize=11)
    axC.legend(loc="upper right", fontsize=8, frameon=False)

    # --- Panel D: transportability sweep (verdict vs nuisance correlation) -----------------------------
    axD.set_title("D  Transportability: verdict is set by corr(anchor, hypothesis)", loc="left",
                  fontweight="bold", fontsize=11)
    if tp is not None:
        axD.plot(tp["obs_corr"], tp["frac_novel"], "-o", color="#2c7fb8", lw=2, ms=3.5,
                 label="fraction NOVEL")
        axD.plot(tp["obs_corr"], tp["frac_redundant"], "-s", color="#d95f5f", lw=1.6, ms=3,
                 label="fraction REDUNDANT")
        obs = {"TCGA": float(d.loc["TCGA", "corr_prolif_er"]) if "corr_prolif_er" in d.columns else 0.19,
               "METABRIC": float(d.loc["METABRIC", "corr_prolif_er"]) if "corr_prolif_er" in d.columns else -0.17}
        for k, c in obs.items():
            i = (tp["obs_corr"] - c).abs().idxmin(); r = tp.loc[i]
            col = "#1a9850" if k == "TCGA" else "#d73027"
            axD.axvline(r["obs_corr"], ls="--", lw=1.2, color=col)
            axD.annotate(f"{k}\ncorr={r['obs_corr']:+.2f}\n{int(round(r['frac_novel']*100))}% NOVEL",
                         xy=(r["obs_corr"], r["frac_novel"]),
                         xytext=(r["obs_corr"], 0.62 if k == "TCGA" else 0.30), ha="center", fontsize=8,
                         color=col, fontweight="bold", arrowprops=dict(arrowstyle="->", color=col, lw=1))
        axD.set_xlabel("corr(anchor, hypothesis) — a covariate-distribution property")
        axD.set_ylabel("fraction of simulated cohorts"); axD.set_ylim(-0.03, 1.03)
        axD.legend(loc="center left", fontsize=8, frameon=False)
    else:
        axD.text(0.5, 0.5, "transportability_sweep.csv not found", ha="center", va="center", fontsize=9)

    legend = [Patch(fc=CV["SUPPORTED"], label="SUPPORTED"),
              Patch(fc=CV["EXPLAINED_BY_TEXTBOOK"], label="EXPLAINED_BY_TEXTBOOK"),
              Patch(fc=CV["REFUTED"], label="REFUTED")]
    axA.legend(handles=legend, loc="lower right", fontsize=7.5, frameon=False)

    for ax in (axA, axB, axC, axD):
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = os.path.join(FIGS, "hypothesis_anchor.png")
    fig.savefig(out, dpi=150, bbox_inches="tight"); print("wrote", out)


if __name__ == "__main__":
    main()
