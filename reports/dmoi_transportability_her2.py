#!/usr/bin/env python3
"""Transportability of the hypothesis verdict on a SECOND endpoint — HER2 status (specificity replicate).

Repeats the transportability analysis with a different textbook anchor and endpoint: anchor = the ERBB2 17q12
amplicon (9 genes), endpoint = HER2+ vs HER2- (TCGA IHC; METABRIC SNP6 GAIN vs NEUTRAL), hypothesis = the ER
signature (ER and HER2 are inversely related, so ER carries a partial HER2 signal whose collinearity with the
amplicon can differ across cohorts). For each cohort it measures the marginal effects (Cohen's d of anchor and
hypothesis vs HER2) and corr(anchor, hypothesis); it then runs the same controlled sweep (fix marginal effects,
vary only the nuisance correlation) and locates the two cohorts on the curve.

Writes transportability_her2.csv (sweep) and transportability_her2_diag.csv (per-cohort marginals) plus
reports/figs/transportability_her2.png.

Run:  BRCA_DIR=/path/to/tcga_brca METABRIC_DIR=/path/to/metabric python reports/dmoi_transportability_her2.py
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
B = os.environ.get("BRCA_DIR", ""); M = os.environ.get("METABRIC_DIR", "")
AMP = ["ERBB2", "GRB7", "STARD3", "PGAP3", "TCAP", "PNMT", "PSMD3", "GSDMB", "ORMDL3"]
ER = ["ESR1", "GATA3", "FOXA1", "XBP1", "TFF1", "PGR", "GREB1", "CA12", "SLC39A6", "NAT1", "AR", "MLPH"]


def _cohen_d(x, g):
    a, b = x[g == 0], x[g == 1]
    sp = np.sqrt((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2) + 1e-9
    return (b.mean() - a.mean()) / sp


def _diag(name, E, y):
    A = mo.signature_score(E.loc[[g for g in AMP if g in E.index]].T, [g for g in AMP if g in E.index])
    H = mo.signature_score(E.loc[[g for g in ER if g in E.index]].T, [g for g in ER if g in E.index])
    cd = mo.commonality_decomposition(A, H, y)
    return dict(cohort=name, n=int(len(y)), her2pos=int(y.sum()),
                cohen_d_amplicon=round(float(_cohen_d(A, y)), 3), cohen_d_er=round(float(_cohen_d(H, y)), 3),
                corr_amplicon_er=round(float(np.corrcoef(A, H)[0, 1]), 3),
                unique_er_r2=cd["unique_hypothesis_r2"], redundancy=cd["redundancy"],
                prop_mediated=cd["prop_mediated"], collinearity_label=cd["collinearity_label"])


def main():
    assert B and os.path.isdir(B) and M and os.path.isdir(M), "set BRCA_DIR and METABRIC_DIR"
    ex = pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0); ex = ex[~ex.index.duplicated()]
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    her2 = cl["HER2_Final_Status_nature2012"].reindex(ex.columns); m = her2.isin(["Positive", "Negative"])
    r1 = _diag("TCGA", ex.loc[:, m], (her2[m] == "Positive").astype(int).values)

    mr = pd.read_csv(os.path.join(M, "mrna_microarray.txt"), sep="\t").drop(
        columns=["Entrez_Gene_Id"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
    cp = pd.read_csv(os.path.join(M, "clinical_patient.txt"), sep="\t", comment="#")
    snp = cp.set_index(cp.columns[0])["HER2_SNP6"]
    s = [c for c in mr.columns if c in snp.index and snp[c] in ("GAIN", "NEUTRAL")]
    y = np.array([1 if snp[c] == "GAIN" else 0 for c in s])
    Em = mr[s].apply(pd.to_numeric, errors="coerce"); Em = Em.T.fillna(Em.mean(1)).T
    r2 = _diag("METABRIC", Em, y)

    diag = pd.DataFrame([r1, r2]); diag.to_csv(os.path.join(REPO, "transportability_her2_diag.csv"), index=False)
    print(diag.to_string(index=False))

    # HER2 is a SPECIFICITY contrast, not a clean same-marginal flip: the ER->HER2 effect itself differs across
    # cohorts (TCGA ER ~null -> INERT; METABRIC ER moderate but amplicon-collinear -> REDUNDANT). Left panel:
    # the real per-cohort labels (unique vs shared R²). Right panel: the sweep at METABRIC's marginals, where it
    # is faithful, reproducing the REDUNDANT valley at the observed correlation.
    mb = diag.set_index("cohort").loc["METABRIC"]
    rho_grid = np.round(np.linspace(-0.4, 0.6, 21), 3)
    rows = mo.transportability_sweep(rho_grid, d_anchor=abs(mb["cohen_d_amplicon"]), d_hyp=mb["cohen_d_er"],
                                     n=600, reps=40, seed=0)
    sweep = pd.DataFrame(rows); sweep["cohort"] = "METABRIC"
    sweep.to_csv(os.path.join(REPO, "transportability_her2.csv"), index=False)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 5.0))
    # left: real unique-vs-common R² per cohort
    d2 = diag.set_index("cohort").reindex(["TCGA", "METABRIC"])
    uniq = d2["unique_er_r2"].clip(lower=0).values
    # common R² = r2H - unique; recover r2H from redundancy (common = redundancy*r2H -> r2H = unique/(1-redundancy))
    comm = []
    for u, red in zip(uniq, d2["redundancy"].values):
        r2h = u / max(1e-6, (1 - red)) if red < 1 else u + 0.02
        comm.append(max(r2h - u, 0.0))
    comm = np.array(comm); xx = np.arange(2)
    axL.bar(xx, uniq, color="#2c7fb8", label="unique to ER (beyond amplicon)")
    axL.bar(xx, comm, bottom=uniq, color="#cfd8e3", label="shared with amplicon (common)")
    for xi, lab in enumerate(d2["collinearity_label"]):
        axL.text(xi, uniq[xi] + comm[xi] + 0.0015, lab, ha="center", fontweight="bold", fontsize=11,
                 color="#1a9850" if lab == "NOVEL" else ("#d95f5f" if lab == "REDUNDANT" else "#888"))
    axL.set_xticks(xx); axL.set_xticklabels(d2.index)
    axL.set_ylabel("R² explained in HER2 status"); axL.legend(loc="upper left", fontsize=8.5, frameon=False)
    axL.set_title("A  ER→HER2: INERT (TCGA) vs REDUNDANT (METABRIC)", loc="left", fontweight="bold", fontsize=11)
    axL.spines[["top", "right"]].set_visible(False)
    # right: METABRIC sweep
    axR.plot(sweep["obs_corr"], sweep["frac_novel"], "-o", color="#2c7fb8", lw=2, ms=3.5, label="fraction NOVEL")
    axR.plot(sweep["obs_corr"], sweep["frac_redundant"], "-s", color="#d95f5f", lw=1.6, ms=3, label="fraction REDUNDANT")
    c = float(mb["corr_amplicon_er"]); i = (sweep["obs_corr"] - c).abs().idxmin(); r = sweep.loc[i]
    axR.axvline(r["obs_corr"], ls="--", lw=1.2, color="#d73027")
    axR.annotate(f"METABRIC\ncorr={r['obs_corr']:+.2f}\nREDUNDANT", xy=(r["obs_corr"], r["frac_redundant"]),
                 xytext=(r["obs_corr"] - 0.04, 0.55), ha="center", fontsize=9, color="#d73027", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="#d73027", lw=1))
    axR.set_xlabel("corr(amplicon, ER) — covariate-distribution property")
    axR.set_ylabel("fraction of simulated cohorts"); axR.set_ylim(-0.03, 1.03)
    axR.legend(loc="center left", fontsize=8.5, frameon=False)
    axR.set_title(f"B  METABRIC sweep (d_ER={mb['cohen_d_er']:+.2f}): collinear → REDUNDANT", loc="left",
                  fontweight="bold", fontsize=11)
    axR.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = os.path.join(HERE, "figs", "transportability_her2.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")

    for _, rr in diag.iterrows():
        print(f"  {rr['cohort']:9} d_ER={rr['cohen_d_er']:+.3f} corr={rr['corr_amplicon_er']:+.3f} "
              f"-> real label {rr['collinearity_label']} (redundancy {rr['redundancy']}, {rr['prop_mediated']*100:.0f}% mediated)")
    print("wrote transportability_her2.csv, transportability_her2_diag.csv, figs/transportability_her2.png")
    print("\nHER2 is a specificity contrast: in TCGA ER is too weak for HER2 (INERT, marginal-driven); in "
          "METABRIC ER is moderate but amplicon-collinear (REDUNDANT). The framework distinguishes 'absent/weak' "
          "(INERT) from 'collinear' (REDUNDANT) -- a different regime from the ER/LumB same-marginal flip.")


if __name__ == "__main__":
    main()
