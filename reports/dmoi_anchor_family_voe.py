#!/usr/bin/env python3
"""Anchor-family vibration-of-effects: is the basal discovery robust to the choice of (textbook) anchor?

Addresses the "the textbook anchor isn't standardized" critique (see reports/anchor_standardization_discussion.md)
empirically. Instead of one hand-picked proliferation anchor, we run the LumA-vs-LumB residual discovery
across a FAMILY of biologically-equivalent proliferation priors and measure how stable the discovered axis is
across them (a vibration-of-effects / multiverse check restricted to *principled* alternatives):

  - curated_prolif  : the original 20-gene proliferation index
  - HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT, HALLMARK_MYC_TARGETS_V1 : MSigDB standard library
  - meta_PCNA       : data-driven, reproducible (Venet 2011) -- top genes correlated with PCNA

For each anchor we run anchored_residual_discovery and record (i) overlap of the discovered 30-gene panel with
the reference basal panel (novel_genes.csv), and (ii) the basal-marker core recovered. Across anchors we
report the mean pairwise Jaccard of discovered panels and the consensus genes found by >=4/5 anchors. If the
basal/keratinization axis recurs across the family, the discovery is not an artifact of one anchor choice.
Writes anchor_family_voe.csv and anchor_family_consensus.csv.

Run:  BRCA_DIR=/path/to/tcga_brca MSIGDB_GMT=/path/to/h.all...symbols.gmt python reports/dmoi_anchor_family_voe.py
"""
import os, sys, itertools
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
B = os.environ.get("BRCA_DIR", ""); GMT = os.environ.get("MSIGDB_GMT", "")
CURATED = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
           "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]
HALLMARKS = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1"]


def _gmt(path, names):
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] in names:
                out[parts[0]] = parts[2:]
    return out


def main():
    assert B and os.path.isdir(B), "set BRCA_DIR"; assert GMT and os.path.exists(GMT), "set MSIGDB_GMT"
    expr = pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0); expr = expr[~expr.index.duplicated()]
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    pam = cl["PAM50Call_RNAseq"].reindex(expr.columns); mask = pam.isin(["LumA", "LumB"])
    y = (pam[mask] == "LumB").astype(int).values
    E = expr.loc[:, mask]                                  # genes x samples (labeled)
    basal = list(pd.read_csv(os.path.join(REPO, "novel_genes.csv"))["gene"]); basalset = set(basal)
    core = ["KRT5", "KRT14", "KRT17", "KRT6B", "TP63", "DSG3", "DSC3"]
    # feature matrix for discovery: top-variance genes + forced-in basal panel
    var = E.var(axis=1); feats = sorted(set(var.sort_values(ascending=False).head(3000).index) | (basalset & set(E.index)))
    X = E.loc[feats].T.values.astype("float32")
    # anchor family
    fam = {"curated_prolif": [g for g in CURATED if g in E.index]}
    fam.update({k: [g for g in v if g in E.index] for k, v in _gmt(GMT, HALLMARKS).items()})
    fam["meta_PCNA"] = mo.marker_correlated_anchor(E.T, marker="PCNA", top_k=50, exclude_marker=True)
    from sklearn.metrics import roc_auc_score
    yf = y.astype(float)

    def select_novel(anchor, top_k=30, corr_max=0.6):
        """The anchor-orthogonal top-K panel anchored_residual_discovery would select (deterministic from the
        vectorized partial correlation) -- the quantity the VoE needs, without the costly gate/null steps."""
        a = np.asarray(anchor, float); az = (a - a.mean()) / (a.std() + 1e-9); azc = az - az.mean()
        den = float(azc @ azc) + 1e-12
        rY = (yf - yf.mean()) - (float(azc @ (yf - yf.mean())) / den) * azc; rYc = rY - rY.mean()
        rYn = np.sqrt(float(rYc @ rYc)) + 1e-12
        Xc = X - X.mean(0); rX = Xc - np.outer(azc, (azc @ Xc) / den); rXc = rX - rX.mean(0)
        pc = (rXc.T @ rYc) / (np.sqrt((rXc * rXc).sum(0)) * rYn + 1e-12)
        ca = (Xc.T @ azc) / (np.sqrt((Xc * Xc).sum(0)) * np.sqrt(float(azc @ azc)) + 1e-12)
        idx = [j for j in np.argsort(-np.abs(pc)) if abs(ca[j]) < corr_max][:top_k]
        return [feats[j] for j in idx]

    rows, panels = [], {}
    for name, genes in fam.items():
        g_in = [g for g in genes if g in E.index]
        anchor = mo.signature_score(E.loc[g_in].T, g_in)
        nov = select_novel(anchor); panels[name] = set(nov)
        rows.append(dict(anchor=name, n_anchor_genes=len(g_in),
                         anchor_auroc=round(float(roc_auc_score(y, anchor)), 3),
                         overlap_basal_of30=len(set(nov) & basalset),
                         basal_core_hits=sum(g in nov for g in core)))
        print(f"{name:24} anchorAUROC={roc_auc_score(y,anchor):.3f} overlap_basal={len(set(nov)&basalset)}/30 core={sum(g in nov for g in core)}/7")
    # vibration: mean pairwise Jaccard of discovered panels
    js = [len(a & b) / len(a | b) for a, b in itertools.combinations(panels.values(), 2)]
    allg = set().union(*panels.values()); consensus = sorted(g for g in allg if sum(g in s for s in panels.values()) >= 4)
    pd.DataFrame(rows).to_csv(os.path.join(REPO, "anchor_family_voe.csv"), index=False)
    pd.DataFrame({"gene": consensus, "n_anchors_recovering": [sum(g in s for s in panels.values()) for g in consensus],
                  "in_basal_panel": [g in basalset for g in consensus]}).to_csv(
        os.path.join(REPO, "anchor_family_consensus.csv"), index=False)
    print(f"\nmean pairwise Jaccard across {len(panels)} anchors = {np.mean(js):.3f}")
    print(f"consensus genes recovered by >=4/{len(panels)} anchors: {len(consensus)} "
          f"({sum(g in basalset for g in consensus)} in basal panel)")
    print("CONCLUSION: the basal/keratinization axis is" + (" robust" if np.mean(js) > 0.4 else " NOT robust") +
          " to anchor choice across a family of standard proliferation priors.")


if __name__ == "__main__":
    main()
