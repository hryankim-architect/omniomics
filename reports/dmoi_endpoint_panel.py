#!/usr/bin/env python3
"""Multi-endpoint x cohort panel: characterise (anchor, hypothesis) pairs across four breast-cancer endpoints
in two cohorts, and flag which characterisations transport.

Generalises the single ER/LumA-B and HER2 results into a systematic grid. For each endpoint we anchor on the
textbook driver and test a secondary hypothesis with hypothesis_anchor_test, recording the gated verdict plus
the gate-free commonality label (NOVEL / REDUNDANT / INERT). Running the same four endpoints in TCGA-BRCA and
METABRIC then shows where the label is concordant (transportable) and where it changes — the practical payload
of the transportability analysis.

Endpoints (endpoint | anchor -> hypothesis):
  LumA_vs_LumB     proliferation -> ER          (incomplete anchor; ER is the orthogonal lineage axis)
  HER2_pos_vs_neg  ERBB2 amplicon -> ER         (ER as a secondary HER2 signal)
  ER_pos_vs_neg    ER signature -> proliferation(near-complete anchor; expect nothing to add)
  Basal_vs_rest    basal/keratinization -> immune(does the immune program add beyond the basal axis?)

A third column (TCGA_Agilent) is added automatically when an Agilent microarray matrix is present in BRCA_DIR
(AgilentG4502A_07_3.gz from UCSC Xena: https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/AgilentG4502A_07_3.gz).
It uses the SAME TCGA patients on a different platform, so its agreement with TCGA_RNAseq is a platform-
transportability check, and its agreement with METABRIC is cohort-transportability. If the file is absent the
runner simply produces the two-cohort panel.

A fourth column (SCAN-B) is added when SCANB_DIR holds scanb_markers.csv + scanb_pheno.csv — a fully
independent Swedish RNA-seq cohort (GEO GSE96058, ~3,400 tumours). These are pre-extracted from the public
GSE96058 gene-expression matrix (marker rows only) and the series-matrix phenotype (pam50/er/her2).

Writes endpoint_panel.csv and reports/figs/endpoint_panel.png.

Run:  BRCA_DIR=/path/to/tcga_brca METABRIC_DIR=/path/to/metabric python reports/dmoi_endpoint_panel.py
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
B = os.environ.get("BRCA_DIR", ""); M = os.environ.get("METABRIC_DIR", "")

PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]
ER = ["ESR1", "GATA3", "FOXA1", "XBP1", "TFF1", "PGR", "GREB1", "CA12", "SLC39A6", "NAT1", "AR", "MLPH"]
AMP = ["ERBB2", "GRB7", "STARD3", "PGAP3", "TCAP", "PNMT", "PSMD3", "GSDMB", "ORMDL3"]
BASAL = ["KRT5", "KRT14", "KRT17", "KRT6B", "TP63", "DSG3", "DSC3", "SOX10", "COL17A1", "FOXC1", "MIA", "SFRP1"]
IMMUNE = ["CD8A", "CD3D", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10", "GZMK", "NKG7", "CCL5", "CD2", "PTPRC"]

# endpoint name -> (anchor genes, anchor label, hypothesis genes, hypothesis label)
ENDPOINTS = {
    "LumA_vs_LumB":    (PROLIF, "proliferation", ER, "ER"),
    "HER2_pos_vs_neg": (AMP, "ERBB2_amplicon", ER, "ER"),
    "ER_pos_vs_neg":   (ER, "ER_signature", PROLIF, "proliferation"),
    "Basal_vs_rest":   (BASAL, "basal", IMMUNE, "immune"),
}


def _score(E, genes):
    g = [x for x in genes if x in E.index]
    return mo.signature_score(E.loc[g].T, g)


def _cell(cohort, endpoint, E, y):
    ag, alab, hg, hlab = ENDPOINTS[endpoint]
    A = _score(E, ag); H = _score(E, hg)
    r = mo.hypothesis_anchor_test(A, H, y, cv=4, random_state=0, inner_repeats=1)
    return dict(endpoint=endpoint, cohort=cohort, n=int(len(y)), pos=int(y.sum()),
                anchor=alab, hypothesis=hlab,
                auroc_anchor=round(max(roc_auc_score(y, A), 1 - roc_auc_score(y, A)), 3),
                auroc_hyp=round(max(roc_auc_score(y, H), 1 - roc_auc_score(y, H)), 3),
                corr_anchor_hyp=r["corr_textbook_hypothesis"], delta_beyond=r["delta_hyp_given_textbook"],
                verdict=r["verdict"], unique_r2=r["unique_hypothesis_r2"], redundancy=r["redundancy"],
                prop_mediated=r["prop_mediated"], collinearity_label=r["collinearity_label"])


def _tcga_cohort(name, ex, cl):
    """Run the four endpoints for one TCGA expression matrix (RNA-seq or Agilent) against the shared clinical."""
    ex = ex[~ex.index.duplicated()]
    pam = cl["PAM50Call_RNAseq"].reindex(ex.columns)
    er = cl["breast_carcinoma_estrogen_receptor_status"].reindex(ex.columns)
    her2 = cl["HER2_Final_Status_nature2012"].reindex(ex.columns)
    out = []
    m = pam.isin(["LumA", "LumB"]); out.append(_cell(name, "LumA_vs_LumB", ex.loc[:, m], (pam[m] == "LumB").astype(int).values))
    m = her2.isin(["Positive", "Negative"]); out.append(_cell(name, "HER2_pos_vs_neg", ex.loc[:, m], (her2[m] == "Positive").astype(int).values))
    m = er.isin(["Positive", "Negative"]); out.append(_cell(name, "ER_pos_vs_neg", ex.loc[:, m], (er[m] == "Positive").astype(int).values))
    m = pam.isin(["LumA", "LumB", "Basal", "Her2", "Normal"]); out.append(_cell(name, "Basal_vs_rest", ex.loc[:, m], (pam[m] == "Basal").astype(int).values))
    return out


def _tcga():
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    out = _tcga_cohort("TCGA_RNAseq", pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0), cl)
    ag = os.environ.get("AGILENT_FILE", os.path.join(B, "AgilentG4502A_07_3.gz"))
    if os.path.exists(ag):   # third column: same patients, Agilent microarray platform (platform-transportability)
        am = pd.read_csv(ag, sep="\t", index_col=0)
        am = am.apply(pd.to_numeric, errors="coerce"); am = am.T.fillna(am.mean(1)).T
        out += _tcga_cohort("TCGA_Agilent", am, cl)
    return out


def _metabric():
    mr = pd.read_csv(os.path.join(M, "mrna_microarray.txt"), sep="\t").drop(
        columns=["Entrez_Gene_Id"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
    cp = pd.read_csv(os.path.join(M, "clinical_patient.txt"), sep="\t", comment="#"); cp = cp.set_index(cp.columns[0])
    sub = cp["CLAUDIN_SUBTYPE"]; er = cp["ER_IHC"].astype(str); snp = cp["HER2_SNP6"]

    def E_for(ids):
        s = [c for c in mr.columns if c in ids]
        Em = mr[s].apply(pd.to_numeric, errors="coerce"); return Em.T.fillna(Em.mean(1)).T, s
    out = []
    ids = sub[sub.isin(["LumA", "LumB"])].index; E, s = E_for(ids)
    out.append(_cell("METABRIC", "LumA_vs_LumB", E, np.array([1 if sub[c] == "LumB" else 0 for c in s])))
    ids = snp[snp.isin(["GAIN", "NEUTRAL"])].index; E, s = E_for(ids)
    out.append(_cell("METABRIC", "HER2_pos_vs_neg", E, np.array([1 if snp[c] == "GAIN" else 0 for c in s])))
    ids = er[er.str.startswith(("Pos", "Neg"))].index; E, s = E_for(ids)
    out.append(_cell("METABRIC", "ER_pos_vs_neg", E, np.array([1 if er[c].startswith("Pos") else 0 for c in s])))
    ids = sub[sub.isin(["LumA", "LumB", "Basal", "Her2", "Normal"])].index; E, s = E_for(ids)
    out.append(_cell("METABRIC", "Basal_vs_rest", E, np.array([1 if sub[c] == "Basal" else 0 for c in s])))
    return out


def _scanb():
    """Fourth column: SCAN-B (GSE96058), a fully independent Swedish RNA-seq cohort (~3,400 tumours).
    Reads the pre-extracted marker-gene matrix (scanb_markers.csv) and phenotype (scanb_pheno.csv) from
    SCANB_DIR; returns [] if absent. ER/HER2 are IHC status (1=pos, 0=neg); pam50 is the gene-expression call."""
    sd = os.environ.get("SCANB_DIR", os.path.join(os.path.dirname(B.rstrip("/")), "scanb"))
    mk = os.path.join(sd, "scanb_markers.csv"); ph = os.path.join(sd, "scanb_pheno.csv")
    if not (os.path.exists(mk) and os.path.exists(ph)):
        return []
    E = pd.read_csv(mk, index_col=0); E = E[~E.index.duplicated()]
    pk = pd.read_csv(ph).set_index("sample")
    pam = pk["pam50"].astype(str).reindex(E.columns)
    er = pd.to_numeric(pk["er"], errors="coerce").reindex(E.columns)      # 1=positive, 0=negative, NaN=unknown
    her2 = pd.to_numeric(pk["her2"], errors="coerce").reindex(E.columns)
    out = []
    m = pam.isin(["LumA", "LumB"]); out.append(_cell("SCAN-B", "LumA_vs_LumB", E.loc[:, m], (pam[m] == "LumB").astype(int).values))
    m = her2.isin([0, 1]); out.append(_cell("SCAN-B", "HER2_pos_vs_neg", E.loc[:, m], (her2[m] == 1).astype(int).values))
    m = er.isin([0, 1]); out.append(_cell("SCAN-B", "ER_pos_vs_neg", E.loc[:, m], (er[m] == 1).astype(int).values))
    m = pam.isin(["LumA", "LumB", "Basal", "Her2", "Normal"]); out.append(_cell("SCAN-B", "Basal_vs_rest", E.loc[:, m], (pam[m] == "Basal").astype(int).values))
    return out


def main():
    assert B and os.path.isdir(B) and M and os.path.isdir(M), "set BRCA_DIR and METABRIC_DIR"
    rows = _tcga() + _metabric() + _scanb()
    df = pd.DataFrame(rows)
    cohorts = [c for c in ["TCGA_RNAseq", "TCGA_Agilent", "METABRIC", "SCAN-B"] if c in set(df["cohort"])]
    # transportability concordance per endpoint (same collinearity_label across ALL cohorts present)
    lab = df.pivot(index="endpoint", columns="cohort", values="collinearity_label")
    df["transports"] = df["endpoint"].map(lambda e: bool(lab.loc[e, cohorts].nunique() == 1))
    df = df.sort_values(["endpoint", "cohort"]).reset_index(drop=True)
    df.to_csv(os.path.join(REPO, "endpoint_panel.csv"), index=False)
    print(df[["endpoint", "cohort", "anchor", "hypothesis", "verdict", "collinearity_label",
              "redundancy", "transports"]].to_string(index=False))

    # label grid figure
    order = list(ENDPOINTS)
    COL = {"NOVEL": "#2c7fb8", "REDUNDANT": "#d95f5f", "INERT": "#bdbdbd"}
    fig, ax = plt.subplots(figsize=(5.2 + 1.9 * len(cohorts), 4.9))
    for i, ep in enumerate(order):
        for j, ch in enumerate(cohorts):
            r = df[(df.endpoint == ep) & (df.cohort == ch)].iloc[0]
            ax.add_patch(plt.Rectangle((j, len(order) - 1 - i), 1, 1, facecolor=COL[r["collinearity_label"]],
                                       edgecolor="white", lw=2, alpha=0.85))
            ax.text(j + 0.5, len(order) - 1 - i + 0.62, r["collinearity_label"], ha="center", va="center",
                    fontsize=10, fontweight="bold", color="white")
            ax.text(j + 0.5, len(order) - 1 - i + 0.30, f"{r['verdict'][:11]}", ha="center", va="center",
                    fontsize=7.5, color="white")
            ax.text(j + 0.5, len(order) - 1 - i + 0.08, f"red={r['redundancy']:.2f}", ha="center", va="center",
                    fontsize=7, color="white")
    for i, ep in enumerate(order):
        ag, alab, hg, hlab = ENDPOINTS[ep]
        tp = df[df.endpoint == ep].iloc[0]["transports"]
        ax.text(-0.06, len(order) - 1 - i + 0.5, f"{ep}\n{alab}→{hlab}  {'✓transports' if tp else '✗differs'}",
                ha="right", va="center", fontsize=8.5)
    for j, ch in enumerate(cohorts):
        ax.text(j + 0.5, len(order) + 0.08, ch.replace("_", "\n"), ha="center", va="bottom", fontsize=9.5, fontweight="bold")
    ax.set_xlim(-2.4, len(cohorts)); ax.set_ylim(0, len(order) + 0.45); ax.axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in COL.values()]
    ax.legend(handles, list(COL), loc="lower left", bbox_to_anchor=(-0.55, -0.1), ncol=3, fontsize=8.5, frameon=False)
    ax.set_title(f"Anchored hypothesis labels across {len(order)} endpoints × {len(cohorts)} cohorts",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figs", "endpoint_panel.png"), dpi=150, bbox_inches="tight")
    nt = int(df.drop_duplicates("endpoint")["transports"].sum())
    print(f"\n{nt}/{len(order)} endpoints transport (same label across all {len(cohorts)} cohorts). "
          "wrote endpoint_panel.csv and figs/endpoint_panel.png")


if __name__ == "__main__":
    main()
