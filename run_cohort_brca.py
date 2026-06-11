#!/usr/bin/env python3
"""Next step — run the omniomics engine on the user's OWN multi-cohort data
(dmoi-brca-poc): TCGA-BRCA RNA-seq + METABRIC microarray, HER2-vs-Luminal.

This is a harder, more realistic harmonization than the mouse cohort: two cohorts on
*different platforms* (RNA-seq vs microarray) — a massive batch effect — with a shared
biological axis (breast-cancer subtype). We show omniomics removes the platform batch and
that the subtype signal transfers across cohorts (independently echoing the dmoi-brca-poc
headline: cross-cohort subtype AUROC ~0.92)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import loaders, harmonize as hz
from omniomics import config
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BR = config.brca_data_dir()
ROOT = os.path.dirname(__file__); FIG = os.path.join(ROOT, "figures"); os.makedirs(FIG, exist_ok=True)

# ---- expression matrices (already log2) ----
tcga = loaders.load_matrix(os.path.join(BR,"tcga_brca","HiSeqV2.gz"), gene_col=0)
meta = loaders.load_matrix(os.path.join(BR,"metabric","mrna_microarray.txt"), gene_col=0,
                           drop_cols=["Entrez_Gene_Id"])
tcga = loaders.gene_symbol_index(tcga); meta = loaders.gene_symbol_index(meta)

# ---- labels from the user's cohort files (HER2 vs Luminal) ----
def labels(path):
    c = pd.read_csv(path, sep="\t")
    c = c[c.get("has_rna", True) != False]
    return c.set_index("sample_id")["group"]
lab_t = labels(os.path.join(BR,"tcga_brca","cohort_v4.tsv"))
lab_m = labels(os.path.join(BR,"metabric","cohort_v4.tsv"))
tcga = tcga[[s for s in tcga.columns if s in set(lab_t.index)]]
meta = meta[[s for s in meta.columns if s in set(lab_m.index)]]
print(f"[brca] TCGA {tcga.shape} | METABRIC {meta.shape}")

# ---- harmonize on common gene-symbol panel (already log2 -> log=False) ----
studies = {"TCGA": tcga, "METABRIC": meta}
genes = hz.common_gene_panel(studies)
X, study_of = hz.assemble(studies, genes=genes, log=False)
X = X.dropna()                      # drop genes with missing values (microarray gaps)
Xqn = hz.quantile_normalize(X)
bv_before, _ = hz.batch_variance_explained(Xqn, study_of)
Xc = hz.combat_lite(Xqn, study_of)
bv_after, _ = hz.batch_variance_explained(Xc, study_of)
print(f"[brca] common genes {len(genes)}; combined {X.shape[1]} samples")
print(f"[brca] PC1 cohort/platform variance  before={bv_before['PC1']:.2f}  after={bv_after['PC1']:.2f}")

# group vector aligned to columns
grp = {}
for c in X.columns:
    sid = c.split("|",1)[1]
    grp[c] = (lab_t.get(sid) if c.startswith("TCGA") else lab_m.get(sid))
grp = pd.Series(grp)
y = (grp == "HER2").astype(int)

# ---- biology preservation: cross-cohort subtype transfer (train TCGA -> test METABRIC) ----
def cross_cohort_auroc(M):
    tr = study_of[study_of=="TCGA"].index; te = study_of[study_of=="METABRIC"].index
    sc = StandardScaler().fit(M[tr].T.values)
    clf = LogisticRegression(max_iter=2000, C=0.1).fit(sc.transform(M[tr].T.values), y[tr].values)
    p = clf.predict_proba(sc.transform(M[te].T.values))[:,1]
    return roc_auc_score(y[te].values, p)
auc_before = cross_cohort_auroc(Xqn)
auc_after  = cross_cohort_auroc(Xc)
print(f"[brca] cross-cohort HER2-vs-Luminal AUROC  before={auc_before:.3f}  after={auc_after:.3f}")

# ---- figures ----
def pca_plot(ax, M, color_by, title, palette):
    pcs,_ = hz.pca(M, n=2)
    for k,c in palette.items():
        m = (color_by.values==k)
        ax.scatter(pcs[m,0], pcs[m,1], s=8, c=c, label=str(k), alpha=.6, edgecolor="none")
    ax.set_title(title, fontsize=10); ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(fontsize=7)
fig,ax = plt.subplots(2,2, figsize=(11,9))
pca_plot(ax[0,0], Xqn, study_of, "Before: colored by COHORT/platform", {"TCGA":"#1f77b4","METABRIC":"#d62728"})
pca_plot(ax[0,1], Xqn, grp,      "Before: colored by SUBTYPE", {"Luminal":"#2ca02c","HER2":"#9467bd"})
pca_plot(ax[1,0], Xc,  study_of, "After ComBat: colored by COHORT/platform", {"TCGA":"#1f77b4","METABRIC":"#d62728"})
pca_plot(ax[1,1], Xc,  grp,      "After ComBat: colored by SUBTYPE", {"Luminal":"#2ca02c","HER2":"#9467bd"})
plt.tight_layout(); plt.savefig(os.path.join(FIG,"brca_harmonization_pca.png"),dpi=140); plt.close()

fig,a = plt.subplots(figsize=(5,4))
a.bar(["before\nharmonization","after\nharmonization"], [auc_before, auc_after],
      color=["#bbbbbb","#2E5A87"])
for i,v in enumerate([auc_before,auc_after]): a.text(i,v+0.01,f"{v:.3f}",ha="center")
a.axhline(0.5,ls="--",c="k",lw=.6); a.set_ylim(0.4,1.0); a.set_ylabel("cross-cohort AUROC (TCGA→METABRIC)")
a.set_title("Subtype signal transfers across cohorts\nonly after platform batch removed")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"brca_cross_cohort_auroc.png"),dpi=140); plt.close()

pd.DataFrame({"metric":["PC1_cohort_var_before","PC1_cohort_var_after",
             "xcohort_AUROC_before","xcohort_AUROC_after"],
             "value":[bv_before['PC1'],bv_after['PC1'],auc_before,auc_after]}).to_csv(
             os.path.join(ROOT,"brca_harmonization_metrics.csv"), index=False)
print("[brca] saved figures + metrics")
