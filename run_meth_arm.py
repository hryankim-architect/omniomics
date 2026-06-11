#!/usr/bin/env python3
"""Methylation arm — true multi-omics integration on matched TCGA-BRCA samples.
Pairs RNA (HiSeqV2) with HM450 promoter methylation. Two readouts:
(A) promoter methylation vs expression anti-correlation (the epigenetic silencing signature);
(B) does adding methylation improve breast-cancer subtype classification over RNA alone?"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import loaders, methylation as me
from omniomics import config
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BR = config.brca_tcga_dir()
MT = config.cache_dir()
ROOT = os.path.dirname(__file__); FIG = os.path.join(ROOT,"figures"); os.makedirs(FIG, exist_ok=True)

# ---- load both omics, match samples/genes ----
rna = loaders.gene_symbol_index(loaders.load_matrix(os.path.join(BR,"HiSeqV2.gz"), gene_col=0))
meth = me.load_hm450_promoter(os.path.join(MT,"promoter_meth.tsv.gz"), os.path.join(MT,"promoter_probe_gene.tsv"))
samples = sorted(set(rna.columns) & set(meth.columns))
genes = sorted(set(rna.index) & set(meth.index))
rna = rna.loc[genes, samples]; meth = meth.loc[genes, samples]
print(f"[meth] matched: {len(samples)} samples x {len(genes)} genes (RNA + promoter methylation)")

# ---- (A) per-gene methylation<->expression correlation across samples ----
expr_ok = (rna.std(axis=1) > 0.1) & (meth.std(axis=1) > 0.02) & (meth.mean(axis=1).between(0.05,0.95))
g_ok = rna.index[expr_ok]
rho = np.array([stats.spearmanr(meth.loc[g], rna.loc[g]).correlation for g in g_ok])
rho = rho[~np.isnan(rho)]
print(f"[meth] genes tested: {len(rho)}; median Spearman(meth,expr) = {np.median(rho):.3f}; "
      f"% negative = {100*(rho<0).mean():.1f}%; % strong-neg (<-0.3) = {100*(rho<-0.3).mean():.1f}%")

# ---- (B) multi-omics subtype classification gain ----
coh = pd.read_csv(os.path.join(BR,"cohort_v4.tsv"), sep="\t").set_index("sample_id")
lab = coh["group"].reindex(samples)
keep = lab.isin(["HER2","Luminal"])
y = (lab[keep]=="HER2").astype(int).values
S = [s for s,k in zip(samples,keep) if k]
# top variable genes for a compact feature set
var = rna[S].var(axis=1).sort_values(ascending=False)
top = var.index[:2000]
Xr = StandardScaler().fit_transform(rna.loc[top,S].T.values)
Xm = StandardScaler().fit_transform(meth.loc[top,S].fillna(meth.loc[top,S].mean()).T.values)
cv = StratifiedKFold(5, shuffle=True, random_state=0)
clf = LogisticRegression(max_iter=3000, C=0.1)
auc_r  = cross_val_score(clf, Xr, y, cv=cv, scoring="roc_auc").mean()
auc_m  = cross_val_score(clf, Xm, y, cv=cv, scoring="roc_auc").mean()
auc_rm = cross_val_score(clf, np.hstack([Xr,Xm]), y, cv=cv, scoring="roc_auc").mean()
print(f"[meth] subtype AUROC  RNA={auc_r:.3f}  METH={auc_m:.3f}  RNA+METH={auc_rm:.3f}")

pd.DataFrame({"metric":["median_meth_expr_rho","pct_negative","AUROC_RNA","AUROC_METH","AUROC_RNA+METH"],
              "value":[np.median(rho),100*(rho<0).mean(),auc_r,auc_m,auc_rm]}).to_csv(
              os.path.join(ROOT,"meth_arm_metrics.csv"), index=False)

# ---- figures ----
fig,ax=plt.subplots(1,2, figsize=(11,4.4))
ax[0].hist(rho, bins=60, color="#4477aa"); ax[0].axvline(0,c="k",lw=.8)
ax[0].axvline(np.median(rho),c="crimson",ls="--",label=f"median={np.median(rho):.2f}")
ax[0].set_xlabel("Spearman(promoter methylation, expression)"); ax[0].set_ylabel("genes")
ax[0].set_title(f"Promoter methylation silences expression\n{100*(rho<0).mean():.0f}% of genes negative"); ax[0].legend(fontsize=8)
ax[1].bar(["RNA","METH","RNA+METH"], [auc_r,auc_m,auc_rm], color=["#999","#bb8844","#2E5A87"])
for i,v in enumerate([auc_r,auc_m,auc_rm]): ax[1].text(i,v+0.005,f"{v:.3f}",ha="center",fontsize=9)
ax[1].set_ylim(0.5,1.0); ax[1].set_ylabel("subtype AUROC (5-fold CV)")
ax[1].set_title("Multi-omics: does methylation add to RNA?")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"meth_arm.png"),dpi=140); plt.close()
print("[meth] saved figure + metrics")
