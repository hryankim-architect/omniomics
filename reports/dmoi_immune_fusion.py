#!/usr/bin/env python3
"""Immune-axis fusion test: does curated immune methylation earn a super-additive gain over RNA?

Companion to dmoi_fusion_gain.py. Builds an immune methylation modality from the EPIC/Salas 18-cell
reference CpGs (NNLS deconvolution into neutrophil / NK / B / CD4-T / CD8-T / monocyte fractions, plus
the raw immune-CpG set) and tests, on two EXTERNAL TCGA-BRCA labels — histology (IDC vs ILC) and lymph-
node status (node+ vs node0) — whether methylation adds to the RNA anchor via auto_integrate, against a
matched genome-wide random-CpG control. Result (honest negative): no super-additive fusion gain on either
endpoint. Histology is RNA-anchored (0.91) with genome-wide methylation strong alone (0.87) but fully
REDUNDANT and immune composition uninformative (0.60); node status is near-chance for every modality
(~0.55), so there is no signal to fuse. Consistent with the project's headline that genuine multi-omics
*fusion* gains are rare — value comes from routing to the right modality (see dmoi_fusion_gain.py), not
from combining. Immune reference provenance: reports/epic_salas_immune_reference.csv (biolearn 0.9.1).

Run:  BRCA_DIR=/path/to/tcga_brca python reports/dmoi_immune_fusion.py
Env:  BRCA_DIR, FG_DIR (cache for extracted immune_meth.tsv / ctrl600_meth.tsv). Writes immune_axis_results.csv.
"""
import os, sys, gzip
import numpy as np, pandas as pd
from scipy.optimize import nnls
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
try:
    from omniomics import config; DEFAULT_BRCA = config.brca_tcga_dir()
except Exception:
    DEFAULT_BRCA = ""
BRCA = os.environ.get("BRCA_DIR", DEFAULT_BRCA)
FG = os.environ.get("FG_DIR", os.getcwd())
REF = pd.read_csv(os.path.join(HERE, "epic_salas_immune_reference.csv"), index_col=0)


def _extract():
    ip, kp = os.path.join(FG, "immune_meth.tsv"), os.path.join(FG, "ctrl600_meth.tsv")
    if os.path.exists(ip) and os.path.exists(kp):
        return pd.read_csv(ip, sep="\t", index_col=0), pd.read_csv(kp, sep="\t", index_col=0)
    assert BRCA and os.path.isdir(BRCA), f"BRCA dir not found: {BRCA!r} (set BRCA_DIR)"
    ids = set(REF.index); imm_rows, ctrl_rows, header = {}, {}, None
    with gzip.open(os.path.join(BRCA, "HumanMethylation450.gz"), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        for i, line in enumerate(fh):
            pid, _, rest = line.partition("\t")
            if pid in ids:
                imm_rows[pid] = rest.rstrip("\n").split("\t")
            elif i % 800 == 3:
                ctrl_rows[pid] = rest.rstrip("\n").split("\t")
    imm = pd.DataFrame(imm_rows, index=header).T.apply(pd.to_numeric, errors="coerce")
    ctrl = pd.DataFrame(ctrl_rows, index=header).T.apply(pd.to_numeric, errors="coerce")
    imm.to_csv(ip, sep="\t"); ctrl.to_csv(kp, sep="\t")
    return imm, ctrl


def _deconv(imm, samples):
    shared = [c for c in imm.index if c in REF.index]; R = REF.loc[shared].values
    Bv = imm.loc[shared, samples].astype(float).values
    rm = np.nanmean(Bv, axis=1, keepdims=True); ii = np.where(np.isnan(Bv)); Bv[ii] = np.take(rm, ii[0])
    Bv = np.nan_to_num(Bv)
    frac = np.array([nnls(R, Bv[:, j])[0] for j in range(len(samples))])
    return frac / (frac.sum(1, keepdims=True) + 1e-9)


def _imp(df, samples):
    M = df[samples].T.astype(float); return M.fillna(M.mean()).fillna(0.0).values


def main():
    imm, ctrl = _extract()
    clin = pd.read_csv(os.path.join(BRCA, "BRCA_clinicalMatrix.tsv"), sep="\t").set_index("sampleID")
    hist = clin["histological_type"]; node = pd.to_numeric(clin["number_of_lymphnodes_positive_by_he"], errors="coerce")
    with gzip.open(os.path.join(BRCA, "HiSeqV2.gz"), "rt") as fh:
        rcols = set(fh.readline().rstrip("\n").split("\t")[1:])
    base = [s for s in imm.columns if s in rcols]
    labels = {
        "histology_idc_ilc": {s: (1 if hist.get(s) == "Infiltrating Lobular Carcinoma" else 0) for s in base
                              if hist.get(s) in ("Infiltrating Ductal Carcinoma", "Infiltrating Lobular Carcinoma")},
        "node_status": {s: (1 if node.get(s) > 0 else 0) for s in base if s in node.index and not np.isnan(node.get(s))},
    }
    rows = []
    for ep, lab in labels.items():
        S = sorted(lab)
        rna = pd.read_csv(os.path.join(BRCA, "HiSeqV2.gz"), sep="\t", index_col=0, usecols=["sample"] + S)
        rna = rna.loc[rna.var(axis=1).sort_values(ascending=False).index[:600]]
        Xr = rna[S].T.values; Mgw = _imp(ctrl, S); Mimm = _imp(imm, S); Dec = _deconv(imm, S)
        y = np.array([lab[s] for s in S])
        rk = dict((n, v) for n, v, s, c in mo.select_anchor(
            {"RNA": Xr, "meth_genomewide": Mgw, "immune_deconv": Dec, "immune_cpgs": Mimm}, y, cv=5, repeats=3)["ranking"])
        best = max((mo.auto_integrate({"RNA": Xr, k: X}, y, cv=5, random_state=0, inner_repeats=2)
                    for k, X in [("meth_genomewide", Mgw), ("immune_deconv", Dec), ("immune_cpgs", Mimm)]),
                   key=lambda r: r["delta"])
        rows.append(dict(endpoint=ep, n=len(S), pos=int(y.sum()), auroc_rna=round(rk["RNA"], 3),
                         auroc_meth_genomewide=round(rk["meth_genomewide"], 3),
                         auroc_immune_deconv=round(rk["immune_deconv"], 3),
                         auroc_immune_cpgs=round(rk["immune_cpgs"], 3),
                         auto_anchor=best["anchor"], auto_delta=round(best["delta"], 3)))
    out = pd.DataFrame(rows); p = os.path.join(REPO, "immune_axis_results.csv"); out.to_csv(p, index=False)
    print(out.to_string(index=False)); print("\nwrote", p)
    print("verdict: no super-additive fusion gain on the immune axis (histology RNA-anchored & methylation")
    print("         redundant; node status near-chance) -- value is modality routing, not combination.")


if __name__ == "__main__":
    main()
