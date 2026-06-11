"""DNA-methylation arm: aggregate HM450 promoter probes to gene-level methylation.
Promoter probes = within +/-1500 bp of a gene TSS (precomputed probe->gene table)."""
import numpy as np, pandas as pd

def aggregate(M, probe_gene_tsv, chunk=40000):
    """Aggregate an already-loaded probe x sample beta matrix M to gene x sample
    mean methylation, using a (probe, gene) table. Reuse across genomic contexts."""
    pg = pd.read_csv(probe_gene_tsv, sep="\t", header=None, names=["probe","gene"])
    pg = pg[pg["probe"].isin(M.index)]
    genes = sorted(pg["gene"].unique()); gcode = {g:i for i,g in enumerate(genes)}
    pidx = {p:i for i,p in enumerate(M.index)}; Mv = M.values; nS = Mv.shape[1]
    gsum = np.zeros((len(genes), nS)); gcnt = np.zeros((len(genes), nS))
    pr = pg["probe"].map(pidx).values; gc = pg["gene"].map(gcode).values
    for s in range(0, len(pr), chunk):
        block = Mv[pr[s:s+chunk]]; gi = gc[s:s+chunk]
        mask = ~np.isnan(block)
        np.add.at(gsum, gi, np.where(mask, block, 0.0)); np.add.at(gcnt, gi, mask.astype(float))
    gm = np.divide(gsum, gcnt, out=np.full_like(gsum, np.nan), where=gcnt>0)
    return pd.DataFrame(gm, index=[g.upper() for g in genes], columns=M.columns).groupby(level=0).mean()

def load_hm450_promoter(meth_tsv, probe_gene_tsv, dtype="float32", chunk=40000):
    """meth_tsv: promoter-probe x sample beta matrix (probe in col 0).
    probe_gene_tsv: 2-col (probe, gene). Returns gene x sample mean-promoter-methylation."""
    M = pd.read_csv(meth_tsv, sep="\t", index_col=0)
    M = M.astype(dtype)
    pg = pd.read_csv(probe_gene_tsv, sep="\t", header=None, names=["probe","gene"])
    pg = pg[pg["probe"].isin(M.index)]
    genes = sorted(pg["gene"].unique())
    gcode = {g:i for i,g in enumerate(genes)}
    pidx = {p:i for i,p in enumerate(M.index)}
    Mv = M.values
    nS = Mv.shape[1]; gsum = np.zeros((len(genes), nS), dtype="float64"); gcnt = np.zeros(len(genes))
    pr = pg["probe"].map(pidx).values; gc = pg["gene"].map(gcode).values
    for s in range(0, len(pr), chunk):
        pi = pr[s:s+chunk]; gi = gc[s:s+chunk]
        block = Mv[pi]                          # pairs x samples
        mask = ~np.isnan(block)
        block0 = np.where(mask, block, 0.0)
        np.add.at(gsum, gi, block0)
        np.add.at(gcnt, gi, mask.sum(axis=1)/nS)  # approx per-gene probe count
        # proper per-sample counts:
        cnt_block = mask.astype("float64")
        # accumulate per-sample valid counts
        if s == 0:
            gcnt_s = np.zeros((len(genes), nS))
        np.add.at(gcnt_s, gi, cnt_block)
    gmeth = np.divide(gsum, gcnt_s, out=np.full_like(gsum, np.nan), where=gcnt_s>0)
    out = pd.DataFrame(gmeth, index=[g.upper() for g in genes], columns=M.columns)
    return out.groupby(level=0).mean()
