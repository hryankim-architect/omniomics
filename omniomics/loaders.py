"""Generic loaders that normalize heterogeneous processed files to a gene x sample matrix."""
import glob, os, re, gzip
import numpy as np, pandas as pd

def load_cufflinks_fpkm_dir(path, pattern="*.fpkm_tracking", sample_from=r"(Set\d)_([A-Z]+)_genes"):
    """GSE57577-style: one Cufflinks genes.fpkm_tracking per sample."""
    expr, names, loci = {}, {}, {}
    for f in sorted(glob.glob(os.path.join(path, pattern))):
        m = re.search(sample_from, os.path.basename(f))
        sid = "_".join(m.groups()[::-1]) if m else os.path.basename(f)
        df = pd.read_csv(f, sep="\t", usecols=["gene_id","gene_short_name","locus","FPKM"])
        g = df.drop_duplicates("gene_id").set_index("gene_id")
        names.update(g["gene_short_name"].to_dict()); loci.update(g["locus"].to_dict())
        expr[sid] = df.groupby("gene_id")["FPKM"].sum()
    mat = pd.DataFrame(expr).fillna(0)
    return mat, pd.Series(names), pd.Series(loci)

def load_matrix(path, gene_col=0, sep="\t", drop_cols=None, gene_is_symbol=True):
    """Generic processed expression matrix (genes x samples). Handles .gz."""
    op = gzip.open if path.endswith(".gz") else open
    df = pd.read_csv(path, sep=sep)
    gc = df.columns[gene_col]
    if drop_cols:
        df = df.drop(columns=[c for c in drop_cols if c in df.columns])
    df = df.rename(columns={gc: "gene"})
    df = df.groupby("gene").mean(numeric_only=True)
    return df  # index = gene (symbol), columns = samples

def gene_symbol_index(mat, names=None):
    """Return a copy indexed by uppercase gene symbol (for cross-study joins)."""
    out = mat.copy()
    if names is not None:
        out.index = [str(names.get(i, i)).upper() for i in out.index]
    else:
        out.index = [str(i).upper() for i in out.index]
    return out.groupby(level=0).mean()
