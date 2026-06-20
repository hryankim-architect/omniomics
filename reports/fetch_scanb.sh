#!/usr/bin/env bash
# Reproducibly build the SCAN-B (GSE96058) inputs for dmoi_endpoint_panel.py without committing large data.
# Downloads the public GSE96058 gene-expression matrix (~592 MB) + per-platform series matrices (phenotype),
# then extracts ONLY the panel marker-gene rows (-> scanb_markers.csv) and the pam50/er/her2 labels
# (-> scanb_pheno.csv). Output goes to SCANB_DIR (default: <BRCA_DIR>/../scanb).
#
# Usage:  SCANB_DIR=/path/to/scanb bash reports/fetch_scanb.sh
set -euo pipefail
SD="${SCANB_DIR:-./scanb}"; mkdir -p "$SD"; cd "$SD"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058"

echo "== download (resumable) =="
curl -sSL -C - -o gene_expression.csv.gz "$BASE/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
curl -sSL -o GSE96058-GPL11154_series_matrix.txt.gz "$BASE/matrix/GSE96058-GPL11154_series_matrix.txt.gz"
curl -sSL -o GSE96058-GPL18573_series_matrix.txt.gz "$BASE/matrix/GSE96058-GPL18573_series_matrix.txt.gz"

echo "== extract marker rows + phenotype =="
python3 - <<'PY'
import gzip, csv, collections, os
GENES=set('''MKI67 PCNA CCNB1 CCNB2 CDK1 AURKA AURKB BUB1 CCNE1 CDC20 TOP2A TYMS RRM2 UBE2C CENPF FOXM1 MELK KIF2C NUSAP1 PTTG1
ESR1 GATA3 FOXA1 XBP1 TFF1 PGR GREB1 CA12 SLC39A6 NAT1 AR MLPH
ERBB2 GRB7 STARD3 PGAP3 TCAP PNMT PSMD3 GSDMB ORMDL3
KRT5 KRT14 KRT17 KRT6B TP63 DSG3 DSC3 SOX10 COL17A1 FOXC1 MIA SFRP1
CD8A CD3D GZMB PRF1 IFNG CXCL9 CXCL10 GZMK NKG7 CCL5 CD2 PTPRC'''.split())
with gzip.open('gene_expression.csv.gz','rt') as fh, open('scanb_markers.csv','w',newline='') as out:
    out.write(fh.readline())                       # sample header
    n=0
    for line in fh:
        if line.split(',',1)[0].strip().strip('"') in GENES:
            out.write(line); n+=1
print('marker rows:', n)
def parse(path):
    titles=None; chars=[]
    with gzip.open(path,'rt') as fh:
        for line in fh:
            if line.startswith('!Sample_title'): titles=next(csv.reader([line.rstrip()],delimiter='\t'))[1:]
            elif line.startswith('!Sample_characteristics_ch1'): chars.append(next(csv.reader([line.rstrip()],delimiter='\t'))[1:])
    keymaps={}
    for col in chars:
        k=next((v.split(':',1)[0].strip() for v in col if ':' in v), None)
        if k in ('pam50 subtype','er status','her2 status'): keymaps.setdefault(k,col)
    out={}
    for i,t in enumerate(titles):
        out[t]={k:(col[i].split(':',1)[1].strip() if i<len(col) and ':' in col[i] else 'NA') for k,col in keymaps.items()}
    return out
d={}
for p in ('GSE96058-GPL11154_series_matrix.txt.gz','GSE96058-GPL18573_series_matrix.txt.gz'): d.update(parse(p))
with open('scanb_pheno.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['sample','pam50','er','her2'])
    for t,r in d.items(): w.writerow([t, r.get('pam50 subtype','NA'), r.get('er status','NA'), r.get('her2 status','NA')])
print('phenotype samples:', len(d), '| pam50:', dict(collections.Counter(v.get('pam50 subtype') for v in d.values())))
PY
echo "done -> $SD/scanb_markers.csv , $SD/scanb_pheno.csv"
