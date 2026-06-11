# omniomics

[![golden](https://github.com/hryankim-architect/omniomics/actions/workflows/ci.yml/badge.svg)](https://github.com/hryankim-architect/omniomics/actions/workflows/ci.yml)
![python](https://img.shields.io/badge/python-3.10%E2%80%933.12-blue)
![license](https://img.shields.io/badge/license-MIT-green)

A reusable, manifest-driven **multi-omics analysis engine** — distilled from a full reproduction of
Noh et al. (*Molecular Cell* 2015, GSE57577) and scaled to cross-study, cross-cohort, cross-platform
and multi-assay data, with a provenance-gated golden-task suite and a swarm drop-in.

```bash
pip install -e .                 # editable install (numpy/pandas/scipy/sklearn/matplotlib)
python run_golden.py             # GSE57577 reproduction — known-answer regression gate (CI runs this)
python run_gse57577_dmoi.py      # reproduce the paper's graphical-abstract pattern with DMOI
```
For the breast-cancer (TCGA/METABRIC) pipeline, set `DMOI_BRCA_DATA` and run `omniomics-prepare-brca`
once — see **`LOCAL_SETUP.md`**. Methods + every result are in **`PROTOTYPE_REPORT.md`**;
the scaling rationale in `SCALING_RESEARCH_ROADMAP.md`.

## Layout
```
omniomics/            reusable engine
  geo.py              GEO supplementary download by accession
  loaders.py          normalize heterogeneous processed files -> gene x sample matrix
  expression.py       paired empirical-Bayes moderated DE (valid at n=2) + ncRNA filter
  harmonize.py        common gene panel, quantile norm, ComBat-lite, PCA, batch-variance score
manifest.yaml         dataset registry (add a study here -> engine ingests it)
run_golden.py         regression gate: reproduce verified GSE57577 numbers, unattended
run_cohort.py         Phase 1: harmonize GSE57577 + GSE77003, correct batch, meta-analyze
```

## Run
```bash
python3 run_golden.py     # -> ALL PASS  (WWD=1888, R=9, TKO=3, targets 6/6, ChIP 1.64x)
python3 run_cohort.py     # -> mouse ESC cross-study harmonization + reproducibility
python3 run_cohort_brca.py # -> TCGA vs METABRIC cross-platform harmonization (user's data)
python3 run_meth_arm.py   # -> RNA + HM450 methylation multi-omics on matched TCGA samples
python3 run_joint_dmoi.py # -> MOFA-style joint embedding + DMOI-lite pathway fusion verdict
python3 run_lumab_dmoi.py # -> LumA/LumB pole-conditioned multi-omics (validates compact-rep thesis)
python3 run_meth_context.py # -> methylation by genomic context (promoter/shore/enhancer)
python3 run_dmoi_enhancer.py # -> enhancer methylation in DMOI fusion; repeated-CV significance
python3 run_combat_benchmark.py    # -> EB-ComBat vs naive: N>2 confounded-batch correction
python3 run_disagreement_regime.py # -> when multi-omics/disagreement earns value (RNA degradation)
python3 run_disagreement_axes.py   # -> multi-omics/disagreement across HER2/LumAB/Basal axes
python3 run_mouse_n3_combat.py     # -> EB-ComBat on 3 real heterogeneous mouse studies (N>2)
python3 run_her2_redesign.py       # -> HER2 amplicon prior redesign (selective methylation routing)
# swarm/ : drop-in /golden command + omniomics-runner agent (see swarm/SWARM_WIRING.md)
python3 golden/run_golden_brca.py --use-cached  # -> golden-task regression gate + audit chain
```

## Results
- **Golden task:** engine reproduces every checked GSE57577 value from the manifest, unattended.
- **Harmonization:** 21,017 common gene symbols across two independent labs (FPKM vs RPKM); study
  batch effect on PC1 drops from **0.94 → 0.00** after ComBat-lite.
- **Cross-study biology:** the shared "restore Dnmt3a2 to TKO ESCs" axis reproduces — **70%
  sign-concordance** on robustly-changed genes (Spearman 0.20 genome-wide).

See `PROTOTYPE_REPORT.md` for full write-up and caveats.
