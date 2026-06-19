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

## Anchored multi-omics integration & knowledge-anchored discovery

A second contribution distilled from this project: a multi-omics integrator that is **never worse than its
best single view**, generalised into a **discovery** tool. Full write-up: `reports/anchored_integration_methods.md`;
standalone preprint: `reports/anchored_integration_preprint.pdf`.

**Idea.** Anchor on the empirically-strongest modality *or a fixed textbook prior* (zero trained parameters)
and add the other data only as a **non-negative gated residual** — never below the anchor, improving on it
only where another view carries orthogonal signal. Mining that residual turns the integrator into a
discovery engine that names the features beating the textbook, with a matched-random-panel noise control.

**Library — `omniomics.multiomics`:**
```python
select_anchor(modalities, y)                    # pick the anchor empirically, per task, by CV
anchored_integrate(Xa, Xs, y)                   # gate a secondary onto the anchor's residual (never below)
auto_integrate(modalities, y)                   # any number of modalities, one call (forward gating)
signature_score(expr, genes)                    # fixed textbook signature -> zero-param anchor score
knowledge_anchored_integrate(score, mods, y)    # anchor on established knowledge, gate data on top
anchored_residual_discovery(score, X, names, y) # discover anchor-orthogonal axes beyond the textbook
```

**Validated results (TCGA-BRCA + METABRIC):**

| result | finding |
|---|---|
| anchored integration | never below the best single modality; gate engages only on a positive control (+0.047) |
| anchor selection (9 labels) | methylation for methylation-defined / RNA for expression-defined labels — not RNA-biased |
| knowledge anchor (LumA/B) | 20-gene proliferation prior (0 params) AUROC 0.919; + data 0.947 — first clean fusion gain |
| residual discovery (LumA/B) | basal/keratinization axis beyond proliferation — verified (held-out, stable, label-specific; enriched p≈8e-11) |
| generalization | HER2 → neuroendocrine/immune axis (verified); ER → nothing (textbook complete: a specificity control) |
| **external validation (METABRIC)** | basal axis **reproduces** (Δ+0.036; independent re-discovery 20/30, p≈7e-27); HER2 cohort-specific (amplicon already complete there) |
| **cross-domain (NSCLC anti-PD1)** | anchor on textbook **TMB** → residual recovers the other IO biomarkers: PD-L1 (orthogonal +), EGFR & STK11 (resistance); Δ+0.061 vs random +0.006, p=0.038 — not breast/expression-specific |
| **cross-cancer (TCGA lung)** | the breast basal/keratinization axis **re-discovers** in lung (LUAD vs LUSC, n=1129): residual overlaps the breast panel 10/30 (KRT5/14/6B, TP63, DSG3/DSC3…), hypergeometric p≈3e-16 — the discovered biology, not a cohort artefact |

**Reproduce:**
```bash
python run_auto_integrate.py                  # anchored integration on TCGA-BRCA
python reports/dmoi_knowledge_anchor.py       # textbook-prior anchoring (proliferation / Horvath clock)
python reports/dmoi_residual_discovery.py     # the basal discovery + 3-way verification
python reports/dmoi_discovery_er.py           # generalization: ER (specificity negative)
python reports/dmoi_discovery_her2.py         # generalization: HER2 (2nd positive axis)
python reports/dmoi_external_metabric.py      # external validation: basal reproduces in METABRIC
python reports/dmoi_external_her2_metabric.py # external validation: HER2 cohort-specific
NSCLC_TABLE=.../patient_table.csv python reports/dmoi_discovery_nsclc_io.py  # cross-domain: NSCLC anti-PD1 (TMB anchor)
python reports/dmoi_external_lung.py          # cross-cancer: basal axis re-discovers in TCGA lung (downloads Xena)
```
Recorded metrics: `{auto_integrate,external_subtype,knowledge_anchor,discovery,discovery_{er,her2},
fusion_gain,immune_axis,discordance_test,external_validation{,_her2}_metabric}_results.csv`; CI guards in
`tests/test_golden.py`, unit tests in `tests/test_dmoi_v2.py`.

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
