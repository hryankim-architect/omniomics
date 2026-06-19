---
title: "Anchored multi-omics integration and knowledge-anchored residual discovery: a never-below-the-best-view integrator that finds what the textbook misses"
author: "H. Ryan Kim"
date: "2026"
---

**Running title:** Knowledge-anchored residual discovery

**Affiliation:** Independent (omniomics-prototype). Correspondence: ryan.kim2112@gmail.com.

**Keywords:** multi-omics integration; prior-informed machine learning; gene signatures; breast cancer; TCGA; METABRIC; immunotherapy biomarkers; interpretable discovery.

---

## Abstract

**Background.** Multi-omics integration is widely assumed to improve prediction, yet on real cancer data a
strong single modality is hard to beat and naive fusion often underperforms it. **Methods.** We define an
*anchored* integrator that anchors on the empirically strongest modality — or on a fixed, zero-parameter
textbook prior — and adds the remaining data only as a non-negative gated residual, so the model is provably
never below its anchor and improves on it only where another view carries orthogonal signal. Mining that
residual yields a discovery procedure that names anchor-orthogonal features, with a matched-random-panel
noise control and held-out / permutation verification. We evaluate on TCGA-BRCA (n up to 1,152) and validate
externally on METABRIC (n = 1,175–1,980). **Results.** Across nine subtype labels the gate adds nothing where
one modality dominates and engages only on a constructed positive control (+0.047 AUROC). Anchoring on a
fixed 20-gene proliferation signature (zero trained parameters) reaches AUROC 0.919 on Luminal A vs B —
nearly the trained 1,500-gene model (0.942) — and gating genome-wide data reaches 0.947, the only clean
super-additive gain observed. The residual of that prior recovers a basal/keratinization axis (Reactome
p ≈ 8×10⁻¹¹), verified by held-out replication (10/10 splits), stability (10/10), and a permuted-label null;
it **reproduces in METABRIC** (panel Δ +0.036; independent re-discovery overlaps 20/30, hypergeometric
p ≈ 7×10⁻²⁷). A second anchor (the ERBB2 amplicon for HER2) yields a distinct neuroendocrine/immune axis,
and an estrogen-receptor signature yields nothing (a specificity control). The HER2 axis does **not**
replicate in METABRIC because the amplicon anchor is near-complete there (0.997). The method is also not
breast-cancer- or expression-specific: on an NSCLC anti-PD-1 immunotherapy cohort (n = 227, mutation/clinical
features), anchored on the textbook biomarker tumour mutational burden, the residual independently recovers
the field's other established biomarkers — PD-L1 (orthogonal to TMB), EGFR and STK11 mutation (resistance)
(panel Δ +0.061 vs +0.006, p = 0.038). **Conclusion.** Anchor on
established biology and let a gate decide, honestly, whether genome-wide data beats it; where it does, the
residual names the new axis. External reproducibility tracks biological coherence.

## 1. Introduction

The prevailing premise that adding omics layers improves prediction is contradicted by careful benchmarks:
a large TCGA survival study found adding data types most often *hurt* performance, with mRNA (±miRNA)
sufficient for most cancers [1]; flexible late-fusion shows the dominant modality is task-specific and
pan-cancer multi-omics barely exceeds the best single view [2]; and naive concatenation routinely
underperforms structured fusion [3]. The honest objective is therefore not "always fuse" but to *never do
worse than the best single modality, and improve on it only where another modality carries signal the leader
cannot.* We formalise such an integrator, generalise its anchor from a data modality to established
knowledge, and repurpose the construction as an interpretable discovery engine — then ask, in an independent
cohort, whether its discoveries reproduce.

## 2. Results

**An anchored gate is never below the best single view.** On Luminal A vs B, RNA alone reaches AUROC 0.947
and methylation 0.745; a naive RNA+methylation stack drops to 0.941 (below RNA), whereas the gated combiner
holds at 0.947 (weight β = 0 in 9/10 repeats). A constructed positive control — a held-out methylation set
the RNA cannot see — engages the gate (β = 4–8) for a significant +0.047. Run blind on nine expert-defined
subtype labels, anchor selection routes to methylation for all five methylation-defined clusters and to RNA
for all four PAM50 calls, never below the leader: the method follows the biology rather than a built-in
preference.

**A zero-parameter textbook prior is a strong anchor, and yields the only clean fusion gain (Figure 1A).**
Luminal A vs B is, by textbook, a proliferation distinction. A 20-gene proliferation index with zero trained
parameters reaches AUROC 0.919 — close to the fully trained 1,500-gene model (0.942) — and gating
genome-wide data onto this fixed prior reaches 0.947 (Δ +0.029), exceeding the pure data model. With the
Horvath epigenetic clock as the anchor for normal-tissue age, the clock alone (0.947) already beats RNA
(0.911) and data adds ≈ 0: the textbook suffices and the gate says so.

**The residual names a verified new axis.** Mining the proliferation-anchored residual surfaces the basal /
squamous-lineage axis (KRT5/14/17/6B, TP63, DSG3/DSC3, SOX10, COL17A1, KLK5/7/8; partial r ≈ 0.46). It is
verified three ways — a train-selected panel beats random panels on a held-out test in 10/10 splits (no
selection leakage), the basal core recurs in 10/10 splits, and under permuted labels the advantage collapses
— and is enriched for cornified-envelope formation / keratinization / epidermis development (Reactome
p ≈ 8×10⁻¹¹).

**Generalisation and specificity (Figure 1A).** Anchoring HER2 status on the ERBB2 amplicon (incomplete in
TCGA, AUROC 0.752) discovers a distinct, verified neuroendocrine/secretory + immune axis (Δ +0.054; panel
vs random p = 0.024; held-out 8/8). Anchoring ER status on a textbook ER/luminal signature (complete, 0.938)
discovers nothing (Δ −0.001): a real-data specificity control showing the method does not manufacture axes.

**External validation (Figure 1B).** In the independent METABRIC cohort the basal axis reproduces: the fixed
TCGA panel adds Δ +0.036 over the proliferation prior (combined 0.960; beats random panels, p = 0.048), and
an unbiased re-discovery independently recovers the same basal genes, overlapping the TCGA panel 20/30
(hypergeometric p ≈ 7×10⁻²⁷). The HER2 axis does **not** reproduce — there the amplicon anchor is
near-complete (0.997 vs 0.752 in TCGA), leaving no residual (panel Δ ≈ 0; overlap 2/30). External
reproducibility tracks biological coherence: the pathway-enriched basal axis reproduces, while the diffuse
HER2 axis was a cohort-specific residual.

**Cross-domain generalisation (a different cancer, question, and feature type).** The method is not specific
to breast cancer or to expression. On NSCLC patients receiving anti-PD-1 checkpoint blockade (Hellmann/MSK
2018, n = 227; mutation and clinical features; endpoint = durable clinical benefit), anchoring on the
textbook immuno-oncology biomarker tumour mutational burden (TMB; anchor AUROC 0.60) and mining the residual
independently recovers the field's *other* established biomarkers: PD-L1 score (positive partial correlation,
and genuinely orthogonal to TMB, corr 0.00), EGFR mutation and STK11/LKB1 mutation (both negative — known
checkpoint-blockade resistance). The discovered panel adds Δ +0.061 versus +0.006 for matched random panels
(p = 0.038). Anchored on the textbook biomarker, the procedure rediscovers the known complementary biomarkers
of a different disease.

## 3. Discussion

The method instantiates and unifies three established traditions — clinical-offset / incremental-value
modelling (an established model as offset, omics on the residual; the closest methodological twin) [4];
biologically-informed models that bake known biology into structure [5,6]; and machine learning on
established gene signatures [7]. Its specific contribution is a packaged, never-below-anchor, margin-gated
residual on a *zero-parameter* prior, plus a residual-discovery framing with an explicit noise control, and
a demonstration that the discovered axes' external reproducibility is predictable from their biological
coherence. The honest scope: predictive gains over a strong anchor are modest — the value is *routing and
discovery*, not large AUROC wins — consistent with the rarity of genuine super-additive multi-omics gains
[1,2]. Discovered axes are candidate hypotheses; the basal axis is externally validated, the HER2 axis is
cohort-specific, and other endpoints await further cohorts.

## 4. Methods

**Anchored integration.** `select_anchor` scores each modality by repeated stratified-CV AUROC (tie-broken by
robustness) and picks the top composite, inside the outer CV. `anchored_integrate` pins the anchor and adds a
secondary on its residual as `logit(anchor) + β·secondary`, with β ≥ 0 chosen on held-out data (β = 0
allowed; a margin gate plus repeated inner CV control small-sample noise). `forward_integrate` / `auto_integrate`
extend this to any number of modalities by greedy forward selection. **Knowledge anchor.** `signature_score`
builds a fixed mean-z score from a curated gene set; `knowledge_anchored_integrate` pins it as the anchor.
**Residual discovery.** `anchored_residual_discovery` ranks features by partial correlation with the label
controlling for the anchor, retains those orthogonal to the anchor (|r| < 0.6), gates the panel, and tests it
against matched random panels and (for verification) held-out splits and permuted labels. Data: TCGA-BRCA
(Xena HiSeqV2, HumanMethylation450, clinical) and METABRIC (microarray, clinical). Pathway enrichment:
g:Profiler. All functions are in `omniomics.multiomics`; runners and recorded metrics are in the repository.

## 5. Data and Code Availability

Code, runners, recorded metrics, unit tests and continuous-integration guards: the `omniomics` repository
(`reports/dmoi_*.py`, `*_results.csv`, `tests/`). Full methods note: `reports/anchored_integration_methods.md`.
Data are public (TCGA via UCSC Xena; METABRIC via cBioPortal).

## Figure

![**Figure 1. Knowledge-anchored discovery and its external validation.** **(A)** Across three TCGA-BRCA
endpoints, the gain in AUROC from adding genome-wide data beyond a zero-parameter textbook anchor: the data
adds a real axis where the textbook is incomplete (LumA/B proliferation → basal, +0.029; HER2 amplicon →
neuroendocrine/immune, +0.054) and nothing where it is complete (ER signature, −0.001; specificity).
**(B)** In the independent METABRIC cohort, an unbiased re-discovery recovers the basal axis (20/30 overlap
with the TCGA panel; direct replication Δ +0.036, p = 0.048) but not the HER2 axis (2/30), because the
amplicon anchor is already complete in METABRIC.](figs/discovery_summary.png)

## References

1. Li et al. Does combining numerous data types improve or hinder survival prediction? *BMC Med Inform Decis Mak* (2024).
2. Nikolaou et al. Flexible late-fusion for multi-omics survival. *Cancer Res* (2023).
3. Makrodimitris et al. Linear vs non-linear joint embedding for multi-omics. *Brief Bioinform* (2023).
4. Volkmann et al. A plea for taking all available clinical information into account when assessing the predictive value of omics data. *BMC Med Res Methodol* (2019).
5. Elmarakeby et al. Biologically informed deep neural network for prostate cancer discovery (P-NET). *Nature* (2021).
6. Liu et al. Pathformer: a pathway-informed transformer for multi-omics. *Bioinformatics* (2024).
7. Sarafidis et al. A machine-learning model on established signatures for early breast cancer. *ESMO Open* (2024).
8. Ding, Li, Narasimhan & Tibshirani. Cooperative learning for multiview analysis. *PNAS* (2022).
9. Horvath. DNA methylation age of human tissues and cell types. *Genome Biology* (2013).
