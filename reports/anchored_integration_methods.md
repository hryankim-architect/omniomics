# Anchored multi-omics integration — methods note

*A small, interpretable, never-below-the-best-single-modality integrator, with the experiments that
justify it. Self-contained; the API lives in `omniomics/multiomics.py`.*

## 1. Problem

Multi-omics fusion is usually sold as "more modalities → better predictions." In practice, when one
modality dominates a task, symmetric fusion (concatenate everything, fit one model) often performs
*worse* than that single modality, because a weak, noisy layer dilutes a strong one. The honest goal
is therefore not "always fuse" but: **never do worse than the best single modality, and improve on it
only where another modality carries signal the leader cannot see.** This note defines an integrator
built to that contract and reports what it finds on real cancer data — including the uncomfortable
finding that genuine fusion gains are rarer than the literature implies.

## 2. Method

Three composable pieces, all in `omniomics.multiomics`:

**(a) Anchor selection — `select_anchor`.** Score each modality by repeated stratified-CV predictive
performance (AUROC / C-index), the primary criterion; tie-break by robustness (subtract a multiple of
the CV standard deviation, weight by coverage/measurement reliability). The anchor is the top
composite. The dominant modality is task-specific, so this is always chosen empirically, per task,
inside the outer CV — never fixed a priori.

**(b) Gated residual combiner — `anchored_gate` / `anchored_integrate`.** Pin the anchor; add each
secondary modality on the anchor's *residual* with a non-negative weight β chosen on held-out data,
with **β = 0 allowed**. The combined score is `logit(anchor) + β·secondary`. β = 0 whenever the
secondary is redundant (so combined ≡ anchor), and β > 0 only where the secondary earns a margin gain.
This is a constrained, interpretable special case of cooperative learning (which chooses the *degree*
of fusion adaptively by CV; Ding et al., PNAS 2022).

**(c) Forward gating for ≥ 3 modalities — `forward_integrate` / `auto_integrate`.** Do not concatenate
everything at once. Add modalities one at a time onto the current model's residual, in `select_anchor`
rank order, keeping each only if its gated contribution clears the margin — forward selection over
modalities, with an interpretable record (`added`) of which entered and at what β. `auto_integrate`
is the one-call entry point: `select_anchor` → forward-gated integration, any number of modalities.

**Small-sample robustness.** A secondary is added only if a non-zero β beats the anchor by more than a
`gate_margin`, and each β's inner-CV AUROC is averaged over `inner_repeats` splits. Without these
guards an n = 84 task spuriously engaged the secondary and fell 0.037 *below* the anchor; with them it
correctly defaults to the anchor (Δ = +0.000).

**(d) Knowledge anchor — `knowledge_anchored_integrate` / `signature_score`.** The anchor need not be a
*data* modality at all. Anchoring instead on a **fixed external biological prior** — a textbook gene
signature, a published clock, a known driver score, a clinical marker — carries **zero trained
parameters** (only an in-fold 1-parameter calibration), so it is maximally stable in small samples and
fully interpretable, and the gate then answers the clinically meaningful question *"does the genome-wide
data add anything beyond established biology?"*. This generalises `auto_integrate` from "anchor = best
empirical modality" to "anchor = established knowledge"; mechanically it is `forward_integrate` with the
prior pinned as the anchor. It is also where the programme finally produces a **clean positive fusion
gain** (see §4).

## 3. Why not a "disagreement" feature (the theory that redirected the design)

An earlier design (DMOI v1) added a hand-built RNA–methylation *difference* feature, `z_RNA − z_meth`,
as the cross-omics signal. That feature is a linear combination of the two main effects, so for any
linear classifier it is **provably redundant** — it lies in the span of the inputs and cannot change
the decision boundary. This was confirmed empirically: adding it moved AUROC/R² by Δ = 0.0000. The
real cross-omics signal is an *interaction* (`z_RNA × z_meth`), accessible only to a nonlinear learner:
a gradient-boosted model with explicit interaction terms beat the linear v1 by +0.009 (Wilcoxon
p < 0.05), with a permuted-label null collapsing to ~0.5. The lesson — fuse via interactions and a
nonlinear learner, or via the gated residual above; do not fake it with a linear difference — is what
motivated the anchored frame over symmetric fusion.

## 4. Results on TCGA-BRCA

**Symmetric fusion underperforms; the gate fixes it (LumA/B, RNA-defined).** RNA alone 0.947;
methylation alone 0.745; naive RNA+methylation stack 0.941 (−0.003, *below* RNA); RNA-anchored gated
combiner 0.947 (β = 0 in 9/10 repeats). The gate removes the symmetric-fusion penalty.

**Knowledge anchor — a zero-parameter textbook prior, and the first clean positive fusion gain.** LumA
vs LumB is, by textbook (St Gallen / Ki67), a *proliferation* distinction. A fixed 20-gene proliferation
index with **zero trained parameters** reaches AUROC **0.919** — nearly the fully trained 1500-gene RNA
model (0.942). Gating genome-wide data onto this fixed prior with `knowledge_anchored_integrate` reaches
**0.947 (Δ = +0.029)** — the gate engages and the combined **beats the pure data model**: the first clean
super-additive gain in the programme, unlocked precisely because the anchor is now a prior rather than the
strongest data view, so the data has orthogonal signal to add. The complementary case confirms the gate's
discipline: with the **Horvath clock** (353 CpGs, zero trained parameters) as the anchor for normal-tissue
age, the clock alone scores 0.947 — already above RNA's 0.911 — and gating RNA adds ≈ 0 (Δ = −0.003, noise):
the textbook suffices and the gate says so. Packaged as `knowledge_anchored_integrate` + `signature_score`;
runner `reports/dmoi_knowledge_anchor.py`; `knowledge_anchor_results.csv`.

**Knowledge-anchored residual discovery — separate the known, keep the real new.** The capstone combines
the whole toolkit (`anchored_residual_discovery`): anchor on the textbook prior, gate the data onto its
residual, test that the residual is signal-not-noise, and surface the features driving it that are
*orthogonal* to the anchor — candidate biology *beyond* the textbook. On LumA/B (anchor = proliferation),
the top anchor-orthogonal genes are unambiguous and coherent: **KRT5/14/17/6B, TP63, DSG3/DSC3, SOX10,
COL17A1, KLK5/7/8** — the classic **basal / squamous-lineage axis** (partial correlations ≈ 0.46 controlling
for proliferation). So beyond the textbook "LumA/B = proliferation" the data recovers a real *lineage* axis;
the discovered 30-gene panel gated onto the prior reaches 0.963 and beats matched random panels of the same
size (Δ +0.045 vs +0.027; only 1/40 random panels match it, p ≈ 0.05). Honest scope: the *predictive gain*
is modest (the candidate pool is itself subtype-enriched, so random panels also add ~+0.027 and the
label-permutation null on the raw gain is non-significant) — the method's value is **interpretable
discovery** (it reliably isolates a coherent, anchor-orthogonal new axis) more than a large AUROC win.
Runner `reports/dmoi_residual_discovery.py`; `discovery_results.csv` + `novel_genes.csv`.

*Verified three ways* (`discovery_verification.csv`): (V1) selecting the panel on a TRAIN split and scoring
it on a held-out TEST still beats random panels in 10/10 splits (Δ +0.051 vs +0.026) — so the gain is not
selection leakage; (V2) the basal core (KRT5/6B/14/17, DSG3, DSC3, TP63, SOX10, COL17A1, KLK7/8, CLDN19,
TRIM29) recurs in 10/10 train splits; (V3) under permuted labels the selected-vs-random advantage collapses
(real +0.021 > all five permutations), so the method finds the axis only when the signal is real.

**The gate is capable, not just safe (positive control).** On a methylation-defined endpoint (mean of
a held-out CpG set RNA cannot see), a disjoint methylation set scores 0.983 vs RNA 0.795; the gate
engages strongly (β = 4–8) for a significant **+0.047** over the anchor. So the gate captures real
orthogonal signal when it exists.

**External validation on nine real subtype labels — not RNA-biased (n = 491).** Run blind on five
expert-defined methylation clusters and four PAM50 expression calls (one-vs-rest):

| label family | endpoints | anchor chosen | leader vs other |
|---|---|---|---|
| methylation clusters 1–5 (methylation-defined) | 5/5 | **methylation** | meth 0.70–0.97 vs RNA 0.63–0.91 |
| PAM50 Basal/LumA/LumB/Her2 (expression-defined) | 4/4 | **RNA** | RNA 0.92–0.99 vs meth 0.82–0.99 |

Methylation wins and is selected on every methylation-defined cluster; RNA wins and is selected on
every PAM50 call; the gate never falls below the leader on any of the nine. This is the first real
(non-synthetic) evidence that anchor selection follows the biology rather than a built-in RNA
preference.

**Where multi-omics does *not* help (honest negative).** Patient age is the textbook methylation case
(the epigenetic clock), yet a 1,500-gene RNA model still edges out a random-3,000-CpG methylation model
on both tumour (0.799 vs 0.773) and normal-adjacent tissue (0.89 vs 0.81), and the gate adds nothing.
A random genome-wide CpG slice is not the 353 curated Horvath clock CpGs; demonstrating methylation
dominance needs *curated* biomarkers, not a flaw in the integrator.

**Where it does — the first real win for methylation (prediction confirmed).** Swapping the random slice
for the externally trained **Horvath (2013) 353-CpG clock** flips the normal-tissue result: the clock
methylation reaches AUROC **0.941 vs RNA 0.908, winning 10/10 CV seeds** (a matched random-CpG set scores
only 0.68, so the gain is the *curation*), and `auto_integrate` **selects methylation as the anchor** —
the first real, non-synthetic endpoint where methylation is the superior modality, and a direct refutation
of any "RNA-biased" reading. Even here there is no super-additive *fusion* gain (RNA adds ≈0 on the clock
residual): one modality dominates, the integrator picks it, and the honest theme holds — multi-omics value
comes from *routing to the right curated modality*, not from blindly combining. Runner
`reports/dmoi_fusion_gain.py`; data `fusion_gain_results.csv`.

The immune axis was then tested for the harder prize — a *super-additive* fusion gain — using a curated
immune methylation modality (EPIC/Salas reference CpGs, NNLS-deconvolved into immune-cell fractions) on
two external labels, histology (IDC vs ILC) and lymph-node status. None appeared: histology is RNA-anchored
with methylation strong-but-redundant, node status is near-chance for every modality. Across the whole
programme exactly one endpoint flips the leader to methylation (normal-tissue age) and none shows Δ > 0
fusion — confirming the thesis that multi-omics value here is modality *selection*, not combination. Runner
`reports/dmoi_immune_fusion.py`; data `immune_axis_results.csv`.

**Scorecard.** Across every real endpoint, RNA is simultaneously the most accurate, most
data-efficient, and most stable modality. The integrator's value is therefore **not** higher accuracy
on these tasks; it is (i) *safety* — never below the leader; (ii) *adaptivity* — picks methylation
when the biology is methylation-driven; (iii) *interpretability* — `added` records exactly which
modality contributed and how much. Used well, it is as much a tool to *avoid false multi-omics claims*
as to capture real ones.

## 5. API & reproducibility

```python
from omniomics import multiomics as mo
res = mo.auto_integrate({"RNA": Xr, "methylation": Xm, "CNV": Xc}, y, cv=5)
# res: anchor, auroc_anchor, auroc_combined, delta, added {mod:(folds,β)}, ranking
```

Primitives: `select_anchor`, `anchored_gate`, `anchored_integrate`, `forward_integrate`. Runners:
`run_auto_integrate.py` (TCGA-BRCA demo → `auto_integrate_results.csv`), `reports/dmoi_external_subtype.py`
(nine subtype labels → `external_subtype_results.csv`), `run_dmoi_v2.py` (interaction-vs-linear).
Tests: `tests/test_dmoi_v2.py` (12 unit tests) and `tests/test_golden.py` (recorded-metric CI guards,
including never-below-anchor and the anchor-selection direction on the nine real labels).

## 6. Limitations

The real-label validation uses labels each defined by their own modality (some recovery-of-self
circularity), and no endpoint here produced a fusion *gain* (Δ = 0 outside the constructed positive
control) — these experiments validate anchor selection and gate safety, not orthogonal-signal capture
in the wild. A truly external multi-omics cohort with a known complementary modality (e.g. methylation
deconvolution for immune content, or curated clock CpGs) is the natural next test. The method also
inherits the usual caveat: it cannot manufacture signal an assay never measured.

## 7. Verdict

Anchor on the leading modality, add the rest as a gated residual, and choose both the anchor and the
fusion degree by cross-validation. The result is an integrator that is never worse than its best single
view, improves on it exactly where a modality carries orthogonal signal, and tells you which modality
did the work — which, on honest real-data evaluation, is the property that actually matters.

## 8. Related work — the multi-omics-vs-single-modality debate

The central finding here — that one modality usually dominates a task, naive fusion rarely beats it, and
value comes from *routing to the right modality* rather than blending — sits within an active, mixed
literature. We did not invent the phenomenon; we reproduce scattered prior observations and package them
into a never-below-anchor tool, with a clean curated-biomarker demonstration of when methylation wins.

| Prior study | Finding | Relation to this work |
|---|---|---|
| Li et al., *BMC Med Inform Decis Mak* 2024 — 31-combination TCGA survival benchmark | Adding omics types most often *hurt* prediction; mRNA (±miRNA) sufficient for most cancers; methylation helped only some | Independent confirmation of our headline: one modality dominates per endpoint; curated methylation helps only on specific tasks |
| Nikolaou et al., *Cancer Res* 2023 — flexible late fusion | Dominant modality is task-specific; weight by individual CV success; pan-cancer multi-omics barely above the best single modality | Empirical basis for `select_anchor` (per-task, CV-chosen anchor) |
| Makrodimitris et al., *Brief Bioinform* 2023; Montesinos-López et al., *Front Genet* 2025 | Naive/early concatenation underperforms; PC-concatenation is hard to beat | Why we anchor and gate the residual rather than concatenate |
| Ding, Li, Narasimhan & Tibshirani, *PNAS* 2022 — cooperative learning | Degree of fusion (early↔late) chosen adaptively by CV via an agreement penalty | `anchored_gate` is a constrained, never-below-anchor special case of this |
| Tong et al., *BMC Med Inform Decis Mak* 2020; Spooner et al., *Brief Bioinform* 2024 | Complementary modalities + structured/ensemble late fusion give modest gains | The Δ > 0 regime the gate is designed to capture when it genuinely exists |
| Morandini et al., *GeroScience* 2023 (ATAC-clock); Meyer et al., *Aging Cell* 2020 | Epigenetic clocks beat transcriptomic clocks; methylation age signal is largely orthogonal to expression | Precedent for our clock > RNA result on normal-tissue age |
| Jonkman et al., *Genome Biology* 2022 | Epigenetic clocks partly track immune cell composition (naive vs activated T/NK) | Links our clock and immune-axis threads; nuance on what the clock encodes |

In short, "does multi-omics beat the best single modality?" is well-studied with mixed answers; our
results align with the skeptical-benchmark camp (esp. Li 2024, Nikolaou 2023) and add a packaged,
never-below-anchor gated integrator plus a clean curated-biomarker (Horvath-clock) case of methylation
genuinely winning.

*References:* Li et al., Does combining numerous data types improve or hinder survival prediction? A
large-scale benchmark, *BMC Med Inform Decis Mak* 2024; Nikolaou et al., Flexible late-fusion for
multi-omics survival, *Cancer Res* 2023; Makrodimitris et al., Linear vs non-linear joint embedding for
multi-omics, *Brief Bioinform* 2023; Ding, Li, Narasimhan & Tibshirani, Cooperative learning for
multiview analysis, *PNAS* 2022 (10.1073/pnas.2202113119); Hauptmann et al., Fair comparison of
multi-omics integration architectures, *BMC Bioinformatics* 2022; Montesinos-López et al., Genomic
prediction powered by multi-omics, *Front. Genet.* 2025; Tong et al., Deep-learning feature-level
integration for breast-cancer survival, *BMC Med Inform Decis Mak* 2020; Spooner et al., Benchmarking
ensemble ML for multi-omics clinical outcome prediction, *Brief Bioinform* 2024; Morandini et al.,
ATAC-clock: an aging clock from chromatin accessibility, *GeroScience* 2023; Meyer et al., BiT age: a
transcriptome-based aging clock, *Aging Cell* 2020; Jonkman et al., T and NK cell activation drives
epigenetic clock progression, *Genome Biology* 2022.
