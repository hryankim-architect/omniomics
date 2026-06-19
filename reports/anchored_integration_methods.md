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

*References:* Ding, Li, Narasimhan & Tibshirani, Cooperative learning for multiview analysis, *PNAS*
2022 (10.1073/pnas.2202113119); Nikolaou et al., Flexible late-fusion for multi-omics survival,
*Cancer Res* 2023; Hauptmann et al., Fair comparison of multi-omics integration architectures, *BMC
Bioinformatics* 2022; Montesinos-López et al., Genomic prediction powered by multi-omics, *Front.
Genet.* 2025.
