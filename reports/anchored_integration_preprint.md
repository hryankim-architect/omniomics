# Anchored multi-omics integration and knowledge-anchored residual discovery

*A never-below-the-best-view integrator, and a tool to find what the textbook misses.*

H. Ryan Kim · omniomics-prototype · 2026

---

## Abstract

Multi-omics integration is usually sold as "more modalities → better predictions," yet on real cancer
data a strong single modality is hard to beat and naive fusion often does worse. We define an **anchored**
integrator that anchors on the empirically strongest (or a fixed prior) modality and adds the others only
as a non-negative gated residual, so it is provably never below its best single view and improves on it
only where another modality carries orthogonal signal. Across nine TCGA-BRCA endpoints the gate adds nothing
where one modality dominates and engages only on a constructed positive control — quantifying, rather than
assuming, where multi-omics helps. Generalising the anchor from "best modality" to **established knowledge**
(a textbook gene signature, a published clock) yields a zero-parameter prior that nearly matches a fully
trained model and gives the programme's only clean super-additive gain. Finally, gating genome-wide data on
a textbook prior and mining the *residual* turns the integrator into a **discovery** tool: it recovers a
verified basal/lineage axis beyond "LumA/B = proliferation" (keratinization, p≈8×10⁻¹¹) and a verified
neuroendocrine/immune axis beyond the ERBB2 amplicon, while correctly finding nothing where the textbook is
already complete (ER). The method is modality-agnostic and, above all, honest: discovery where warranted,
silence where not.

## 1. Background

The prevailing assumption that adding omics layers improves prediction is contradicted by careful
benchmarks: a large TCGA survival study found that adding data types most often *hurt* performance and that
mRNA (±miRNA) sufficed for most cancers (Li et al., 2024); late-fusion work shows the dominant modality is
task-specific and pan-cancer multi-omics barely exceeds the best single view (Nikolaou et al., 2023); and
naive concatenation routinely underperforms structured fusion (Makrodimitris et al., 2023). The honest goal
is therefore not "always fuse" but **never do worse than the best single modality, and improve on it only
where another modality carries signal the leader cannot.** This note defines an integrator built to that
contract and then repurposes it as a discovery engine.

## 2. Method

All functions live in `omniomics.multiomics`.

**Anchor selection (`select_anchor`).** Score each modality by repeated stratified-CV performance (primary),
tie-break by robustness (subtract a multiple of the CV s.d., weight by coverage). The anchor is the top
composite, chosen empirically per task inside the outer CV.

**Gated residual (`anchored_gate` / `anchored_integrate`).** Pin the anchor; add each secondary on the
anchor's residual with a non-negative weight β chosen on held-out data, β = 0 allowed. The combined score is
`logit(anchor) + β·secondary`. β = 0 ⇒ combined ≡ anchor (never below); β > 0 only where a secondary earns a
margin gain. A constrained, interpretable special case of cooperative learning (Ding et al., 2022).

**Forward gating (`forward_integrate` / `auto_integrate`).** For ≥ 3 modalities, add one at a time onto the
current residual in rank order, keeping each only if it clears the margin — forward selection over
modalities, with an interpretable record of which entered and at what β.

**Knowledge anchor (`knowledge_anchored_integrate` / `signature_score`).** The anchor need not be a data
modality: a fixed textbook signature or published clock (zero trained parameters) becomes the anchor and the
data is gated onto its residual. The gate then answers a clinically meaningful question — *does the
genome-wide data add anything beyond established biology?*

**Residual discovery (`anchored_residual_discovery`).** Anchor on the prior, gate the data, test that the
residual is signal-not-noise, and surface the features driving it that are *orthogonal* to the anchor —
candidate biology beyond the textbook — with a matched-random-panel null as the noise control.

Small-sample robustness throughout: a secondary is added only if a non-zero β beats the anchor by more than
a `gate_margin`, with each β's inner-CV score averaged over repeated splits.

## 3. Results (TCGA-BRCA)

**Anchored integration never falls below the leader.** On LumA/B, RNA alone scores 0.947 and methylation
0.745; a naive stack dips to 0.941 (below RNA) while the gated combiner stays at 0.947 (β = 0 in 9/10
repeats). A constructed positive control — a held-out methylation set RNA cannot see — flips this: the gate
engages (β = 4–8) for a significant +0.047. So the gate is both protective and capable.

**Not RNA-biased.** Run blind on nine expert-defined subtype labels, `auto_integrate` anchors on
methylation for all five methylation-defined clusters and on RNA for all four PAM50 calls — anchor
selection follows the biology, never below the leader on any of the nine.

**Knowledge anchoring — the first clean super-additive gain.** A 20-gene proliferation index with zero
trained parameters scores 0.919 on LumA/B (nearly the trained 1500-gene RNA model, 0.942); gating
genome-wide data reaches **0.947 (Δ +0.029)**, beating the pure data model. With the Horvath clock as the
anchor for normal-tissue age the clock alone (0.947) already beats RNA (0.911) and data adds ≈ 0 — the
textbook suffices and the gate says so.

**Residual discovery, verified.** Mining the proliferation-anchored residual on LumA/B surfaces the basal/
squamous-lineage axis (KRT5/14/17/6B, TP63, DSG3/DSC3, SOX10, COL17A1, KLK5/7/8; partial r ≈ 0.46). It is
verified three ways: a train-selected panel beats random panels on a held-out test in 10/10 splits (no
selection leakage), the basal core recurs in 10/10 splits, and under permuted labels the advantage collapses
(real +0.021 > all permutations). g:Profiler confirms a coherent program — cornified-envelope formation /
keratinization / epidermis development (Reactome p ≈ 8×10⁻¹¹).

**Generalises, reproduces, and discriminates.**

| endpoint | textbook anchor (0 params) | anchor AUROC | Δ data adds | verdict |
|---|---|---|---|---|
| LumA vs LumB | proliferation index (20 genes) | 0.919 | +0.029 | incomplete → basal/lineage axis (verified) |
| HER2 status | ERBB2 amplicon (9 genes) | 0.752 | +0.054 | incomplete → neuroendocrine/immune axis (verified, ≠ basal) |
| ER status | ER/luminal signature (20 genes) | 0.938 | −0.001 | complete → no hidden axis (specificity) |

The HER2 discovery is a second, *different* verified axis (held-out +0.068 in 8/8 splits; panel p = 0.024;
label-specific), a diffuse neuroendocrine/secretory + immune program distinct from basal. ER is a real-data
specificity negative: the textbook signature already matches genome-wide RNA.

**Modality-agnostic.** The discovery runs unchanged on DNA-methylation features; on LumA/B it correctly
finds no methylation axis (random and basal-targeted CpGs both Δ ≈ 0) — the basal axis is transcriptional,
not predictively methylation-encoded here.

**External validation (independent cohort).** The basal discovery reproduces in METABRIC (microarray,
n = 1175 LumA/B): the fixed TCGA-discovered panel adds Δ +0.036 over the proliferation prior (combined
0.960; beats random panels, p = 0.048), and an unbiased re-discovery in METABRIC independently recovers the
same basal genes (KRT5/14/17/6B, TP63, SOX10, DSG3/DSC3 …), overlapping the TCGA panel 20/30 (hypergeometric
p ≈ 7×10⁻²⁷). The discovered axis is a reproducible biological signal, not a single-cohort artefact. The
HER2 axis, by contrast, does *not* replicate in METABRIC — there the amplicon anchor is near-complete
(0.997 vs 0.752 in TCGA), leaving no residual (panel Δ ≈ 0, re-discovery overlap 2/30). Reproducibility
tracks biological coherence: the pathway-enriched basal axis reproduces; the diffuse HER2 axis was a
cohort-specific residual.

**Hypothesis-as-anchor.** The frame also inverts: a *hypothesis* is expressed as a candidate anchor and
gated onto the textbook anchor's residual for a three-way verdict — SUPPORTED (adds beyond the textbook),
EXPLAINED_BY_TEXTBOOK (predicts alone but redundant once the dominant prior is controlled), or REFUTED —
operationalising the observation that most signatures predict outcome and a hypothesis is novel only if it
survives adjustment for the dominant prior. On LumA/B against the proliferation anchor, a basal/keratinization
hypothesis is SUPPORTED (Δ +0.039), immune is REFUTED, and a random gene set is EXPLAINED_BY_TEXTBOOK
(predicts at 0.63 alone but adds ≈ 0). Scaled across the 50 MSigDB Hallmark sets, the screen is
self-validating — proliferation hallmarks (E2F/G2M/MYC) add exactly 0 beyond the anchor — and, requiring
agreement across two proliferation anchors plus FDR, the robust SUPPORTED hits are the estrogen-response
lineage programs (matching LumA/B biology). This pathway-level ranking is cohort/platform-sensitive (it does
not reproduce in METABRIC microarray, where the anchor-family agreement filter correctly nulls
proliferation hallmarks instead), so it is a within-cohort hypothesis-ranking tool, not a cross-platform
claim — unlike the gene-level basal discovery, which reproduces. To keep this honest, `hypothesis_anchor_test`
returns a gate-free commonality/mediation decomposition (`collinearity_label` ∈ {NOVEL, REDUNDANT, INERT}) that
separates an *absent* hypothesis from a merely *redundant* (collinear) one: the estrogen-response axis is NOVEL
in TCGA (unique R² 0.038) but REDUNDANT in METABRIC (redundancy 1.00, 96 % mediated through proliferation). A
controlled transportability sweep — fixing both marginal effects and varying only the anchor–hypothesis
correlation — shows the *same* ER effect is 100 % NOVEL at TCGA's correlation (+0.19) yet collapses into the
collinear/suppression valley at METABRIC's (−0.17): the verdict is a covariate-distribution (transportability)
property, not a change in biology.

## 4. Related work

The approach instantiates three established traditions: clinical-offset / incremental-value modelling, where
an established model is an offset and omics are selected on the residual (Volkmann et al., 2019, a near-exact
methodological twin); biologically-informed models that bake known biology into structure (Elmarakeby et
al., 2021, P-NET; Liu et al., 2024, Pathformer); and ML on established signatures (Sarafidis et al., 2024).
Our contribution is the packaged, **never-below-anchor, margin-gated** residual on a zero-parameter prior,
plus the residual-discovery framing with a matched-random-panel noise control.

## 5. Limitations

Predictive gains over a strong anchor are modest; the method's value is *routing and discovery*, not large
AUROC wins. Genuine super-additive multi-omics gains are rare on real endpoints (one constructed positive
control; one real knowledge-anchored gain). The basal discovery is externally validated in METABRIC
(replication + independent re-discovery); the HER2 axis tested negative there (cohort-specific — the amplicon
anchor is near-complete in METABRIC), so external reproducibility tracks biological coherence. Discovered
axes are candidate hypotheses, not causal claims.

## 6. Conclusion

Anchor on the leading view (or the textbook), add the rest as a gated residual, and let the gate tell you —
honestly — whether the genome-wide data beats established biology. Where it does, the same residual names the
new axis (basal for LumA/B; neuroendocrine/immune for HER2); where it does not (ER), the method stays silent.
Discovery where warranted, silence where not.

## References

Li et al., *BMC Med Inform Decis Mak* 2024 · Nikolaou et al., *Cancer Res* 2023 · Makrodimitris et al., *Brief
Bioinform* 2023 · Ding, Li, Narasimhan & Tibshirani, *PNAS* 2022 · Volkmann et al., *BMC Med Res Methodol*
2019 · Elmarakeby et al., *Nature* 2021 · Liu et al., *Bioinformatics* 2024 · Sarafidis et al., *ESMO Open*
2024 · Morandini et al., *GeroScience* 2023 · Horvath, *Genome Biology* 2013.

*Code, runners, recorded metrics and CI guards: the `omniomics-prototype` repository. Full method note:
`reports/anchored_integration_methods.md`.*
