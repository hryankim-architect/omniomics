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

*Biological characterization.* g:Profiler over-representation confirms the basal panel is a coherent
program — **cornified-envelope formation / keratinization / epidermis development** (Reactome p ≈ 8×10⁻¹¹;
9/30 genes in the cornified-envelope set) — i.e. a textbook basal/squamous-lineage signature, not a
statistical fluke. The HER2 panel, by contrast, is verified-predictive but shows **no single enriched
pathway** (a diffuse, multi-program set — neuroendocrine/secretory and immune by individual-gene function),
an honest reminder that a verified predictive axis need not map to one annotated pathway. `discovery_enrichment.csv`.

*Generalises, reproduces, and discriminates (three anchors/endpoints).* The same machinery was applied to
two further textbook anchors. **HER2 status**, anchored on the ERBB2 17q12 amplicon (ERBB2, GRB7, STARD3, …;
9 genes, 0 trained parameters, n = 808): the amplicon is *incomplete* (AUROC 0.752), so the method finds a
real, verified, anchor-orthogonal axis (Δ +0.054; discovered panel +0.085 vs random +0.023, p = 0.024;
held-out positive in 8/8 splits; label-specific) — a **neuroendocrine/secretory + immune** program (DOC2A,
CAMK2B, NPW, RAMP1, BEX1, PCSK1N, SULT4A1, MS4A8B, AZU1, SEMA3D) **distinct from the LumA/B basal axis**.
**ER status**, anchored on a textbook ER/luminal signature (ESR1, GATA3, FOXA1, …; 20 genes, n = 1152): the
signature is *complete* (0.938) and genome-wide RNA adds ≈ 0 (Δ −0.001) — the method correctly returns **no
hidden axis**.

| endpoint | textbook anchor (0 params) | anchor AUROC | data adds (Δ) | verdict |
|---|---|---|---|---|
| LumA vs LumB | proliferation index (20 genes) | 0.919 | **+0.029** | incomplete → **basal / lineage axis** (verified) |
| HER2 status | ERBB2 amplicon (9 genes) | 0.752 | **+0.054** | incomplete → **neuroendocrine/secretory + immune axis** (verified, ≠ basal) |
| ER status | ER/luminal signature (20 genes) | 0.938 | −0.001 | complete → **no hidden axis** (specificity) |

So the method is **not a one-off**: it makes a real, verified discovery a *second* time and finds a
*different*, context-appropriate axis (HER2 → neuroendocrine/immune; LumA/B → basal), while correctly finding
**nothing** where the textbook is already complete (ER) — discovery where warranted, silence where not.
Runners `reports/dmoi_discovery_er.py`, `reports/dmoi_discovery_her2.py`; `discovery_{er,her2}_results.csv`,
`discovery_her2_verification.csv`.

*Modality-agnostic.* The discovery operates on any feature matrix, not just RNA. Run with the DATA modality
set to DNA methylation (450K CpGs) and the same RNA proliferation prior, on LumA/B, it runs unchanged but
correctly finds **no methylation axis**: random genome-wide CpGs add ≈ 0 (Δ −0.001) and even basal-gene-
targeted CpGs add ≈ 0 (Δ +0.000; basal methylation alone only 0.756). So the basal/lineage axis is a
*transcriptional* signal, not predictively methylation-encoded for this endpoint — and the method does not
manufacture a methylation discovery. `reports/dmoi_discovery_methylation.py`; `methylation_discovery_results.csv`.

*External validation (independent cohort).* The basal discovery reproduces in **METABRIC** (microarray,
n = 1175 LumA/B). The TCGA-discovered basal panel — fixed, so no selection leakage in METABRIC — gated onto
the same proliferation prior adds **Δ +0.036 (combined 0.960)** and beats matched random panels (0/20 beat
it, p = 0.048); and running the discovery *unbiased* on METABRIC independently recovers the basal axis (top
genes KRT5/14/17/6B/16, COL17A1, TP63, SOX10, DSG3/DSC3, KLK5/6/7), overlapping the TCGA panel **20/30**
(hypergeometric p ≈ 7×10⁻²⁷). So the discovered axis is a reproducible biological signal, not a TCGA
artefact. `reports/dmoi_external_metabric.py`; `external_validation_metabric.csv`.

The **HER2 axis, in contrast, does *not* replicate** in METABRIC — for an interpretable reason. There the
ERBB2 amplicon anchor is near-complete (AUROC 0.997 vs 0.752 in TCGA: METABRIC's HER2 calls are far more
amplification-concordant), so essentially no residual exists; the TCGA HER2 panel adds Δ +0.000 (p = 1.0)
and an unbiased re-discovery overlaps the TCGA panel only 2/30 (vs 20/30 for basal). So the discovery's
reproducibility **tracks its biological coherence**: the clean, pathway-enriched basal/keratinization axis
reproduces across cohorts, while the diffuse HER2 axis was a *cohort-specific* residual that vanishes where
the textbook anchor happens to be complete. The method stays well-behaved throughout — it finds nothing
where there is no residual (METABRIC HER2, like TCGA ER). `reports/dmoi_external_her2_metabric.py`;
`external_validation_her2_metabric.csv`.

*Cross-domain (a completely different dataset and question).* The method is neither breast-cancer- nor
expression-specific. On **NSCLC patients receiving anti-PD-1 checkpoint blockade** (Hellmann/MSK 2018,
n = 227; mutation/clinical features, endpoint = durable clinical benefit), anchoring on the textbook
immuno-oncology biomarker **tumour mutational burden (TMB)** and mining the residual independently recovers
the field's *other* established biomarkers: **PD-L1 score** (positive, and genuinely orthogonal to TMB:
corr 0.00), **EGFR mutation** and **STK11/LKB1 mutation** (both negative — known checkpoint-blockade
resistance). Discovered panel Δ +0.061 vs random +0.006 (p = 0.038). So anchored on the textbook biomarker,
the method recovers the known complementary biomarkers in a new domain. `reports/dmoi_discovery_nsclc_io.py`;
`discovery_nsclc_io_results.csv`.

*Cross-cancer replication of a discovered axis.* The strongest test that a discovered axis is real biology
rather than a cohort artefact is to look for it in a *different cancer*. The breast basal/keratinization
residual should re-appear wherever the squamous lineage does. On **TCGA lung carcinoma** (LUAD adeno vs LUSC
squamous, n = 1,129; UCSC Xena), anchoring on the *same* zero-parameter proliferation signature (AUROC 0.77
for histology — incomplete) and mining the residual recovers the squamous/keratinization axis, and it
overlaps the breast basal panel **10/30** (hypergeometric p ≈ 3×10⁻¹⁶): the same genes (KRT5, KRT14, KRT6B,
TP63, DSG3, DSC3, FAT2, CALML3, ANXA8, TRIM29) discovered in a different cancer. Honest caveat: here the
panel-vs-random delta control *saturates* (p ≈ 0.18) because squamous-vs-adeno is a near-trivial
transcriptomic split — random high-variance panels also separate it — so on this easy endpoint the
informative metric is the gene-level replication, which is decisive. `reports/dmoi_external_lung.py`;
`external_validation_lung.csv`.

*A complement for the saturating null — selection stability.* The panel-vs-random null saturates whenever an
endpoint is broadly predictable (lung histology above; ER status in the genome-wide methylation run), because
many panels predict it. `anchored_residual_discovery(stability_reps=N)` adds a robust, difficulty-independent
signal: the mean recurrence of the discovered panel across 50% subsamples, measured against a permuted-label
null (subsamples overlap, so raw recurrence is only meaningful relative to label-shuffled data). It returns
`stability`, `stability_null`, and `stability_gain = stability − stability_null`. On the lung basal axis,
where panel-vs-random is uninformative (p ≈ 0.22), stability is decisive: the keratinization panel re-selects
in ~90% of subsamples versus ~1% under permuted labels (**stability_gain ≈ +0.89**), confirming a real,
reproducible structured axis. The metric ranks real axes far above spurious in-sample noise (≈ +0.7 vs ≈ +0.25
on synthetic controls), complementing — not replacing — the FDR count and biological coherence.

*Anchor standardization — robustness to the choice of textbook anchor.* The discovery's honesty depends on
the anchor's provenance, which was hand-picked. Following Venet et al. (2011, *PLoS Comput Biol*) — most
random signatures predict outcome, and adjusting for a proliferation meta-gene removes the effect — we (a)
added a reproducible, data-driven anchor recipe (`marker_correlated_anchor`, the meta-PCNA construction: top
genes correlated with a canonical marker), and (b) ran the LumA/B discovery across a *family* of
biologically-equivalent proliferation priors (our 20-gene curated index, MSigDB `HALLMARK_E2F_TARGETS` /
`G2M_CHECKPOINT` / `MYC_TARGETS_V1`, and meta-PCNA) as a vibration-of-effects check. The basal/keratinization
axis is robust to anchor choice: **every** anchor re-recovers the basal core (KRT5/14/17/6B, TP63, DSG3, DSC3)
7/7 and 21–26/30 of the reference panel; mean pairwise Jaccard of the discovered panels **0.69**; 24 genes
are consensus across ≥4/5 anchors (21 in the basal panel). So the discovery is not an artifact of a single
hand-picked anchor. `reports/dmoi_anchor_family_voe.py`; `anchor_family_voe.csv`, `anchor_family_consensus.csv`;
full rationale in `reports/anchor_standardization_discussion.md`.

*Hypothesis-as-anchor — confirm / explain-away / refute a hypothesis with the textbook anchor + real data.*
The same frame turns a hypothesis into a *candidate anchor* and tests it against the established textbook
anchor: `hypothesis_anchor_test` gates the hypothesis onto the textbook anchor's residual and returns a
3-way verdict — **SUPPORTED** (adds signal beyond the textbook → candidate novel mechanism),
**EXPLAINED_BY_TEXTBOOK** (predicts on its own but is redundant once the textbook prior is controlled), or
**REFUTED** (neither). This operationalizes Venet et al. (2011): a signature is a novel mechanism only if it
survives adjustment for the dominant prior. On LumA-vs-LumB with the proliferation textbook anchor, the
**basal/keratinization** hypothesis is SUPPORTED (Δ +0.039 beyond the anchor, β = 8; AUROC 0.78, orthogonal
to proliferation), **immune/cytotoxic** is REFUTED (AUROC 0.51 ≈ chance, adds 0.000 — immune does not
separate LumA/B), and a **random 30-gene** set is EXPLAINED_BY_TEXTBOOK (predicts alone at 0.63 but adds
≈ 0 once proliferation is controlled — exactly the Venet artifact). `reports/dmoi_hypothesis_anchor.py`;
`hypothesis_anchor_results.csv`.

*Hypothesis screen — a whole standard library at once.* `rank_hypotheses` batches the test over a library:
we scored all 50 MSigDB Hallmark gene sets as candidate hypothesis anchors and ranked them by signal added
beyond the proliferation textbook anchor on LumA/B. The screen is internally self-validating: the
proliferation-type hallmarks (`E2F_TARGETS` 0.914, `G2M_CHECKPOINT` 0.894, `MYC_TARGETS_V1/V2`) predict
strongly *alone* but add **exactly 0.0** beyond the anchor → correctly EXPLAINED_BY_TEXTBOOK (they *are* the
anchor's axis). The SUPPORTED hits (11/50) are the genuinely orthogonal lineage programs — led by
**ESTROGEN_RESPONSE (early & late)**, plus P53, WNT, NOTCH, apical-surface/junction, myogenesis — matching
known biology (LumA vs LumB differ in proliferation *and* ER signalling). `reports/dmoi_hypothesis_screen.py`;
`hypothesis_screen_results.csv`.

A **robust** variant (`reports/dmoi_hypothesis_screen_robust.py`) adds two safeguards: anchor-family
averaging (require support beyond *both* the curated and meta-PCNA proliferation anchors) and Benjamini-
Hochberg FDR. Requiring agreement across anchors tightens the set from 11 to **4** robust hits —
ESTROGEN_RESPONSE (late & early), APICAL_SURFACE, P53 — i.e. the ER lineage axis is the anchor-robust biology
for LumA/B, while WNT/NOTCH/myogenesis were anchor-dependent and drop out. (The empirical-null FDR is liberal
here because the null is very tight, so anchor-family agreement is the binding robustness gate.)
`hypothesis_screen_robust.csv`.

Repeating the screen on **other endpoints** is a specificity check (`reports/dmoi_hypothesis_screen_endpoints.py`).
For HER2 (anchored on the ERBB2 17q12 amplicon, AUROC 0.92 here) and ER (anchored on a textbook ER/luminal
signature, 0.94) the textbook anchor is already *complete*, and the screen correctly surfaces **0/50**
SUPPORTED hallmarks — it does not manufacture hypotheses beyond a saturated anchor, the complement of the
incomplete-anchor LumA/B case where the ER lineage axis is SUPPORTED. `hypothesis_screen_her2.csv`,
`hypothesis_screen_er.csv`.

*External reproduction in METABRIC — an honest negative with a methodological lesson*
(`reports/dmoi_hypothesis_screen_metabric.py`; `hypothesis_screen_metabric.csv`). Repeating the robust LumA/B
screen on METABRIC (microarray, n = 1,175) does **not** reproduce the TCGA pathway-level result: estrogen-
response is REFUTED there (AUROC 0.58, adds ≈ 0) and the SUPPORTED sets do not overlap. Two things are
instructive. (i) On microarray the 20-gene curated proliferation anchor is a weaker proxy, so E2F/G2M/MYC
hallmarks show *large* deltas beyond it (0.08–0.12) — but the **anchor-family agreement filter correctly nulls
them** (`both_anchors_support = False`: they add ≈ 0 beyond meta-PCNA, which fully captures proliferation), so
no proliferation pathway is falsely called SUPPORTED. (ii) The *pathway-level* "what adds beyond the textbook"
question is therefore more cohort/platform-sensitive than the *gene-level* discovery — the basal/keratinization
**genes** reproduced robustly across cohorts and cancers (METABRIC 20/30, lung, HNSC), but the Hallmark
hypothesis ranking does not transfer cleanly across RNA-seq and microarray. The honest reading: trust the
anchor-family-robust, gene-level discoveries; treat the pathway screen as a within-cohort hypothesis-ranking
tool, not a cross-platform claim.

*Root cause of the non-reproduction* (`reports/dmoi_hypothesis_metabric_diagnosis.py`;
`hypothesis_metabric_diagnosis.csv`). The screen measures a *conditional* (partial) association — "does the
hypothesis add **beyond** the anchor?" — which depends on the joint correlation structure of anchor and
hypothesis. Estrogen-response actually separates LumA from LumB in *both* cohorts (Cohen d ≈ −0.2). The
difference is the proliferation–ER correlation: **+0.19 in TCGA vs −0.17 in METABRIC (a sign flip)**. In
METABRIC, LumB's lower ER is collinear with its higher proliferation, so adjusting for proliferation removes
the ER signal (partial corr ER∣prolif = **−0.007** ≈ 0 → REFUTED); in TCGA the two are not collinear, so ER
keeps orthogonal information (partial corr **−0.243** → SUPPORTED). So ER biology is present in both — it is
*redundant with proliferation* in METABRIC, not absent. This is an intrinsic property of conditional
inference (the verdict is anchor-relative), and a precise reason the pathway screen is best read
within-cohort.

To make this honest rather than misleading, `hypothesis_anchor_test` now also returns a **commonality
decomposition** (unique vs shared variance; Tonidandel & LeBreton 2011) and a **mediation split** (direct vs
indirect through the anchor), with a `collinearity_label` ∈ {NOVEL, REDUNDANT, INERT}. This distinguishes a
hypothesis that is *absent* from one that is merely *redundant* (collinear) with the anchor. Applied to the
ER case: in **TCGA** ER is **NOVEL** (unique R² = 0.038); in **METABRIC** it is **REDUNDANT** — its signal is
not absent but **shared with proliferation (redundancy = 1.00, 96 % mediated through the proliferation
anchor)**. So the cross-cohort difference is correctly characterised as a change in collinearity structure,
not a disappearance of ER biology.

*Transportability of the verdict (a controlled sweep).* The framing above implies the verdict is a
**transportability** quantity (Degtiar & Rose 2023): it depends on each cohort's covariate distribution
(here, corr(prolif, ER)), not on the hypothesis's marginal effect. `transportability_sweep` makes this
quantitative and fully reproducible without any external data: holding both marginal separations fixed at the
BRCA values (Cohen's d_prolif ≈ 1.85, d_ER ≈ −0.22) and varying *only* the residual anchor–hypothesis
correlation, it simulates many cohorts per correlation and records the fraction labelled NOVEL vs REDUNDANT.
The two real cohorts then sit on the curve at their measured correlations: at **TCGA's +0.19 the ER hypothesis
is 100 % NOVEL**, while at **METABRIC's −0.17 it falls into the collinear/suppression valley (≈ 2 % NOVEL,
≈ 42 % REDUNDANT)** — the *same* ER effect, opposite verdicts, driven only by the nuisance correlation. This is
the clearest statement of why an anchored hypothesis screen must be read within-cohort and re-characterised
(not just re-run) when ported across cohorts.

*A second endpoint (HER2) — specificity of the labels.* Repeating the analysis on HER2 status (anchor = the
ERBB2 17q12 amplicon, hypothesis = the ER signature) gives a *different* regime, which the labels capture
correctly. Here the ER→HER2 marginal effect itself differs across cohorts, so HER2 is not a clean
same-marginal flip: in **TCGA** ER barely separates HER2± (Cohen's d ≈ −0.05) and is labelled **INERT** (no
appreciable signal to be novel or redundant), whereas in **METABRIC** ER is a moderate signal (d ≈ −0.29) that
is collinear with the amplicon and is labelled **REDUNDANT** (redundancy 0.89, 66 % mediated). The
transportability sweep at METABRIC's own marginals reproduces the REDUNDANT valley at the observed
correlation. So across two endpoints the framework distinguishes all three failure modes — *novel* (ER in
TCGA-LumA/B), *redundant/collinear* (ER in METABRIC-LumA/B and HER2), and *absent/weak* (ER for HER2 in TCGA)
— rather than collapsing them into a single "fails to add."

*A multi-endpoint × cohort panel.* Scaling the characterisation to four breast-cancer endpoints across four
columns — TCGA RNA-seq, **TCGA Agilent microarray (the same patients on a different platform)**, METABRIC
(an independent microarray cohort), and **SCAN-B (GSE96058, a fully independent Swedish RNA-seq cohort of
~3,400 tumours)** — yields a compact map of where the label transports. Two endpoints transport, two do not.
The transportable ones are biologically robust: **Basal-vs-rest** with a basal/keratinization anchor and an
immune hypothesis is **NOVEL in all four columns** (immune infiltration adds a real axis beyond the basal
lineage program), and **ER-status** with the ER-signature anchor and a proliferation hypothesis carries a small
but real unique proliferation slice (NOVEL) in all four — robust across platform *and* independent cohort. The
two that do *not* transport are exactly the ER-collinearity cases studied above — **LumA-vs-LumB**
(proliferation→ER) and **HER2** (amplicon→ER). The four-column view exposes the mechanism cleanly: for
LumA-vs-LumB the label tracks **measurement technology**, NOVEL on *both* RNA-seq cohorts (TCGA RNA-seq and the
independent SCAN-B) yet REDUNDANT on *both* microarrays (TCGA Agilent and METABRIC) — it even flips on the same
TCGA patients between their RNA-seq and Agilent measurements, before any change of cohort. So the
transportability caveat is specific and predictable: it attaches to hypotheses whose collinearity with the
anchor is itself measurement-dependent, while genuinely anchor-orthogonal axes (basal→immune) transport cleanly
across platform and cohort. Runner: `reports/dmoi_endpoint_panel.py` (`endpoint_panel.csv`,
`figs/endpoint_panel.png`); the Agilent column auto-adds when AgilentG4502A_07_3.gz (UCSC Xena) is present and
the SCAN-B column when the pre-extracted GSE96058 marker matrix + phenotype are in SCANB_DIR.

*Tissue-independence (a third cancer, head & neck).* HNSC is uniformly squamous, so it offers no
within-cohort histology contrast — but that allows a confound control. Scoring TCGA HNSC (head & neck
squamous), LUSC (lung squamous) and LUAD (lung adeno) with the breast 30-gene basal panel, the score
separates squamous from adeno at AUROC 0.96, and HNSC (head & neck) and LUSC (lung) — two *different
tissues* — both score high while LUAD scores low (median 0.63, 0.18 vs −0.69). So the discovered axis tracks
squamous *lineage*, not tissue of origin. Within HNSC it also tracks the textbook biology: well-differentiated
(G1, keratinizing) tumours score higher than poorly differentiated (G3) (one-sided p ≈ 2×10⁻⁴).
`reports/dmoi_external_hnsc.py`; `external_validation_hnsc.csv`.

*Clinical significance — an honest separation of identity from outcome.* A reproducible axis need not be
prognostic. In TCGA-BRCA (n = 866, 132 events) the basal/keratinization score is **not** associated with
overall survival (univariate Cox HR 1.01, p = 0.89; adjusted for the proliferation anchor p = 0.36; KM
median-split log-rank p = 0.29), whereas the proliferation score is (p = 0.001) — so the discovered axis
captures lineage *identity*, not outcome. It is clinically coherent on identity: the score marks ER-negative
/ basal-like disease (AUROC 0.70). Reporting this negative honestly matters — the method finds real, robust
biology, and that biology happens to be a differentiation marker rather than a survival driver.
`reports/dmoi_clinical_survival.py`; `clinical_basal_survival.csv`.

*Genome-wide scale (end-to-end on real data).* The scale-prep stack was exercised on a real
**485,577-probe × 829-sample** TCGA-BRCA HumanMethylation450 matrix that does not comfortably fit in memory.
A one-pass streaming SIS screened all 485k probes to the top 5,000 by association with ER status in ≈ 24 s
(low memory), and the vectorized partial-correlation discovery with BH-FDR and parallel nulls then ran on
the survivors, anchored on the textbook ESR1-methylation prior (ER− tumours show ESR1 promoter
hypermethylation). The anchor predicts ER status (AUROC 0.65, incomplete); genome-wide CpGs lift it to 0.90;
and the top CpG discovered *beyond* ESR1 maps to **PGR** (progesterone receptor — the canonical
ER-coregulated gene), with 4,837/5,000 screened CpGs FDR-significant (ER status has a broad methylation
footprint). As with the easy histology split, the panel-vs-random null saturates, so the informative
outputs are the anchor AUROC, the FDR count, and the biologically coherent top hits. This confirms the
pipeline runs end-to-end at genome-wide scale and recovers textbook biology. `reports/dmoi_meth_genomewide.py`;
`meth_genomewide_results.csv`, `meth_genomewide_novel.csv`.

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
| Ellen, Nikolaou et al., *Sci Rep* 2023 — autoencoder multimodal NSCLC survival | Best multimodal model rarely beats a small two-modality subset; useful combination is task-specific and narrow | Empirical basis for `select_anchor` (per-task, CV-chosen anchor) |
| Makrodimitris et al., *Brief Bioinform* 2023; Montesinos-López et al., *Front Genet* 2025 | Naive/early concatenation underperforms; PC-concatenation is hard to beat | Why we anchor and gate the residual rather than concatenate |
| Ding, Li, Narasimhan & Tibshirani, *PNAS* 2022 — cooperative learning | Degree of fusion (early↔late) chosen adaptively by CV via an agreement penalty | `anchored_gate` is a constrained, never-below-anchor special case of this |
| Tong et al., *BMC Med Inform Decis Mak* 2020; Spooner et al., *Brief Bioinform* 2024 | Complementary modalities + structured/ensemble late fusion give modest gains | The Δ > 0 regime the gate is designed to capture when it genuinely exists |
| Morandini et al., *GeroScience* 2023 (ATAC-clock); Meyer et al., *Aging Cell* 2020 | Epigenetic clocks beat transcriptomic clocks; methylation age signal is largely orthogonal to expression | Precedent for our clock > RNA result on normal-tissue age |
| Jonkman et al., *Genome Biology* 2022 | Epigenetic clocks partly track immune cell composition (naive vs activated T/NK) | Links our clock and immune-axis threads; nuance on what the clock encodes |

In short, "does multi-omics beat the best single modality?" is well-studied with mixed answers; our
results align with the skeptical-benchmark camp (esp. Li 2024, Nikolaou 2023) and add a packaged,
never-below-anchor gated integrator plus a clean curated-biomarker (Horvath-clock) case of methylation
genuinely winning.

> **Citation note.** The canonical, PubMed-verified, DOI-linked reference list for this work is in
> `reports/anchored_integration_manuscript.md` (§References). Two earlier entries were corrected on
> re-verification: the late-fusion NSCLC-survival paper is **Ellen, Jacob, Nikolaou & Markuzon, *Sci Rep*
> 2023** (10.1038/s41598-023-42365-x) — not *Cancer Res* — and the joint-embedding comparison is
> **Makrodimitris et al., *Brief Bioinform* 2023** (10.1093/bib/bbad416). A third earlier entry
> ("Sarafidis, *ESMO Open* 2024") was a mis-citation and is replaced in the manuscript by Boscolo Bielo
> et al., *ESMO Open* 2026 (10.1016/j.esmoop.2026.106064).

*References:* Li et al., Does combining numerous data types improve or hinder survival prediction? A
large-scale benchmark, *BMC Med Inform Decis Mak* 2024 (10.1186/s12911-024-02642-9); Ellen, Jacob, Nikolaou
& Markuzon, Autoencoder-based multimodal prediction of non-small cell lung cancer survival, *Sci Rep* 2023
(10.1038/s41598-023-42365-x); Makrodimitris et al., An in-depth comparison of linear and non-linear joint
embedding methods for bulk and single-cell multi-omics, *Brief Bioinform* 2023 (10.1093/bib/bbad416); Ding,
Li, Narasimhan & Tibshirani, Cooperative learning for multiview analysis, *PNAS* 2022 (10.1073/pnas.2202113119); Hauptmann et al., Fair comparison of
multi-omics integration architectures, *BMC Bioinformatics* 2022; Montesinos-López et al., Genomic
prediction powered by multi-omics, *Front. Genet.* 2025; Tong et al., Deep-learning feature-level
integration for breast-cancer survival, *BMC Med Inform Decis Mak* 2020; Spooner et al., Benchmarking
ensemble ML for multi-omics clinical outcome prediction, *Brief Bioinform* 2024; Morandini et al.,
ATAC-clock: an aging clock from chromatin accessibility, *GeroScience* 2023; Meyer et al., BiT age: a
transcriptome-based aging clock, *Aging Cell* 2020; Jonkman et al., T and NK cell activation drives
epigenetic clock progression, *Genome Biology* 2022.
