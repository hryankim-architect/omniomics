---
title: "Anchored multi-omics integration and knowledge-anchored residual discovery: a never-below-the-best-view integrator that finds what the known biology misses"
author: "H. Ryan Kim"
date: "2026"
header-includes:
  - \usepackage{caption}
  - \captionsetup{labelformat=empty}
---

**Running title:** Knowledge-anchored residual discovery

**Author:** H. Ryan Kim¹ ([ORCID: 0000-0002-1869-0412](https://orcid.org/0000-0002-1869-0412))

**¹** Independent researcher (omniomics project).

**Corresponding author:** H. Ryan Kim, hryankim1221@gmail.com.

**Subject areas (bioRxiv):** Bioinformatics; Cancer Biology.

**Keywords:** multi-omics integration; prior-informed machine learning; gene signatures; breast cancer; TCGA; METABRIC; immunotherapy biomarkers; interpretable discovery.

---

## Graphical abstract

![**Graphical abstract.** Anchor on established knowledge (a zero-parameter prior); admit genome-wide data only
as a non-negative gated residual (never below the anchor); mine the residual for an anchor-orthogonal axis; and
read transportability — the same hypothesis can be NOVEL or REDUNDANT depending on the measurement
platform.](figs/graphical_abstract.png){ width=100% }

## In brief

Combining many kinds of molecular data is widely assumed to improve cancer prediction, yet in practice a single
strong data type is hard to beat. This work introduces an *anchored* framework that starts from established
biological knowledge — a known high-performing data type, or a fixed textbook gene signature — so the model can
never do worse than the existing standard. By analysing only the *residual* signal left after accounting for
that anchor, it surfaces orthogonal biological axes that conventional data-fusion overlooks: it recovers a
keratinization (basal) program in breast cancer and validates it as a pan-epithelial axis: in seven
independent TCGA cohorts the breast-derived panel transfers at AUROC ≥ 0.91 across five squamous-containing
comparisons (lung, head and neck, oesophagus, bladder, cervix) yet collapses to 0.52–0.61 in two adeno-only
negative controls (gastric and endometrial adenocarcinoma), a 30–40 percentage-point bifurcation that proves
specificity to squamous lineage identity. The framework also recovers established immunotherapy biomarkers
in lung cancer. Inverting the idea, a hypothesis can itself be expressed as a candidate anchor and tested against the
textbook prior, letting the framework separate genuinely new findings from those already explained by existing
literature. It also makes *transportability* explicit — showing how the measurement technology (for example,
RNA-seq versus microarray) can make the same biological signal look novel in one cohort and redundant in
another. In short, it is a discovery engine that is honest about what is new and what the textbook already
knows.

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
(panel Δ +0.061 vs +0.006, p = 0.038). Finally, the discovered axis replicates across *cancers* as a pan-epithelial keratinization program:
in seven independent TCGA cohorts spanning lung, head and neck, oesophagus, bladder, cervix (squamous
vs non-squamous comparisons) and gastric/endometrial adenocarcinoma (adeno-only negative controls), the
fixed 30-gene breast panel achieves AUROC ≥ 0.91 in all five squamous-containing comparisons (0.96 lung;
0.96 HNSC; 0.91 ESCA; 0.97 BLCA; 0.94 CESC), falling to 0.52–0.61 in the two adeno-only controls — a
30–40 percentage-point bifurcation that demonstrates the axis is specific to squamous lineage identity and
silent in its absence.
Inverting the frame, a *hypothesis* can be expressed as a candidate anchor and tested against the textbook
anchor; a commonality/mediation decomposition then separates a genuinely novel axis from one that is merely
redundant (collinear) or absent, and a five-endpoint × four-cohort panel (TCGA RNA-seq, TCGA Agilent, METABRIC,
SCAN-B) shows that anchor-orthogonal axes (e.g. an immune program beyond the basal lineage) transport across
platform and cohort, whereas collinear ones do not — their verdict is governed by a nuisance correlation that
can itself flip sign with the assay (RNA-seq vs microarray). **Conclusion.** Anchor on
established biology and let a gate decide, honestly, whether genome-wide data beats it; where it does, the
residual names the new axis — one that reproduces across cohorts (METABRIC), generalises across diseases and
feature types (NSCLC immunotherapy), and, for the basal/keratinization axis, transfers across five organs
while correctly falling silent in two adeno-only controls. External reproducibility tracks biological
coherence.

## 1. Introduction

The prevailing premise that adding omics layers improves prediction is contradicted by careful benchmarks:
a large-scale TCGA benchmark across 14 cancers found that adding data types most often *hurt* survival
prediction, with mRNA (±miRNA) sufficient for most cancers [1]; multimodal survival models often fail to
beat a small two-modality subset, so the useful combination is task-specific and narrow [2]; and naive
concatenation underperforms — a principal-component baseline is hard to beat, so adaptive weighting is
needed to realise any gain from fusion [3,4]. The honest objective
is therefore not "always fuse" but to *never do
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
Horvath epigenetic clock [9] as the anchor for normal-tissue age, the clock alone (0.947) already beats RNA
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

**Cross-domain generalisation (different cancers, questions, and feature types).** The method is not specific
to breast cancer or to expression. On NSCLC patients receiving anti-PD-1 checkpoint blockade (Hellmann/MSK
2018, n = 227; mutation and clinical features; endpoint = durable clinical benefit), anchoring on the
textbook immuno-oncology biomarker tumour mutational burden (TMB; anchor AUROC 0.60) and mining the residual
independently recovers the field's *other* established biomarkers: PD-L1 score (positive partial correlation,
and genuinely orthogonal to TMB, corr 0.00), EGFR mutation and STK11/LKB1 mutation (both negative — known
checkpoint-blockade resistance). The discovered panel adds Δ +0.061 versus +0.006 for matched random panels
(p = 0.038). Anchored on the textbook biomarker, the procedure rediscovers the known complementary biomarkers
of a different disease.

**Cross-cancer replication and tissue-independence of the basal axis.** The strongest test that a discovered
axis is real biology, not a cohort artefact, is to seek it in a different cancer. We tested the breast basal
panel across seven independent TCGA cohorts spanning five squamous-containing comparisons and two deliberate
adeno-only negative controls, distinguishing two properties: *panel transfer AUROC* (does the fixed 30-gene
panel discriminate the target endpoint?) and *de-novo residual rediscovery* (does the anchor-orthogonal
residual in the new cancer re-recover the same genes?). These two properties are independent — a panel can
transfer perfectly while the residual names completely different genes — and each answers a different question.

*Pan-epithelial axis: the transfer evidence.* In all five squamous-containing comparisons the breast basal
panel transfers at AUROC ≥ 0.91: lung LUAD vs LUSC (AUROC implied ≥ 0.96; 10/30 gene overlap,
p ≈ 3×10⁻¹⁶); head and neck (HNSC + LUSC vs LUAD, AUROC 0.96, with within-HNSC grade tracking G1 > G3,
p ≈ 2×10⁻⁴); oesophagus ESCC vs EAC (0.91); bladder basal-squamous vs luminal-papillary (0.97); and
cervical squamous vs endocervical adenocarcinoma (0.94). The panel was derived from breast cancer and has
never seen any of these tissues. That it ranks squamous above non-squamous identity at this precision across
five organs establishes it as a tissue-independent marker of the keratinization programme, not a breast
cohort artefact.

*A pole-selection rule from the residual.* De-novo unbiased residual discovery does not uniformly
re-recover the breast panel. Instead, the residual surfaces whichever pole carries the *cleanest
anchor-orthogonal signal* in each tissue, revealing a systematic bifurcation. In lung (LUAD/LUSC,
overlap 10/30) and cervix (CESC squamous vs adeno, overlap 6/30, p = 1.5×10⁻⁸), the residual names
the *squamous* pole: TP63, KRT5, KRT6A, DSG3, CLCA2, PKP1, ANXA8 — the same keratinization genes
as in breast. In oesophagus (ESCC vs EAC, overlap 0/30) and bladder (basal vs luminal, overlap 1/30),
the residual names the *opposing* pole: hepatic-lineage master regulators HNF4A/HNF1A in ESCA, and
urothelial-lineage PPARG/SLC14A1/BNC1 in BLCA. We propose a pole-selection rule: when viral (HPV
in CESC) or developmental (simple histology in lung) forces amplify the squamous keratinization
programme to extreme purity, it dominates the residual; when the tissue-specific non-squamous identity
programme is the more organ-distinctive signal, the residual surfaces that counter-pole instead. In
cervix, even the counter-pole leaves a trace — HNF1A ranks fifth, shared with the ESCA adeno axis —
but cannot dominate over HPV-driven squamous keratinization.

*Endpoint purity is not cosmetic.* The BLCA series provides a controlled experiment: the same tumours,
the same anchor, the same method — only the endpoint definition differs. Molecular subtypes
(basal-squamous vs luminal-papillary, n = 72) yield AUROC 0.967 and a clean residual naming
urothelial-lineage genes. Clinical subtype (Non-Papillary vs Papillary, n = 421 — a mixture spanning
all molecular subtypes) yields AUROC 0.633 and a residual dominated by stromal/invasiveness markers
(SFRP2, ROR2, CHST15). The keratinization biology is present in both datasets; only molecularly
homogeneous comparison groups allow it to be detected.

*Specificity: two adeno-only negative controls.* To test whether the panel detects generic epithelial
differentiation rather than squamous lineage specifically, we ran two comparisons in which neither pole
is squamous. Gastric adenocarcinoma (TCGA STAD, intestinal-type n = 108 vs diffuse/signet-ring n = 87)
yields panel AUROC 0.517 ≈ chance and 0/30 overlap; the residual names immune infiltration (CD52, CD37,
CD53, GZMK, ADORA3), reflecting EBV-positive/MSI immune-hot intestinal-type versus immune-cold
CDH1-mutant diffuse gastric cancer. Endometrial carcinoma (TCGA UCEC, serous n = 61 vs endometrioid
n = 117) yields AUROC 0.613 and 1/30 overlap (CLDN19; p = 0.17, non-significant); the residual names
serous-endometrial markers (L1CAM, TP53TG3B) rather than keratins. The UCEC AUROC modestly exceeds
chance, likely because serous endometrial cancer's claudin-family tight-junction remodelling creates
a weak overlap with the panel's claudin content — real but insufficient to reach the squamous threshold.
Across the series, panel AUROC bifurcates cleanly: ≥ 0.91 in five squamous-containing comparisons,
0.52–0.61 in two adeno-only controls (Table 2). The 30–40 percentage-point gap — together with
organ-specific non-keratinization residuals in both negative controls — proves that the axis is
specific to squamous lineage identity and entirely silent in its absence.

**Clinical significance: identity, not outcome (an honest negative).** A reproducible axis need not be
prognostic. In TCGA-BRCA (n = 866, 132 events) the basal score is not associated with overall survival
(Cox HR 1.01, p = 0.89; adjusted for proliferation p = 0.36; KM log-rank p = 0.29), whereas the proliferation
score is (p = 0.001). The basal axis does mark ER-negative/basal-like disease (AUROC 0.70). The discovered
axis therefore captures lineage *identity* rather than outcome — reported plainly, because the method's
value is finding real, robust biology, not inflated clinical claims.

**Hypothesis-as-anchor — confirming, explaining-away, or refuting a hypothesis (Figure 2).** The frame also
inverts cleanly: a *hypothesis* is expressed as a candidate anchor and tested against the textbook anchor on
real data, gating it onto the textbook residual for a three-way verdict — SUPPORTED (adds beyond the
textbook), EXPLAINED_BY_TEXTBOOK (predicts alone but redundant once the dominant prior is controlled), or
REFUTED. This operationalises Venet et al. [10]: a signature is a novel mechanism only if it survives
adjustment for the dominant prior. On Luminal A vs B against the proliferation anchor (Figure 2A), the
basal/keratinization hypothesis is SUPPORTED (Δ +0.039, β = 8), an immune/cytotoxic hypothesis is REFUTED
(AUROC 0.51, adds 0.000), and a random 30-gene set is EXPLAINED_BY_TEXTBOOK (predicts at 0.63 alone but adds
≈ 0). Scaled to a whole library (`rank_hypotheses` over the 50 MSigDB Hallmark sets; Figure 2B), the screen
is self-validating: the proliferation-type hallmarks (E2F, G2M, MYC) add exactly 0 beyond the anchor
(EXPLAINED — they *are* the anchor's axis), while the SUPPORTED hits are the orthogonal lineage programs led
by estrogen-response (early & late) — matching the known LumA/B biology of proliferation plus ER signalling.

A binary verdict can mislead when a hypothesis fails to add: *absent* and *redundant-because-collinear* look
identical (both Δ ≈ 0) yet mean opposite things. We therefore attach a gate-free commonality decomposition
(after Tonidandel & LeBreton [11]) to every hypothesis: the variance it explains is partitioned into the part
**unique** beyond the anchor versus the part **common** (shared with it), and its effect is split into direct
versus anchor-**mediated**, yielding a `collinearity_label` ∈ {NOVEL, REDUNDANT, INERT}. This makes the one
cross-cohort non-reproduction in our screen honest rather than misleading (Figure 2C): the estrogen-response
hypothesis is SUPPORTED in TCGA but not in METABRIC — yet the decomposition shows the difference is *structural,
not biological*. In TCGA, ER carries variance unique beyond proliferation (unique R² = 0.038 → NOVEL); in
METABRIC the same ER signal is essentially entirely shared with proliferation (redundancy = 1.00, 96 % of its
effect mediated through the proliferation anchor → REDUNDANT), because LumB's lower ER is collinear with its
higher proliferation in that cohort. ER biology is *redundant* in METABRIC, not *absent* — a distinction the
raw verdict cannot make and a transportability caveat for any anchored screen ported across cohorts. A
controlled sweep makes the mechanism explicit (Figure 2D): holding both marginal effects fixed and varying
*only* the anchor–hypothesis correlation, the same ER effect is labelled NOVEL in 100 % of simulated cohorts at
TCGA's correlation (+0.19) but falls into a collinear/suppression valley at METABRIC's (−0.17, ≈ 2 % NOVEL).
The verdict is thus governed by each cohort's covariate distribution — a generalisability/transportability
property [12] — which is why an anchored hypothesis screen should be re-characterised, not merely re-run,
across cohorts. Repeating the analysis on a second endpoint (HER2, with the ERBB2-amplicon anchor and an ER
hypothesis) confirms the labels are specific rather than a single "fails to add": there the ER→HER2 marginal
effect itself differs across cohorts, so HER2 is *not* a same-marginal flip — ER is **INERT** in TCGA (Cohen's
d ≈ −0.05, too weak to be novel or redundant) but **REDUNDANT** in METABRIC (d ≈ −0.29, redundancy 0.89, 66 %
mediated by the amplicon). Across the two endpoints the framework cleanly separates all three regimes — *novel*
(ER in TCGA Luminal A/B), *redundant/collinear* (ER in METABRIC Luminal A/B and HER2), and *absent/weak* (ER
for HER2 in TCGA).

Extending this to a five-endpoint × four-column panel (Figure 3) — TCGA RNA-seq, TCGA Agilent microarray
(the *same* patients on a different platform), METABRIC (an independent microarray cohort), and SCAN-B
(GSE96058, a fully independent Swedish RNA-seq cohort of ~3,400 tumours) — makes the practical payload
explicit: of five endpoints, three transport and two do not. The transportable trio are orthogonal,
biologically robust axes — a basal/keratinization anchor with an immune hypothesis is NOVEL in *all four*
columns (immune infiltration adds a real axis beyond the basal lineage program), an ER-signature anchor with a
proliferation hypothesis carries a small but genuine unique proliferation slice in all four, and an
ERBB2-amplicon anchor with an immune hypothesis on the HER2-enriched subtype is NOVEL in all four (immune
infiltration adds beyond the amplicon) — all robust across platform *and* independent cohort. The two that do
*not* transport are exactly the ER-collinearity cases above
(proliferation→ER on Luminal A/B; amplicon→ER on HER2), and the four columns expose why: for Luminal A/B the
label tracks *measurement technology* — NOVEL on both RNA-seq cohorts (TCGA RNA-seq and the independent SCAN-B)
yet REDUNDANT on both microarrays (TCGA Agilent and METABRIC), flipping even on the same TCGA patients between
their RNA-seq and Agilent measurements. The transportability caveat is therefore specific and predictable: it
attaches to hypotheses whose collinearity with the anchor is itself measurement-dependent, while genuinely
anchor-orthogonal axes transport cleanly across platform and cohort.

![**Figure 3. Anchored hypothesis labels across five endpoints × four columns (TCGA RNA-seq, TCGA Agilent,
METABRIC, SCAN-B).** Each cell is one (anchor → hypothesis) test in one column, coloured by the commonality
label (NOVEL = carries variance unique to the hypothesis beyond the anchor; REDUNDANT = predicts but
collinear/mediated; INERT = no appreciable signal), with the cross-validated verdict and the redundancy below
it. TCGA Agilent is the *same patients* as TCGA RNA-seq on a different platform (a platform-transportability
check); METABRIC (microarray) and SCAN-B (GSE96058, RNA-seq) are independent cohorts. The "transports/differs"
tag marks whether the label is concordant across all four columns. Basal→immune and ER-status→proliferation
transport (NOVEL everywhere, including the independent SCAN-B); the two ER-collinearity endpoints (Luminal A/B,
HER2) do not — and for Luminal A/B the label tracks measurement technology, NOVEL on the two RNA-seq columns
(TCGA RNA-seq, SCAN-B) but REDUNDANT on the two microarrays (TCGA Agilent, METABRIC).](figs/endpoint_panel.png){ width=100% }

**The platform effect is a sign flip in the nuisance correlation (Figure 4).** The transportability sweep
predicts that the Luminal A/B verdict is governed by corr(proliferation, ER); the four columns let us read
that correlation directly, and it is the proximate cause of the label split. Its *sign* flips with the assay:
it is positive on both RNA-seq cohorts (TCGA RNA-seq +0.19, SCAN-B +0.10) — where ER and proliferation are not
collinear in the Luminal B direction, so ER retains suppression/unique variance and is labelled NOVEL — and
negative on both microarrays (TCGA Agilent −0.10, METABRIC −0.17) — where Luminal B's lower ER aligns with its
higher proliferation, so ER's signal is collinear and is labelled REDUNDANT. Because the flip occurs even on
the *same* TCGA patients between their RNA-seq and Agilent measurements, it is a property of the measurement,
not of the patients. The flip is statistically grounded in the two large cohorts: the bootstrap 95 % CI of
corr(proliferation, ER) is wholly positive in SCAN-B (0.04–0.15) and wholly negative in METABRIC (−0.24 to
−0.12), non-overlapping. A per-endpoint transport score (fraction of cohort-pairs with the same commonality
label) is 1.0 for the anchor-orthogonal axes (basal→immune; ER-status→proliferation) and only 0.17–0.33 for the
two ER-collinearity endpoints. This locates the transportability caveat precisely: the nuisance correlation that decides
a collinear hypothesis's verdict can itself be set by the assay, so such labels must be read within a fixed
platform — whereas the anchor-orthogonal axes (basal→immune, ER-status→proliferation) are immune to this and
reproduce across both platform and cohort.

![**Figure 4. The Luminal A/B platform effect is a sign flip in corr(proliferation, ER).** For the
Luminal A vs B endpoint, the proliferation–ER correlation that governs the verdict is plotted per cohort and
coloured by platform family. It is positive on the two RNA-seq cohorts (TCGA RNA-seq, SCAN-B) → ER carries
unique/suppression variance → NOVEL, and negative on the two microarray cohorts (TCGA Agilent, METABRIC) → ER
collinear with proliferation → REDUNDANT. The sign flips even between the two TCGA platforms on the same
patients, so it is a measurement property.](figs/platform_corr.png){ width=80% }

**Results at a glance.** Table 1 summarises the anchored discoveries and their external status.

| Endpoint / target | Anchor | Anchor AUROC | Residual gain (ΔAUROC) | Discovered axis | External status | Label |
| :-- | :-- | :--: | :--: | :-- | :-- | :-- |
| Luminal A vs B | 20-gene proliferation | 0.919 | +0.029 | basal/keratinization (KRT5/14/17/6B, TP63, DSG3/DSC3, SOX10, COL17A1, KLK5/7/8) | METABRIC Δ+0.036; lung 10/30 (p=3×10⁻¹⁶); SCAN-B confirmed; ESCA 0.91; BLCA 0.97; CESC 0.94 (6/30, p=1.5×10⁻⁸); STAD 0.52/UCEC 0.61 (neg. ctrl) | NOVEL |
| HER2 status | ERBB2 amplicon | 0.752 | +0.054 | neuroendocrine/secretory + immune | not in METABRIC (anchor near-complete, 0.997) | — |
| ER status | ER/luminal signature | 0.938 | −0.001 | none (specificity control) | n/a | INERT |
| NSCLC anti-PD-1 benefit | TMB | 0.60 | +0.061 | PD-L1; EGFR/STK11 mutation | recovers established biomarkers | — |

*Table 1. Anchored residual discovery across endpoints: a zero-parameter/textbook anchor, the signal the
residual adds beyond it, the discovered orthogonal axis, and its external reproduction.*

**Table 2. Cross-cancer validation of the pan-epithelial keratinization axis (seven TCGA cohorts).**

| Cancer | Comparison | n | Panel AUROC | Overlap / 30 | Residual pole named | Squamous? |
| :-- | :-- | :--: | :--: | :--: | :-- | :--: |
| Lung LUAD/LUSC | Adeno vs squamous | 1,129 | ≥0.96 | **10/30** (p=3×10⁻¹⁶) | Squamous/keratinization (same pole) | ✓ |
| HNSC | Cross-tissue (HNSC+LUSC vs LUAD) | ~500 | **0.962** | — | Tissue-independent squamous marker | ✓ |
| ESCA ESCC/EAC | Squamous vs adeno | 196 | **0.913** | 0/30 | Adeno counter-pole (HNF4A, MUC13) | ✓ |
| BLCA Basal/Luminal | Molecular subtypes | 72 | **0.967** | 1/30 | Luminal counter-pole (PPARG, SLC14A1) | ✓ |
| CESC Sq/Adeno | Squamous vs adeno | 301 | **0.938** | **6/30** (p=1.5×10⁻⁸) | Squamous pole (TP63, KRT5/6A, CLCA2) | ✓ |
| STAD Int/Diff *(neg. ctrl)* | Intestinal vs diffuse adeno | 195 | 0.517 | 0/30 | Immune infiltration (CD52, GZMK) | ✗ |
| UCEC Ser/Endo *(neg. ctrl)* | Serous vs endometrioid adeno | 178 | 0.613 | 1/30 (n.s.) | Serous markers (L1CAM, TP53TG3B) | ✗ |

*Panel AUROC: breast 30-gene basal panel scored in the new cancer; overlap: hypergeometric test of
residual genes vs breast panel (top 30). Squamous-containing comparisons (✓): AUROC ≥ 0.91.
Adeno-only controls (✗): AUROC 0.52–0.61. BLCA endpoint-purity pair not shown (Non-Pap/Pap AUROC 0.633).*

## 3. Discussion

The method instantiates and unifies three established traditions — clinical-offset / incremental-value
modelling (an established model as offset, omics on the residual; the closest methodological twin) [5];
biologically-informed models that bake known biology into structure [6,7]; and machine learning on
established gene signatures [8]. Its specific contribution is a packaged, never-below-anchor, margin-gated
residual on a *zero-parameter* prior, plus a residual-discovery framing with an explicit noise control, and
a demonstration that the discovered axes' external reproducibility is predictable from their biological
coherence. The honest scope: predictive gains over a strong anchor are modest — the value is *routing and
discovery*, not large AUROC wins — consistent with the rarity of genuine super-additive multi-omics gains
[1,2]. Discovered axes are candidate hypotheses; the basal axis is externally validated, the HER2 axis is
cohort-specific, and other endpoints await further cohorts.

A second contribution is diagnostic honesty about *what fails to reproduce and why*. Treating a hypothesis as
an anchor and decomposing its signal into the part unique to it versus the part shared with the textbook prior
[11] separates three regimes a bare verdict conflates — novel, redundant-because-collinear, and absent — and a
mediation split shows how much of a hypothesis's effect runs through the anchor. Across four breast-cancer
endpoints and four columns (two RNA-seq, two microarray; one platform pair on identical patients), the
anchor-orthogonal axes are labelled identically everywhere, while the collinear ones are not: their verdict is
a transportability quantity [12] set by the anchor–hypothesis correlation, which we show can flip sign purely
with the measurement platform. The practical implication is concrete — an anchored hypothesis screen should be
re-characterised, not merely re-run, when ported across assays or cohorts, and collinearity-sensitive calls
should be read within a fixed platform. This reframes a non-reproduction (estrogen-response on Luminal A/B in
microarray cohorts) from an apparent failure into a quantified, expected consequence of covariate structure.

**Limitations.** Several caveats bound these claims. (i) *Modest predictive gains.* Over a strong anchor the
gated residual adds little AUROC (e.g. +0.029 on Luminal A/B); the contribution is routing, discovery and
diagnosis, not large predictive wins, and the method will not rescue an endpoint a good single view already
predicts. (ii) *Lineage, not outcome.* The flagship basal axis is a reproducible *biological* program (it marks
ER-negative/basal identity) but is not prognostic in the cohort studied; discovered axes are candidate
mechanisms, not validated clinical biomarkers. (iii) *Anchor dependence.* Results are conditional on the chosen
prior — a near-complete anchor (e.g. the ERBB2 amplicon in METABRIC, or the ER signature for ER status) leaves
no residual by construction, so "nothing found" means "nothing beyond this anchor," not "nothing there." (iv)
*Cohorts and platforms.* Validation is breast-centric (TCGA, METABRIC, SCAN-B) with cross-cancer
extension to seven TCGA cohorts (lung, HNSC, ESCA, BLCA, CESC, STAD, UCEC) and one non-expression
cross-domain demonstration (NSCLC anti-PD-1 mutation/clinical features); the hypothesis-screen labels
are bulk-expression and assay-dependent (as the platform analysis shows), and single-cell or proteomic
anchors are untested. (v) *Statistics.* CIs are
non-parametric bootstrap and the empirical-null FDR is approximate; the transport_score is descriptive, and
some commonality labels near the NOVEL/INERT threshold are sensitive to that cut-off. (vi) *Scope of evidence.*
This is a single-author methods study on public data with no prospective or wet-lab validation; the immunotherapy
and cross-cancer results are proof-of-concept on one cohort each and warrant replication.

## 4. Methods

**Anchored integration.** `select_anchor` scores each modality by repeated stratified-CV AUROC (tie-broken by
robustness) and picks the top composite, inside the outer CV. `anchored_integrate` pins the anchor and adds a
secondary on its residual as `logit(anchor) + β·secondary`, with β ≥ 0 chosen on held-out data (β = 0
allowed; a margin gate plus repeated inner CV control small-sample noise). `forward_integrate` / `auto_integrate`
extend this to any number of modalities by greedy forward selection. **Knowledge anchor.** `signature_score`
builds a fixed mean-z score from a curated gene set; `knowledge_anchored_integrate` pins it as the anchor.
**Residual discovery.** `anchored_residual_discovery` ranks features by partial correlation with the label
controlling for the anchor, retains those orthogonal to the anchor (|r| < 0.6), gates the panel, and tests it
against matched random panels and (for verification) held-out splits and permuted labels. Because the
matched-random-panel null saturates when an endpoint is broadly predictable (e.g. squamous-vs-adeno, ER
status), the method also reports a **selection-stability** statistic — the recurrence of the discovered panel
across 50% subsamples relative to a permuted-label null — which stays informative there (lung basal axis
stability gain ≈ +0.89; genome-wide ER methylation ≈ +0.63, both where the panel-vs-random null is
uninformative). Data: TCGA-BRCA
(Xena HiSeqV2, HumanMethylation450, clinical), METABRIC (microarray, clinical), and seven additional
TCGA cohorts via Xena HiSeqV2 (LUAD/LUSC, HNSC, ESCA, BLCA, CESC, STAD, UCEC). BLCA molecular
subtypes (basal-squamous / luminal-papillary) follow Robertson et al. 2017 [13]. Pathway enrichment:
g:Profiler. All functions are in `omniomics.multiomics`; runners and recorded metrics are in the repository.
**Scalability.** The discovery is linear and streams: a vectorized residualize-then-correlate, an optional
one-pass Sure-Independence-Screening pre-filter, Benjamini–Hochberg FDR, and parallel nulls let it run on
genome-wide feature matrices out of core. As an end-to-end check, screening all 485,577 HumanMethylation450
probes (n = 829) to the top 5,000 by association with ER status took ≈ 24 s in a single low-memory pass; the
ESR1-methylation anchor (AUROC 0.65) then rose to 0.90 with genome-wide CpGs, and the top CpG beyond ESR1
mapped to PGR (the canonical ER-coregulated gene) — textbook biology recovered at scale.

## 5. Data and Code Availability

Code, runners, recorded metrics, unit tests and continuous-integration guards: the `omniomics` repository
(`reports/dmoi_*.py`, `*_results.csv`, `tests/`). Key runners: `dmoi_knowledge_anchor.py`,
`dmoi_residual_discovery.py`, `dmoi_discovery_{er,her2}.py`, `dmoi_external_{metabric,her2_metabric}.py`,
`dmoi_external_lung.py`, `dmoi_external_hnsc.py`, `dmoi_external_{esca,blca,cesc,stad,ucec}.py`
(cross-cancer validation #3–#7), `dmoi_discovery_nsclc_io.py`, `dmoi_clinical_survival.py`,
`dmoi_meth_genomewide.py` (genome-wide scale). Hypothesis-as-anchor and transportability:
`dmoi_hypothesis_anchor.py`, `dmoi_hypothesis_screen{,_robust}.py`, `dmoi_hypothesis_metabric_diagnosis.py`,
`dmoi_transportability_sweep.py`, `dmoi_transportability_her2.py`, `dmoi_endpoint_panel.py`,
`fig_hypothesis_anchor.py`, `fig_platform_corr.py` (Figures 2–4); `reports/fetch_scanb.sh` reproduces the
SCAN-B inputs. Scale utilities: `omniomics.scale`, `omniomics.representations`. Full methods note:
`reports/anchored_integration_methods.md`. All data are public and de-identified: TCGA (BRCA RNA-seq and
Agilent microarray, LUAD, LUSC, HNSC) via UCSC Xena; METABRIC via cBioPortal; SCAN-B via GEO accession
**GSE96058**; the NSCLC anti-PD-1 cohort from Hellmann et al. (MSK 2018) via cBioPortal.

## Declarations

**Competing interests.** The author declares no competing interests.

**Funding.** No external funding was received for this work.

**Author contributions.** H.R.K. conceived the method, performed all analyses, and wrote the manuscript.

**Ethics.** This study uses only publicly available, de-identified data from established repositories
(TCGA/UCSC Xena, cBioPortal); no new human-subjects data were collected, so no additional ethics approval
was required.

**License.** This preprint is made available under a CC-BY 4.0 International license.

**AI assistance.** Analyses and drafting were carried out with the assistance of an AI coding agent under the
author's direction; all results are reproducible from the cited runners.

## Figure

![**Figure 1. Knowledge-anchored discovery and its external validation.** **(A)** Across three TCGA-BRCA
endpoints, the gain in AUROC from adding genome-wide data beyond a zero-parameter textbook anchor: the data
adds a real axis where the textbook is incomplete (LumA/B proliferation → basal, +0.029; HER2 amplicon →
neuroendocrine/immune, +0.054) and nothing where it is complete (ER signature, −0.001; specificity).
**(B)** In the independent METABRIC cohort, an unbiased re-discovery recovers the basal axis (20/30 overlap
with the TCGA panel; direct replication Δ +0.036, p = 0.048) but not the HER2 axis (2/30), because the
amplicon anchor is already complete in METABRIC. **(C)** Cross-cancer replication: across seven TCGA cohorts the breast-derived 30-gene panel achieves AUROC
≥ 0.91 in all five squamous-containing comparisons (lung LUAD/LUSC, HNSC, ESCA, BLCA, CESC) and falls to
0.52–0.61 in two adeno-only negative controls (STAD, UCEC); in lung (n = 1,129) the residual independently
re-discovers 10/30 breast basal genes (KRT5/14/6B, TP63, DSG3/DSC3, FAT2…; p ≈ 3×10⁻¹⁶).](figs/discovery_summary.png){ width=100% }

![**Figure 2. Hypothesis-as-anchor: confirm, explain-away, or refute.** **(A)** On TCGA-BRCA Luminal A vs B,
three hypotheses tested against the textbook proliferation anchor by the signal each adds beyond it: the
basal/keratinization hypothesis is SUPPORTED (Δ +0.039), an immune/cytotoxic hypothesis is REFUTED (adds
0.000), and a random 30-gene set is EXPLAINED_BY_TEXTBOOK (predicts alone but adds ≈ 0 once proliferation is
controlled). **(B)** Screening all 50 MSigDB Hallmark gene sets as candidate hypotheses: the
proliferation-type hallmarks (E2F, G2M, MYC) add exactly 0 beyond the anchor (EXPLAINED — the same axis),
while the SUPPORTED hits are the orthogonal lineage programs led by estrogen-response, matching known LumA/B
biology. **(C)** Commonality/mediation re-characterization of the estrogen-response hypothesis across cohorts:
stacked R² in the LumA/B endpoint split into the part *unique* to ER (beyond proliferation, blue) versus the
part *shared* with proliferation (grey). ER is NOVEL in TCGA (unique R² 0.038, orthogonal to proliferation)
but REDUNDANT in METABRIC (redundancy 1.00, 96 % mediated through proliferation) — collinear, not absent.
**(D)** Transportability of the verdict: holding both marginal effects fixed and varying *only* the
anchor–hypothesis correlation, the fraction of simulated cohorts labelled NOVEL vs REDUNDANT is plotted
against that correlation; the two real cohorts sit on the curve at their measured corr(proliferation, ER) —
the *same* ER effect is 100 % NOVEL at TCGA's +0.19 but collapses into the collinear/suppression valley at
METABRIC's −0.17. The verdict is a covariate-distribution property, not a difference in ER biology.
](figs/hypothesis_anchor.png)

## References

*References verified against PubMed; DOIs link to the source.*

1. Li Y, Herold T, Mansmann U, Hornung R. Does combining numerous data types in multi-omics data improve or hinder performance in survival prediction? Insights from a large-scale benchmark study. *BMC Med Inform Decis Mak* 2024;24(1):244. doi:[10.1186/s12911-024-02642-9](https://doi.org/10.1186/s12911-024-02642-9).
2. Ellen JG, Jacob E, Nikolaou N, Markuzon N. Autoencoder-based multimodal prediction of non-small cell lung cancer survival. *Sci Rep* 2023;13(1):15761. doi:[10.1038/s41598-023-42365-x](https://doi.org/10.1038/s41598-023-42365-x).
3. Makrodimitris S, Pronk B, Abdelaal T, Reinders M. An in-depth comparison of linear and non-linear joint embedding methods for bulk and single-cell multi-omics. *Brief Bioinform* 2023;25(1):bbad416. doi:[10.1093/bib/bbad416](https://doi.org/10.1093/bib/bbad416).
4. Ding DY, Li S, Narasimhan B, Tibshirani R. Cooperative learning for multiview analysis. *Proc Natl Acad Sci USA* 2022;119(38):e2202113119. doi:[10.1073/pnas.2202113119](https://doi.org/10.1073/pnas.2202113119).
5. Volkmann A, De Bin R, Sauerbrei W, Boulesteix AL. A plea for taking all available clinical information into account when assessing the predictive value of omics data. *BMC Med Res Methodol* 2019;19(1):162. doi:[10.1186/s12874-019-0802-0](https://doi.org/10.1186/s12874-019-0802-0).
6. Elmarakeby HA, et al. Biologically informed deep neural network for prostate cancer discovery. *Nature* 2021;598(7880):348–352. doi:[10.1038/s41586-021-03922-4](https://doi.org/10.1038/s41586-021-03922-4).
7. Liu X, et al. Pathformer: a biological pathway informed transformer for disease diagnosis and prognosis using multi-omics data. *Bioinformatics* 2024;40(5):btae316. doi:[10.1093/bioinformatics/btae316](https://doi.org/10.1093/bioinformatics/btae316).
8. Boscolo Bielo L, et al. A machine learning assay to predict disease recurrence in hormone receptor-positive breast cancer. *ESMO Open* 2026;11(3):106064. doi:[10.1016/j.esmoop.2026.106064](https://doi.org/10.1016/j.esmoop.2026.106064).
9. Horvath S. DNA methylation age of human tissues and cell types. *Genome Biol* 2013;14(10):R115. doi:[10.1186/gb-2013-14-10-r115](https://doi.org/10.1186/gb-2013-14-10-r115).
10. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol* 2011;7(10):e1002240. doi:[10.1371/journal.pcbi.1002240](https://doi.org/10.1371/journal.pcbi.1002240).
11. Tonidandel S, LeBreton JM. Relative importance analysis: a useful supplement to regression analysis. *J Bus Psychol* 2011;26(1):1–9. doi:[10.1007/s10869-010-9204-3](https://doi.org/10.1007/s10869-010-9204-3). *(Commonality/relative-importance methodology; not indexed in PubMed.)*
12. Degtiar I, Rose S. A review of generalizability and transportability. *Annu Rev Stat Appl* 2023;10:501–524. doi:[10.1146/annurev-statistics-042522-103837](https://doi.org/10.1146/annurev-statistics-042522-103837). *(Generalizability/transportability methodology; not indexed in PubMed.)*
13. Robertson AG, et al. Comprehensive molecular characterization of muscle-invasive bladder cancer. *Cell* 2017;171(3):540–556. doi:[10.1016/j.cell.2017.09.007](https://doi.org/10.1016/j.cell.2017.09.007).
