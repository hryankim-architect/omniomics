---
title: "Anchored Multi-Omics Integration and Knowledge-Anchored Residual Discovery"
subtitle: "A never-below-the-best-view integrator, and a tool to find what the known biology misses — a white paper"
author: "H. Ryan Kim"
date: "2026"
---

*A white paper for a technical-but-general audience. It distils the manuscript and its companion NotebookLM
materials (briefing, study guide, data table, slide deck, infographic, and the transcribed audio overview and
explainer video) into a single accessible document. All quantitative claims match the recorded result files and
the manuscript.*

## Executive summary

Combining many kinds of molecular data is widely assumed to improve cancer prediction, yet in practice a single
strong data type is hard to beat and naive fusion often does worse. This white paper presents an **anchored
integrator**: it starts from established biological knowledge — a known high-performing data type, or a fixed,
zero-parameter gene signature — and admits genome-wide data only as a non-negative gated residual, so the model
is *provably never below its anchor* and improves on it only where another view carries orthogonal signal.

Three ideas follow from that one move. First, mining the **residual** (what the anchor cannot explain) turns the
integrator into a discovery engine — a targeted search for the biology the known prior misses. Second, running
it **in reverse** lets a hypothesis be tested against the prior and labelled novel, already-explained, or
absent. Third, the framework makes **transportability** explicit: it shows precisely when, and why, the same
biological signal looks novel in one cohort and redundant in another — often a property of the measurement
technology rather than the biology. The result is a discovery engine that is honest about what is new and what
the known biology already explains.

## 1. The problem: more data is not automatically better

The prevailing instinct in bioinformatics is to gather more — more scans, deeper sequencing, more molecular
layers — on the assumption that fusing them yields clearer answers. The data say otherwise. A large benchmark
across 14 cancers found that naively concatenating data types most often *hurt* survival prediction: stacking
layers adds statistical variance that drowns the clear signal in a strong modality, and the model spends its
capacity fitting noise. The honest objective is therefore not "always fuse," but *never do worse than the best
single view, and improve on it only where another view carries signal the leader cannot.*

A kitchen analogy from the companion audio captures it: a great broth does not come from tipping every leftover
into the pot; you start from a strong foundation and add an ingredient only if it genuinely improves the
flavour.

## 2. The method: anchor, gate, mine the residual

- **Anchor** on the empirically strongest modality, or on a fixed **zero-parameter** prior — a textbook
  gene-signature score that has had *no* training on the current dataset.
- **Gate.** Admit secondary data as `logit(anchor) + β·secondary`, with the weight `β ≥ 0` chosen on held-out
  data. If the new layer is redundant or noisy, the gate stays shut (`β = 0`). This guarantees the model is
  never below its anchor and yields the framework's golden rule: *never below the best single view.*
- **Mine the residual.** Rank features by partial correlation with the outcome while controlling for the
  anchor, surfacing **anchor-orthogonal** axes, with a matched-random-panel noise control and held-out /
  permutation verification.

On TCGA breast cancer, blind anchor selection routed correctly (methylation for methylation-defined clusters,
RNA for PAM50 calls), never falling below the leader; a constructed positive control (signal invisible to RNA)
engaged the gate for +0.047 AUROC. The companion explainer frames the gate as a veteran clinician's established
knowledge, with genome-wide data an eager newcomer allowed to speak only when it adds something genuinely new.

## 3. Discovery: a verified basal/keratinization axis

Luminal A vs B breast cancer is, by convention, a proliferation distinction. A **20-gene proliferation index
(zero trained parameters)** reaches AUROC **0.919** — close to a fully trained 1,500-gene model (0.942) — and
gating genome-wide data onto this fixed prior reaches **0.947**, the cleanest super-additive gain in the study.
Strikingly, a textbook rule with no training nearly matches a 1,500-gene model; the framework is, in effect,
teaching the model to read the textbook first so its compute hunts only for what is unknown.

Mining that prior's residual surfaces a **basal / squamous-lineage axis** (KRT5/14/17/6B, TP63, DSG3/DSC3,
SOX10, COL17A1, KLK5/7/8; keratinization, Reactome p ≈ 8×10⁻¹¹), verified by held-out replication (10/10),
recurrence of the core panel (10/10), and collapse under permuted labels. It reproduces externally in METABRIC
(Δ +0.036; 20/30 gene overlap, hypergeometric p ≈ 7×10⁻²⁷), is confirmed in the independent SCAN-B cohort
(Section 6), and recovers the same squamous/keratinization program in TCGA lung cancer (LUAD vs LUSC; 10/30
overlap, p ≈ 3×10⁻¹⁶) — the same genes in a different organ.

## 4. Generalisation across anchors, diseases, and feature types

- **A second knowledge anchor (HER2).** Anchoring on the ERBB2 amplicon (incomplete in TCGA, 0.752) discovers
  a verified neuroendocrine/secretory + immune axis (Δ +0.054); it does not reproduce in METABRIC because the
  amplicon anchor is near-complete there (0.997), leaving no residual — an interpretable negative.
- **A specificity control (ER).** Anchoring ER status on a complete ER/luminal signature (0.938) discovers
  nothing (Δ −0.001): the method refuses to manufacture an axis when the prior already suffices — an honest
  negative, the opposite of a hallucination.
- **A different disease and feature type (NSCLC immunotherapy).** Anchored on the textbook biomarker tumour
  mutational burden (TMB) and given mutation/clinical features, the residual independently recovers the field's
  other established biomarkers — PD-L1 (orthogonal to TMB), and EGFR and STK11 mutation (resistance) — Δ +0.061
  vs +0.006 for random panels (p = 0.038).

## 5. Hypothesis-as-anchor: confirm, explain-away, or refute

Run in reverse, a candidate hypothesis is gated onto the textbook anchor's residual for a three-way verdict —
**SUPPORTED** (adds beyond the prior), **EXPLAINED_BY_TEXTBOOK** (predicts alone but redundant once the prior
is controlled), or **REFUTED** — operationalising a mathematical form of peer review. A gate-free **commonality
decomposition** (after Tonidandel & LeBreton) then partitions a hypothesis's explained variance into the part
*unique* to it versus the part *shared* with the anchor, and a mediation split shows how much of its effect runs
through the anchor. This yields a label — **NOVEL**, **REDUNDANT** (collinear), or **INERT** — that separates a
signal that is *absent* from one that is merely *redundant*. A 50-set hallmark screen is self-validating:
proliferation-type hallmarks add ≈ 0 beyond the proliferation anchor, while estrogen-response programs are the
supported orthogonal hits.

## 6. Transportability: what fails to reproduce, and why

The framework's diagnostic contribution is to explain non-reproduction rather than merely report it. The
estrogen-response hypothesis is **NOVEL** in TCGA (unique R² = 0.038) but **REDUNDANT** in METABRIC (redundancy
= 1.00, 96 % mediated through proliferation) — collinear, not absent. A controlled sweep (hold both marginal
effects fixed, vary only the anchor–hypothesis correlation) reproduces this: the same ER effect is 100 % NOVEL
at TCGA's correlation (+0.19) yet collapses into a collinear/suppression valley at METABRIC's (−0.17). A second
endpoint, HER2, confirms specificity (INERT in TCGA, REDUNDANT in METABRIC).

Across four endpoints and four columns — TCGA RNA-seq; TCGA Agilent (the *same patients*, different platform);
METABRIC (independent microarray); and SCAN-B (independent RNA-seq, ~3,400 tumours) — two endpoints transport
and two do not:

| Endpoint (anchor → hypothesis) | TCGA RNA-seq | TCGA Agilent | METABRIC | SCAN-B | Transports? |
| :-- | :-- | :-- | :-- | :-- | :-- |
| Basal vs rest (basal → immune) | NOVEL | NOVEL | NOVEL | NOVEL | ✓ |
| ER status (ER → proliferation) | NOVEL | NOVEL | NOVEL | NOVEL | ✓ |
| LumA vs LumB (proliferation → ER) | NOVEL | REDUNDANT | REDUNDANT | NOVEL | ✗ |
| HER2 (amplicon → ER) | INERT | INERT | REDUNDANT | NOVEL | ✗ |

For Luminal A/B the label tracks **measurement technology**: corr(proliferation, ER) is positive on both
RNA-seq cohorts (+0.19, +0.10 → NOVEL) and negative on both microarrays (−0.10, −0.17 → REDUNDANT), flipping
even on the *same* TCGA patients between their RNA-seq and Agilent measurements. As the companion media puts it,
the same tumour viewed through "RNA-seq glasses" versus "microarray glasses" casts an opposite statistical
shadow — so a model can end up *learning the platform, not the patient*. Genuinely anchor-orthogonal axes
(basal→immune; ER-status→proliferation) are immune to this and transport across both platform and cohort.

## 7. Implications and conclusion

Predictive gains over a strong anchor are modest by design; the value is **routing and discovery** plus
**diagnostic honesty**. The framework reports honest negatives (the basal axis marks lineage identity, not
survival, in the studied cohort), scales efficiently (485,577 methylation probes vs ER status in ~24 s,
recovering the canonical *PGR* gene), and shows that external reproducibility tracks biological coherence:
pathway-enriched, anchor-orthogonal axes reproduce across cohorts, platforms, and cancers, while collinear,
measurement-dependent signals do not — and the framework now says which is which.

The broader caution lands beyond any single result: if a calibrated model's verdict can flip with the assay,
some apparent replication failures across the literature may be platform artefacts rather than biological
disagreements. An anchored, residual-aware frame gives a rigorous way to tell the two apart.

*Core philosophy: anchor on established biology and let a gate decide, honestly, whether genome-wide data beats
it — and when a hypothesis fails to add, say whether it is novel-elsewhere, merely redundant, or truly absent.*

---

## Appendix A — Results at a glance

| Endpoint / target | Anchor | Anchor AUROC | Residual gain (ΔAUROC) | Discovered axis | External status | Label |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| Luminal A vs B | 20-gene proliferation | 0.919 | +0.029 | basal/keratinization (KRT5/14/17/6B, TP63, DSG3/DSC3, SOX10, COL17A1, KLK5/7/8) | METABRIC Δ+0.036; lung 10/30; SCAN-B confirmed | NOVEL |
| HER2 status | ERBB2 amplicon | 0.752 | +0.054 | neuroendocrine/secretory + immune | not in METABRIC (anchor near-complete 0.997) | — |
| ER status | ER/luminal signature | 0.938 | −0.001 | none (specificity control) | n/a | INERT |
| NSCLC anti-PD-1 benefit | TMB | 0.60 | +0.061 | PD-L1; EGFR/STK11 mutation | recovers established biomarkers | — |

## Appendix B — Companion media (NotebookLM Studio, transcribed)

Both generated from the same manuscript source and, on transcription, faithful to it (automatic-speech
slips corrected: PD-L1, STK11, LeBreton, gene names). Full transcripts accompany this paper.

- **Audio overview** (~24 min, two-host deep dive): data-overload imagery; the broth analogy; the gate as a
  "taste test at the door"; the basal discovery and its METABRIC/lung/NSCLC validation; the ER specificity
  control as an honest negative; and the transportability trap (the proliferation–ER correlation flipping sign
  on the *same patients* across RNA-seq vs microarray, framed as two coloured glasses).
- **Explainer video** (~7.6 min), five parts — the multi-omics myth, the knowledge anchor, mining the residual,
  cross-cancer discoveries, and the platform trap; memorable framings include "never below the best single
  view," a "search engine for the dark matter of the tumour's biology," and "learning the platform, not the
  patient."

## Appendix C — Data and code

All data are public and de-identified: TCGA (BRCA RNA-seq and Agilent, LUAD, LUSC, HNSC) via UCSC Xena;
METABRIC and the NSCLC anti-PD-1 cohort via cBioPortal; SCAN-B via GEO accession GSE96058. Code, reproduce
runners, recorded metrics, unit tests and CI guards live in the `omniomics` repository; the formal manuscript
and preprint accompany this white paper.
