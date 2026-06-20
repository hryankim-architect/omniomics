# Standardizing the "textbook anchor": problem, literature, and options

*Internal discussion memo — for reading + a joint decision. Not yet a committed methods change.*
*Sources verified via PubMed and Consensus (links/DOIs inline). Date: 2026-06.*

---

## 1. The problem (the critique)

In `anchored_residual_discovery`, the **predictive** guarantee (never-below-the-anchor) holds for *any*
anchor and is rigorous. But the **discovery** interpretation — "what beats the textbook" — is only as
principled as the anchor's provenance. In this project every anchor was **hand-picked**:

| Endpoint | Anchor used | Source type |
|---|---|---|
| LumA vs LumB | 20-gene proliferation signature | curated gene set |
| HER2 status | ERBB2 17q12 amplicon (9 genes) | locus / curated |
| ER status | ER/luminal signature | curated gene set |
| normal-tissue age | Horvath clock | published clock |
| NSCLC anti-PD1 | tumour mutational burden (TMB) | clinical biomarker |
| ER (methylation) | ESR1 promoter probes | single-gene locus |

Source, granularity, and curation differ wildly, and there is **no rule** mapping an endpoint to "the"
anchor. This creates three concrete risks:

1. **Researcher degrees of freedom (cherry-picking).** A chosen anchor can make a residual look real or
   absent; "objective textbook" framing hides a subjective choice.
2. **Non-reproducibility.** Two analysts may pick different anchors → different residuals → different
   "discoveries."
3. **Circularity.** A data-trained anchor smuggles fitting into a claimed "zero-parameter prior."

---

## 2. What the literature says

### 2a. The smoking gun — and, remarkably, the prescription
According to PubMed: **Venet, Dumont & Detours (2011), *PLoS Comput Biol* — "Most random gene expression
signatures are significantly associated with breast cancer outcome"**
([DOI: 10.1371/journal.pcbi.1002240](https://doi.org/10.1371/journal.pcbi.1002240)).
- Signatures *unrelated to cancer* (postprandial laughter; mouse social defeat) significantly predict breast
  cancer outcome. Of 47 published signatures, **60% were no better than random** sets of the same size;
  **>90% of random signatures >100 genes** were significant predictors.
- **The fix they propose is exactly our idea:** build a proliferation meta-gene (**meta-PCNA** = the top 1%
  of genes correlated with PCNA in a *normal-tissue* compendium) and **adjust for it** — doing so abrogates
  almost all signature–outcome association. I.e. *anchor on the dominant confounder, read the residual.*
- **Two gifts to us:** (i) meta-PCNA is a **reproducible recipe**, not hand-curation → a template for a
  standardized anchor; (ii) a warning — proliferation is so pervasive that **dropping cell-cycle genes does
  not remove it**; you must *regress out* a meta-gene (which is what our gate does), not just exclude genes.

### 2b. There is no single "universal" anchor
- **Michiels et al. (2005), *Lancet*** — signature gene lists are **highly unstable** and depend on the
  training-set sampling; 5 of 7 landmark studies classified no better than chance. Advocates **repeated
  random validation**.
  [Consensus](https://consensus.app/papers/details/5c64938c77ab54c1a9b313a773cfe511/)
- **Nagy et al. (2020), *Sci Rep*** — even cancer-hallmark genes show **tissue-specific, heterogeneous**
  prognostic effects → no universal anchor across cancers.
  [Consensus](https://consensus.app/papers/details/52cbbe35de8e5267aedd830e1e4ed6da/)

### 2c. The meta-problem is named and tooled
- **Multiverse / researcher degrees of freedom.** Del Giudice & Gangestad (2021) give the most useful
  frame: decide whether alternatives are *truly arbitrary* (Type E) or *principled-nonequivalent* (Type N) —
  and warn that lumping a **biologically justified** choice into an undifferentiated multiverse can **bury a
  real effect**. Our anchor choice is principled-nonequivalent, not arbitrary.
  [Del Giudice 2021](https://consensus.app/papers/details/42c0716ece6d552986802d64c1bd1608/) ·
  [Olsson-Collentine 2023 (preregister + report multiverse)](https://consensus.app/papers/details/dc59b4da186451dfa4d7da5358e4e26e/) ·
  [Götz 2024 (multiverse tutorial)](https://consensus.app/papers/details/c8b2375373955a8cb203d5d209030bf1/)
- **Vibration of Effects (VoE).** Quantify how a result moves across model/adjustment choices, including
  sign flips ("Janus effect"). The natural tool for the anchor problem: report the VoE of the discovered
  axis across an **anchor family**.
  Patel & Ioannidis 2015, *J Clin Epi*
  [Consensus](https://consensus.app/papers/details/ea7feebd47855391ac7ccee3410f61ee/) ·
  Tierney 2021, *PLoS Biol*
  [Consensus](https://consensus.app/papers/details/4868511550705a6d9c0c58e8715b9632/)
- **Informed machine learning.** Prior-misspecification sensitivity is an established, evaluable axis.
  von Rueden 2019 (taxonomy: source / representation / integration)
  [Consensus](https://consensus.app/papers/details/4c5ecab02cf2528a99fca77a179bad89/) ·
  Oneto 2025, *PLOS Comp Biol* (best-practice tips)
  [Consensus](https://consensus.app/papers/details/f15f818586a25d558412d67da844e67a/)

### 2d. Standardized anchor sources already exist
- **MSigDB Hallmark** (Liberzon et al., 2015, *Cell Systems*) — 50 curated, de-redundified, coherent gene
  sets distilled from many "founder" sets. A community-standard library to draw anchors from **by rule**
  (proliferation → `HALLMARK_E2F_TARGETS` / `HALLMARK_G2M_CHECKPOINT`, etc.).
  [Consensus](https://consensus.app/papers/details/f377c0d90c665d7daed4eec58d60df5f/)

---

## 3. Key reframe

Our method sits on the **solution** side of this debate: knowledge-anchored residual discovery is a
generalization of Venet's "adjust for the dominant prior, then look at the residual." We should (a) position
it that way explicitly, and (b) adopt the field's standardization devices rather than invent our own.

The honest reframing of every claim: not "*the* textbook anchor," but "**given declared prior X**, here is
the X-orthogonal residual" — conditional, transparent, and testable.

---

## 4. A standardization recipe (what the literature points to)

1. **Anchor provenance (rule, not hand-pick).** Draw the anchor from a standard library (MSigDB Hallmark) or
   a reproducible data recipe (meta-PCNA style: top-k genes correlated with a canonical marker in a
   normal-tissue compendium). Record the exact gene list + source.
2. **Pre-specification.** Declare the anchor *before* inspecting the discovery outcome (multiverse /
   preregistration ethos).
3. **Anchor-family robustness / VoE.** Run discovery across a *family* of plausible standard anchors; accept
   only axes stable across the family (and stable in sign). Per Del Giudice, the family must be
   **biologically equivalent candidates**, not "any gene set."
4. **Replication as the ultimate guard.** Cross-cohort / cross-cancer replication (already achieved for the
   basal axis, PGR, IO biomarkers) is the strongest control against a cherry-picked anchor.

**Net:** prediction is anchor-agnostic and needs no standardization; **discovery** is standardized by
*standard source + anchor-family VoE + conditional framing + replication* — and Venet's meta-PCNA is a
concrete, citable recipe we can adopt directly.

---

## 5. Caveats to keep honest
- **Don't over-multiverse** a principled choice (Del Giudice): a real axis can be hidden among poorly
  justified alternatives. Keep the anchor family small and biologically motivated.
- **Proliferation is sticky** (Venet): excluding cell-cycle genes ≠ removing proliferation; only a
  regressed-out meta-gene works — which our gate already does.
- **VoE is descriptive, not a license:** large vibration → caution; small vibration → robustness. It informs,
  it doesn't "prove."

---

## 6. Decision — options to choose from (for our next step)

| Option | What it is | Effort | What it buys |
|---|---|---|---|
| **A. meta-PCNA recipe helper** | A function that builds a Venet-style data-driven anchor (top-k genes correlated with a canonical marker in a normal-tissue compendium) — standardized, reproducible | Small–Med | Removes hand-curation for the *definition* of an anchor; directly citable to Venet 2011 |
| **B. Anchor-family VoE analysis** | Run the basal discovery across several standard proliferation anchors (MSigDB hallmarks + meta-PCNA) and quantify robustness / sign-stability of the discovered axis | Med | Empirically shows the discovery isn't an artifact of one hand-picked anchor — the core remedy |
| **C. Both (A feeds B)** | Build the recipe, then use it as one member of the anchor family in the VoE check | Med–Large | The complete standardized story |
| **D. Document-only** | Adopt the *reframing + provenance + pre-specification* conventions in the methods note; no new code yet | Tiny | Honest framing now; defer empirical robustness |

My lean: **B (or C)** is the highest-value because it *empirically* answers the cherry-picking worry on our
own flagship result; **A** is the clean way to make anchors reproducible; **D** is worth doing regardless as
the framing layer. None of these touches the (already rigorous) predictive guarantee.

---

## 7. Update — Option C implemented (2026-06)

Both pieces are now built and run:

- **A. meta-PCNA recipe** → `omniomics.multiomics.marker_correlated_anchor(expr, marker="PCNA", top_frac=0.01)`:
  a reproducible, hand-curation-free anchor definition (top genes correlated with a canonical marker), unit-
  tested.
- **B. anchor-family VoE** → `reports/dmoi_anchor_family_voe.py`: the LumA/B discovery run across a family of
  proliferation priors {curated 20-gene index, MSigDB `HALLMARK_E2F_TARGETS` / `G2M_CHECKPOINT` /
  `MYC_TARGETS_V1`, meta-PCNA}.

**Result — the basal axis is robust to anchor choice** (`anchor_family_voe.csv`, `anchor_family_consensus.csv`):

| Anchor | anchor AUROC | overlap w/ basal panel (/30) | basal core (/7) |
|---|---|---|---|
| curated_prolif | 0.919 | 26 | 7 |
| HALLMARK_E2F_TARGETS | 0.914 | 24 | 7 |
| HALLMARK_G2M_CHECKPOINT | 0.894 | 23 | 7 |
| HALLMARK_MYC_TARGETS_V1 | 0.802 | 21 | 7 |
| meta_PCNA (data-driven) | 0.926 | 26 | 7 |

Mean pairwise Jaccard of the discovered panels = **0.69**; **24 consensus genes** recovered by ≥4/5 anchors
(21 in the basal panel). Every anchor — including the weakest (MYC, AUROC 0.80) and the data-driven
meta-PCNA — re-recovers the basal core 7/7. This empirically answers the cherry-picking concern: the
discovery is not an artifact of the hand-picked anchor. Guard: `test_anchor_family_voe_guard`.

## References (verified)
1. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol* 2011;7(10):e1002240. [doi:10.1371/journal.pcbi.1002240](https://doi.org/10.1371/journal.pcbi.1002240) *(PubMed)*
2. Michiels S, Koscielny S, Hill C. Prediction of cancer outcome with microarrays: a multiple random validation strategy. *Lancet* 2005. [Consensus](https://consensus.app/papers/details/5c64938c77ab54c1a9b313a773cfe511/)
3. Nagy Á, et al. Pancancer survival analysis of cancer hallmark genes. *Sci Rep* 2020. [Consensus](https://consensus.app/papers/details/52cbbe35de8e5267aedd830e1e4ed6da/)
4. Del Giudice M, Gangestad SW. A Traveler's Guide to the Multiverse… *Adv Methods Pract Psychol Sci* 2021. [Consensus](https://consensus.app/papers/details/42c0716ece6d552986802d64c1bd1608/)
5. Olsson-Collentine A, et al. Meta-analyzing the multiverse… *Psychol Methods* 2023. [Consensus](https://consensus.app/papers/details/dc59b4da186451dfa4d7da5358e4e26e/)
6. Götz M, et al. The multiverse of universes: a tutorial… *Int J Psychol* 2024. [Consensus](https://consensus.app/papers/details/c8b2375373955a8cb203d5d209030bf1/)
7. Patel CJ, Burford B, Ioannidis JPA. Assessment of vibration of effects due to model specification… *J Clin Epidemiol* 2015. [Consensus](https://consensus.app/papers/details/ea7feebd47855391ac7ccee3410f61ee/)
8. Tierney B, et al. Leveraging vibration of effects analysis for robust discovery… *PLoS Biol* 2021. [Consensus](https://consensus.app/papers/details/4868511550705a6d9c0c58e8715b9632/)
9. von Rueden L, et al. Informed Machine Learning — a taxonomy and survey… *IEEE TKDE* 2019. [Consensus](https://consensus.app/papers/details/4c5ecab02cf2528a99fca77a179bad89/)
10. Oneto L, et al. Eight quick tips for biologically and medically informed machine learning. *PLOS Comput Biol* 2025. [Consensus](https://consensus.app/papers/details/f15f818586a25d558412d67da844e67a/)
11. Liberzon A, et al. The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Systems* 2015. [Consensus](https://consensus.app/papers/details/f377c0d90c665d7daed4eec58d60df5f/)
