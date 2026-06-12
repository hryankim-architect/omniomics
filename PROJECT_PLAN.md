# Project Plan — *Re-reading an epigenetics classic with DMOI*

**A modern, reproducible reanalysis of Noh et al. (Molecular Cell 2015, GSE57577) that re-runs the
three assays through modern pipelines + statistics, adds modern mechanistic/foundation-model layers,
and demonstrates DMOI (Dialectical Multi-Omics Integration) reconstructing the paper's central
pattern — all packaged as the installable, golden-tested `omniomics` engine.**

---

## 1. Motivation

Noh et al. engineered the Dnmt3a ADD histone-reading domain (WT / WWD / R) in Dnmt-TKO mouse ESCs and
read out the consequences with bulk RNA-seq (Cufflinks FPKM), RRBS (methylKit), and ChIP-seq density —
hand-integrating them into the now-iconic graphical-abstract matrix (ADD variant × chromatin context →
localization, methylation, phenotype). The biology is excellent; the *analysis* is 2015-era. Eleven
years on, we can (a) re-run it with reproducible Nextflow pipelines + replicate-based statistics,
(b) add sequence-to-epigenome and foundation-model layers, and (c) show that a single interpretable
integration method — **DMOI** — reconstructs the paper's central pattern from the data, faithfully
marking what the data can and cannot support. The output is both a clean reproducible-research artifact
and a concrete demonstration of DMOI's value (and limits).

## 2. Specific aims

**Aim 1 — Modern reproduction.** Re-derive the paper's three-assay results with current best practice:
nf-core pipelines from raw reads, count-based differential expression, proper DMR calling,
spike-in-aware ChIP quantification — versus the original processed files — and quantify where modern
methods change the conclusions.

**Aim 2 — DMOI reconstruction of the graphical abstract (core demo).** Formalize and extend the
already-working DMOI reconstruction: the Dnmt3a localization matrix, the binding↔methylation
*disagreement* layer (WWD binds more but methylates less), and the ESC-differentiation phenotype —
with rigorous statistics and an honest "not assayed" mark for the H3T3ph/centromere axis.

**Aim 3 — Modern mechanistic + foundation-model layer.** Predict the WWD/R redistribution and its
methylation output from *sequence* (Enformer/Borzoi/EPInformer), embed the system in
foundation-model space, and address the missing centromere axis *in silico* / via a modern-assay design.

**Aim 4 — Reproducible, installable, provenance-gated delivery.** Ship it as `omniomics`
(`pip install -e`, golden-task pytest CI, hash-chained audit, swarm drop-in) — already largely built.

## 3. What already exists (leverage — do not rebuild)

| Asset | Status |
|-------|--------|
| GSE57577 reproduction (RNA/ChIP/RRBS, moderated DE, CGI, GO, aneuploidy) | ✅ done, golden-tested |
| DMOI graphical-abstract reconstruction (`run_gse57577_dmoi.py`) | ✅ localization + disagreement + phenotype |
| `omniomics` engine (loaders, EB-ComBat, MOFA-lite, DMOI representation) | ✅ pip-installable, CI green |
| Cross-cohort / multi-omics scaling (TCGA-METABRIC, methylation arm) | ✅ done |
| Provenance: golden_tasks.yaml + hash-chained audit + swarm wiring | ✅ done |

The project is therefore **~50% complete**; the plan focuses the remaining effort on Aims 1 and 3 and
on polishing Aim 2 into a publishable artifact.

## 4. Technical approach

**Aim 1 — modern reproduction**
- Pull raw FASTQ from SRA (SRP041894/95) via the existing GEO/SRA tooling.
- RNA-seq → **nf-core/rnaseq** (STAR + Salmon) → **DESeq2/edgeR on counts** (replicate-based, vs the
  paper's FPKM eyeballing and our paired limma-moderated test at n=2).
- Methylation → **nf-core/methylseq** (Bismark) → **DSS/methylKit DMRs** with proper test statistics;
  consider **EM-seq**-style coverage modeling.
- ChIP-seq → **nf-core/chipseq** (or CUT&Tag reanalysis) with **spike-in normalization** so
  redistribution is a quantitative, not relative, claim.
- Deliverable: a "2015 vs 2026" concordance table — which findings hold, which sharpen, which were
  artifacts (e.g., the WWD chrY-loss clonal confound we already found).

**Aim 2 — DMOI core demo (polish + rigor)**
- Generalize `omniomics.multiomics.dmoi_representation` to arbitrary context poles; compute the
  variant × context localization matrix and the per-context binding↔methylation disagreement with
  bootstrap CIs.
- Tie phenotype to RNA (differentiation program) with a proper signature score + CI.
- Produce a single publication-quality figure that *is* the modernized graphical abstract, plus a
  methods note framing DMOI as interpretable structured fusion (vs black-box concatenation).

**Aim 3 — mechanistic + foundation layer**
- **Sequence→epigenome**: run Enformer/Borzoi on the H3K4me3-CGI / gene-body / enhancer loci to test
  whether the *sequence* predicts where Dnmt3a redistributes and where methylation is suppressed;
  EPInformer for expression. Frame as "does the engineered readout act through predictable cis-context."
- **Foundation embeddings**: place the WT/WWD/R ESC states in scGPT/Geneformer space (or a methylation
  FM) to see if the perturbation has a coherent latent direction.
- **The missing centromere axis**: the H3T3ph/mitotic-centromere column needs mitotic CUT&Tag +
  nanopore satellite methylation — out of scope to generate, so deliver (i) a precise experimental
  design (already in `PROTOCOL_design.md`) and (ii) an *in-silico* placeholder analysis of
  pericentromeric/satellite methylation from any available long-read data, clearly labeled as a gap.

**Aim 4 — delivery** (mostly done): keep golden CI green; add a reproducible `Makefile`/`nextflow`
entry that runs the whole reanalysis end-to-end; publish the repo + a short write-up.

## 5. Deliverables

1. **Modernized figure** — the DMOI-reconstructed graphical abstract (localization + disagreement +
   phenotype), publication-quality.
2. **"2015 vs 2026" reproduction report** — concordance table + what modern methods change.
3. **`omniomics` v1.0** — tagged release, `pip install`, golden CI badge, DOI via Zenodo.
4. **Short paper / preprint (bioRxiv) or technical blog post** — "An interpretable multi-omics
   re-reading of an ADD-engineering classic," honest about where multi-omics/DMOI helps and where it
   doesn't.
5. **Swarm golden task** — the reanalysis registered as a known-answer regression in the swarm.
6. **Talk / portfolio piece** — 10-slide deck of the arc.

## 6. Milestones & timeline (8 weeks, ~part-time)

| Week | Milestone |
|------|-----------|
| 1 | Pull SRA raw data; stand up nf-core/rnaseq + methylseq on the swarm; smoke run |
| 2 | Aim 1 RNA-seq: DESeq2 vs FPKM concordance; lock the count pipeline |
| 3 | Aim 1 methylation + ChIP: DMRs + spike-in ChIP; "2015 vs 2026" table v1 |
| 4 | Aim 2: DMOI representation generalized + bootstrap CIs; modernized figure v1 |
| 5 | Aim 3: Enformer/Borzoi on the loci; sequence-predicts-redistribution test |
| 6 | Aim 3: foundation-model embedding + centromere-gap design/in-silico note |
| 7 | Write-up draft (report + preprint/blog); omniomics v1.0 release + Zenodo DOI |
| 8 | Polish figures, deck, golden task; submit preprint / publish post |

(Buffer assumed; Aims 1 and 3 are the long poles. Aim 2/4 are mostly polish.)

## 7. Data & compute

- **Public**: GSE57577 (RNA/ChIP/RRBS), SRA SRP041894/95 raw reads; mm9/mm10 + UCSC annotation;
  optional public scNMT-seq / nanopore mESC data for the centromere axis.
- **Compute**: the existing 5-node swarm (chi-mac-p brain + nf-core on a Linux node). Enformer/Borzoi
  inference is GPU-friendly but runnable on CPU for a handful of loci.
- **No private data needed** for the core (unlike the BRCA arm), so the whole project is openly
  reproducible — a deliberate contrast to the TCGA work.

## 8. Success criteria

- Aim 1: ≥ the paper's headline findings reproduced from raw reads; a documented concordance table.
- Aim 2: the modernized figure reconstructs WWD→CGI localization (≈1.6×), the disagreement layer, and
  the differentiation phenotype, with CIs — and the golden test asserts the key numbers.
- Aim 3: a quantified statement on whether sequence models predict the redistribution (positive *or*
  negative is publishable if rigorous).
- Aim 4: `pip install omniomics` + green CI + one-command end-to-end rerun; preprint/blog out.

## 9. Risks & mitigations (honest)

| Risk | Mitigation |
|------|-----------|
| Raw-read pipelines are heavy / slow | Run on the Linux swarm node; cache; subset chromosomes for dev |
| H3T3ph/centromere axis not in the data | Deliver design + in-silico note, label the gap (don't fabricate) |
| Sequence models may *not* predict the redistribution | A clean negative is a valid, honest result |
| Multi-omics/DMOI gains are modest (we showed this) | Frame DMOI as *interpretable reconstruction*, not an accuracy contest |
| Scope creep across 4 aims | Aims 2 & 4 are ~done; protect time for Aim 1 (the long pole) |

## 10. Stretch goals

- scNMT-seq-style single-cell reanalysis if public mESC data exist (true joint methylation+transcriptome).
- A small "engineer-the-reader" in-silico screen: perturb ADD specificity in a sequence model and
  predict the resulting redistribution.
- Generalize the DMOI graphical-abstract reconstruction into a reusable `omniomics` report template for
  *any* perturbation-with-localization study.

## 11. Key references (modern methods)

- nf-core pipelines (rnaseq / methylseq / chipseq) — https://nf-co.re
- DESeq2 — https://bioconductor.org/packages/DESeq2 ; DSS DMRs — https://bioconductor.org/packages/DSS
- Enformer (Avsec 2021) / Borzoi — https://www.nature.com/articles/s41592-021-01252-x ; https://github.com/calico/borzoi
- EPInformer — https://pmc.ncbi.nlm.nih.gov/articles/PMC11312614/
- scGPT — https://www.nature.com/articles/s41592-024-02201-0 ; MOFA+ / scvi-tools (MultiVI)
- scNMT-seq (Clark 2018) — joint single-cell methylation + transcriptome + accessibility
- Companion docs in this repo: `PROTOCOL_design.md` (assay design incl. the centromere axis),
  `SCALING_RESEARCH_ROADMAP.md`, `PROTOTYPE_REPORT.md`.
