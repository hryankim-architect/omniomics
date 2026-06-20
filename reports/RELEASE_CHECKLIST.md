# Release & submission checklist

Status as of this revision: package **omniomics 0.3.0** builds cleanly (`twine check` PASSED, clean-venv
import verified); manuscript is PubMed-referenced with real ORCID and declarations. Upload and submission are
**your actions** (they need credentials / a logged-in browser); everything below is prepared.

## A. PyPI release (omniomics 0.3.0)

Build artifacts (regenerate locally; `dist/` is gitignored):

```bash
cd ~/Downloads/molcell_epigenetics/omniomics-prototype
python -m build                 # -> dist/omniomics-0.3.0-{tar.gz,-py3-none-any.whl}
python -m twine check dist/*     # both should PASS
```

Before uploading:
- **Check the name is free** on PyPI: https://pypi.org/project/omniomics/ . If taken, rename in
  `pyproject.toml` (e.g. `omniomics-anchored`) and rebuild.
- Test on **TestPyPI** first: `twine upload --repository testpypi dist/*` then
  `pip install -i https://test.pypi.org/simple/ omniomics`.

Upload (your action — needs a PyPI API token):
```bash
python -m twine upload dist/*
```

Version history: 0.1.0 (engine) → 0.2.0 (anchored integration + knowledge-anchored discovery) → **0.3.0**
(out-of-core `omniomics.scale`, gated learned-representation adoption `omniomics.representations`, plus
scale-prep in `anchored_residual_discovery`: vectorized partial-corr, SIS `screen_top`, `n_jobs`, BH-FDR
`n_fdr`, `leave_one_cohort_out`). Optional extras: `[scale]` (pyarrow, h5py), `[deep]` (scvi-tools),
`[enrichment]`, `[modern]`, `[test]`.

## B. bioRxiv preprint submission

Manuscript: `reports/anchored_integration_manuscript.pdf` (~10 pp, xelatex). It already contains:
- Title, author **H. Ryan Kim** with ORCID 0000-0002-1869-0412, corresponding email
  (**hryankim1221@gmail.com**), subject areas (Bioinformatics; Cancer Biology).
- Structured abstract; Introduction; Results (anchored gate; knowledge anchor; residual discovery; HER2/ER
  generalization; METABRIC/lung/HNSC external validation; NSCLC cross-domain; clinical-significance negative;
  **hypothesis-as-anchor + commonality/mediation; 4-endpoint × 4-cohort transportability panel; platform
  sign-flip**); Discussion; Methods (with a Scalability note); Data & Code Availability; **Declarations**
  (competing interests: none; funding: none; author contributions; ethics: public de-identified data;
  **CC-BY 4.0**; AI-assistance disclosed); **12 references, DOI-linked** (10 PubMed-verified; Tonidandel &
  LeBreton 2011 and Degtiar & Rose 2023 are stats-methodology refs, flagged not-in-PubMed).
- **Four figures** (all present in `reports/figs/`): Figure 1 `discovery_summary.png` (anchored gate / knowledge
  anchor / generalization), Figure 2 `hypothesis_anchor.png` (hypothesis-as-anchor + screen + commonality +
  transportability sweep, panels A–D), Figure 3 `endpoint_panel.png` (4-endpoint × 4-cohort label grid),
  Figure 4 `platform_corr.png` (corr(prolif,ER) sign flip by platform).

Submission steps (your action — at https://www.biorxiv.org/submit-a-manuscript):
1. Confirm corresponding-author ORCID and the CC-BY 4.0 license selection match the manuscript.
2. Upload the PDF (or convert `anchored_integration_manuscript.md` to the journal's preferred format).
3. Upload the four figure files above if separate figure files are requested.
4. Category: **Bioinformatics** (cross-list Cancer Biology). Declare AI assistance per bioRxiv policy.
5. Add the code-availability URL: the `omniomics` GitHub repository.
6. Data: all public — TCGA (BRCA RNA-seq + Agilent, LUAD, LUSC, HNSC) via UCSC Xena; METABRIC + NSCLC anti-PD-1
   via cBioPortal; **SCAN-B via GEO GSE96058**. No data upload needed (SCAN-B inputs regenerate via
   `reports/fetch_scanb.sh`).

The preprint companion `reports/anchored_integration_preprint.pdf` (shorter, narrative) carries the same
results and an aligned abstract.

## C. Repository state

Clean working tree; all results committed with reproduce-runners + skip-safe CI guards
(suite: **78 passed, 1 skipped**). Large public data live outside the repo (referenced by env vars
BRCA_DIR/METABRIC_DIR/SCANB_DIR), matching TCGA/METABRIC. Per-round node-handoff scripts (`commit-*.sh`) are
gitignored scaffolding, not deliverables.
