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

Manuscript: `reports/anchored_integration_manuscript.pdf` (5 pp, xelatex). It already contains:
- Title, author **H. Ryan Kim** with ORCID 0000-0002-1869-0412, corresponding email, subject areas
  (Bioinformatics; Cancer Biology).
- Structured abstract; Introduction; Results (anchored gate; knowledge anchor; residual discovery; HER2/ER
  generalization; METABRIC/lung/HNSC external validation; NSCLC cross-domain; clinical-significance negative);
  Discussion; Methods (with a Scalability note); Data & Code Availability; **Declarations** (competing
  interests: none; funding: none; author contributions; ethics: public de-identified data; **CC-BY 4.0**;
  AI-assistance disclosed); **7 PubMed-verified, DOI-linked references**; Figure 1 (3 panels).

Submission steps (your action — at https://www.biorxiv.org/submit-a-manuscript):
1. Confirm corresponding-author ORCID and the CC-BY 4.0 license selection match the manuscript.
2. Upload the PDF (or convert `anchored_integration_manuscript.md` to the journal's preferred format).
3. Upload Figure 1 (`reports/figs/discovery_summary.png`) if a separate figure file is requested.
4. Category: **Bioinformatics** (cross-list Cancer Biology). Declare AI assistance per bioRxiv policy.
5. Add the code-availability URL: the `omniomics` GitHub repository.
6. Data: all public (TCGA via UCSC Xena; METABRIC + NSCLC anti-PD-1 via cBioPortal) — no data upload needed.

## C. Repository state

Clean working tree; all results committed with reproduce-runners + skip-safe CI guards (suite: 61 passed,
1 skipped). Per-round node-handoff scripts (`commit-*.sh`) are gitignored scaffolding, not deliverables.
