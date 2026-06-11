# omniomics — local / chi-mac-p setup

`omniomics` is now a `pip install -e` package; all scripts resolve data paths via `omniomics.config`
(env vars), so they run unchanged on the sandbox, your Mac, or the swarm host.

## 1. Install (editable)

```bash
cd ~/Downloads/molcell_epigenetics/omniomics-prototype   # or wherever the repo lives
pip install -e .            # installs omniomics + deps (numpy/pandas/scipy/sklearn/matplotlib)
# optional: pip install -e ".[enrichment]"   # adds gseapy for run_go.py
```
After this, `import omniomics` works from anywhere, and `omniomics-prepare-brca` is on PATH.

## 2. Point at your data

Two env vars (set in `~/.zshrc` to persist):

```bash
export DMOI_BRCA_DATA=~/Downloads/AI/dmoi-brca-poc/data    # contains tcga_brca/ metabric/ msigdb/
export OMNIOMICS_CACHE=~/.omniomics_cache                  # where methylation artifacts go (optional)
```
`config.brca_data_dir()` also auto-detects `~/Downloads/AI/dmoi-brca-poc/data` and a few common
locations, so the env var is only needed if your path differs.

## 3. Build the methylation cache (once)

The methylation-context scripts need probe→context maps + filtered HM450 matrices. Build them once:

```bash
omniomics-prepare-brca          # or:  python -m omniomics.prepare
```
This downloads hg19 refGene + CpG-islands, classifies HM450 probes (promoter / CGI-shore /
distal-enhancer), and writes `pg_*.tsv`, `promoter_meth.tsv.gz`, `enh_meth_all.tsv.gz`,
`ctx_meth_lumab.tsv.gz` into the cache. Idempotent — skips anything already built. (~2–4 min; it
streams the 485k-probe matrix in pure Python, no awk/shell.)

## 4. Run

```bash
# self-contained (downloads its own GEO data):
python run_golden.py                 # GSE57577 reproduction golden task
python run_gse57577_dmoi.py          # graphical-abstract DMOI reconstruction
python run_cohort.py                 # mouse cross-study harmonization
python run_mouse_n3_combat.py        # EB-ComBat on 3 real mouse studies

# uses DMOI_BRCA_DATA (+ cache from step 3):
python run_cohort_brca.py            # TCGA vs METABRIC cross-platform harmonization
python run_meth_arm.py               # RNA + methylation multi-omics
python run_dmoi_enhancer.py          # enhancer-DMOI fusion (repeated-CV significance)
python run_combat_benchmark.py       # EB-ComBat vs naive (confounded N=4)
python golden/run_golden_brca.py --use-cached   # full pipeline regression gate + audit chain
```

## 5. Swarm integration (chi-mac-p)

```bash
cp -r ~/Downloads/molcell_epigenetics/omniomics-prototype ~/omniomics && cd ~/omniomics && pip install -e .
cp swarm/commands/golden.md   ~/swarm/.claude/commands/
cp swarm/agents/omniomics-runner.md ~/swarm/.claude/agents/
# merge swarm/settings.snippet.json into ~/swarm/.claude/settings.json
```
Then inside `claude` on the swarm: `@omniomics-runner add a dataset` → `/golden` (regression gate,
blocks merge on FAIL) → `@swarm-reviewer`. See `swarm/SWARM_WIRING.md`.

## Notes
- Scripts keep a `sys.path.insert(__file__)` shim so they also run from the repo folder *without*
  installing — `pip install -e` just makes `omniomics` importable everywhere (needed by the swarm).
- macOS default `python3` (or conda) works; tested on Python 3.11/3.12.
