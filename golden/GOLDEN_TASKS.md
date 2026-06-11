# Golden tasks — provenance-gated regression for the omniomics engine

Known-answer reproductions the swarm must keep passing as it evolves. This is the concrete
realization of the roadmap's Phase 5 (evaluation + provenance) and mirrors the
`bixbench_runner` / `genotex_runner` pattern already in `swarm-starter/rag/`.

## Tasks (`golden_tasks.yaml`)

| Task | Command | What it asserts |
|------|---------|-----------------|
| `gse57577_reproduction` | `python run_golden.py` | WWD=1888 / R=9 / TKO=3 DE genes, 6/6 named targets, ChIP 1.64× |
| `brca_multiomics_pipeline` | `python golden/run_golden_brca.py` | cross-cohort batch removal (PC1 0.99→0), transfer AUROC 0.75→0.87, promoter-methylation 74% negative, RNA subtype 0.91, compact pole 0.91 vs 2000-gene 0.95 |

Both currently **PASS** unattended.

## Provenance (hash-chained audit)

Every golden run appends a record to `audit.log`, each carrying the SHA-256 of
`prev_hash + record` — a tamper-evident chain (PROV-AGENT-style), identical in spirit to the
swarm's `post-edit-audit.sh`. Verified:

```
chain valid: True
```

So any scaled analysis carries a verifiable trail: which datasets, which metrics, which result,
linked to the previous run. Provenance becomes a quality signal, not just compliance.

## Wiring into the swarm

```
/golden                      # slash-command alongside /eval (run all golden tasks)
  → engine runs pipelines unattended
  → metrics asserted vs golden_tasks.yaml tolerances
  → swarm-reviewer agent must confirm before merge
  → critic.py / approval.py BLOCK on any FAIL (on_fail: block_merge)
  → hash-chained record appended to golden/audit.log
```

This closes the loop the roadmap describes: as the swarm adds datasets, swaps ComBat-lite for
scVI/EB-ComBat, or changes the integrator, these golden tasks catch silent regressions, and every
run is reproducible and audited.

## Run
```bash
python golden/run_golden_brca.py             # re-runs the full pipeline, then asserts
python golden/run_golden_brca.py --use-cached # assert against existing metric CSVs
```
