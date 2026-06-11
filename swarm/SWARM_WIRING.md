# Wiring omniomics into the Agentic Bioinformatics Swarm

Drop-in bundle that registers the omniomics engine + golden tasks into `swarm-starter`'s Claude Code
layer, alongside the existing RAG `/eval` loop. Realizes the roadmap's Phase 5 (evaluation +
provenance) inside the swarm.

## Install (on the swarm host)

```bash
# 1. place the engine next to the swarm
cp -r omniomics-prototype ~/omniomics

# 2. add the slash-command + agent to the swarm's Claude config
cp omniomics-prototype/swarm/commands/golden.md   ~/swarm/.claude/commands/
cp omniomics-prototype/swarm/agents/omniomics-runner.md ~/swarm/.claude/agents/

# 3. merge the permission + audit-hook additions into ~/swarm/.claude/settings.json
#    (see settings.snippet.json — add the Bash allows and keep the existing PostToolUse audit hook)
```

## How it plugs into the existing swarm

| Swarm component (swarm-starter) | omniomics use |
|---------------------------------|----------------|
| `/eval` slash-command | parallel `/golden` — known-answer regression for the analysis engine |
| `swarm-engineer` / `swarm-reviewer` agents | add `omniomics-runner`; reviewer still gates merges |
| `post-edit-audit.sh` (hash-chained `audit.log`) | golden runs append their own hash-chained `golden/audit.log` (same PROV-AGENT pattern) |
| `critic.py` / `approval.py` | `on_fail: block_merge` — a golden FAIL blocks, HITL approves overrides |
| `bixbench_runner.py` / `genotex_runner.py` | `golden_tasks.yaml` is the same idea: curated known-answer tasks; add `gse57577_reproduction` + `brca_multiomics_pipeline` to the benchmark rotation |
| BioMCP (`biomcp_*`) | expose GEO/SRA fetch so `omniomics.geo` ingestion is a tool call |

## Daily loop (post-install)

```
cd ~/swarm && claude
> @omniomics-runner add GSE_xxxxx to the cohort and harmonize
> /golden            # regression gate — must be green
> @swarm-reviewer review the last diff
```

Every Edit/Write fires the audit hook; every `/golden` appends a provenance record. As the engine
scales (more datasets, scVI/EB-ComBat swaps, new assays), the golden suite catches silent regressions
and the chain keeps each run reproducible and auditable.

## Governance contract

- A golden FAIL **blocks merge** (`critic.py`); only a human-approved override (`approval.py`) proceeds.
- Tolerances live in `golden/golden_tasks.yaml` and are **not** loosened to make a run pass.
- `omniomics-runner` must run `/golden` and verify the audit chain before declaring work done.
