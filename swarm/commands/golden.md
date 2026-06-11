---
description: Run the omniomics golden-task suite (known-answer reproductions) and report regressions. Blocks merge on any FAIL.
allowed-tools: Bash
---

Run the omniomics golden tasks — known-answer reproductions the swarm must keep passing as it
evolves. Mirrors `/eval` (RAG regression) but for the multi-omics analysis engine.

**Stage 1 — GSE57577 reproduction (fast):**

```bash
cd ~/omniomics && python run_golden.py
```

Expected: `ALL PASS` (WWD=1888 / R=9 / TKO=3 DE genes, 6/6 named targets, ChIP 1.64×).
If any check FAILs, **stop and report** — a DE or ChIP regression must be diagnosed before the
multi-omics pipeline, which depends on the same engine.

**Stage 2 — BRCA multi-omics pipeline (slow, re-runs cross-cohort + methylation + joint embedding):**

```bash
cd ~/omniomics && python golden/run_golden_brca.py --use-cached   # assert against current outputs
# or: python golden/run_golden_brca.py                            # full unattended re-run
```

Expected: all 13 checks PASS — cross-cohort batch removal (PC1 0.99→0), transfer AUROC 0.75→0.87,
promoter methylation 74% negative, RNA subtype 0.91, compact pole 0.91 vs 2000-gene 0.95,
enhancer-DMOI ≥0.91 and > promoter > RNA-pole, EB-ComBat+covariate preserves biology (≥0.85) while
naive ComBat-lite destroys it (<0.6).

Each run appends a SHA-256 hash-chained record to `golden/audit.log`. Verify chain integrity:

```bash
cd ~/omniomics && python - <<'PY'
import json,hashlib; prev="0"*64; ok=True
for l in open("golden/audit.log"):
    if not l.strip(): continue
    r=json.loads(l); h=r.pop("hash"); p=r.pop("prev")
    ok &= (p==prev and hashlib.sha256((p+json.dumps(r,sort_keys=True)).encode()).hexdigest()==h); prev=h
print("audit chain valid:", ok)
PY
```

## Report format

```
GSE57577 golden:   PASS/FAIL  (DE 1888/9/3, targets 6/6, ChIP 1.64×)
BRCA pipeline:     X/13 PASS
audit chain:       valid / BROKEN

[if any FAIL]
Failed: <check name> = <value> (golden <expected±tol>)
```

If all pass: **"✅ no regression — golden tasks green, provenance intact."**
On any FAIL, surface the metric delta and **block the merge** (hand to `@swarm-reviewer`).
