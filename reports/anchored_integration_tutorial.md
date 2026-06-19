# Anchored multi-omics integration — quickstart

A small, interpretable, **never-below-the-best-single-modality** integrator. Anchor on the most
predictive/robust modality, then add the others *only where they earn it*. Full method write-up:
[`DMOI_method_assessment.md`](DMOI_method_assessment.md). API lives in `omniomics/multiomics.py`.

## Why

On real cancer endpoints a strong single modality (usually RNA) is hard to beat, and symmetric fusion
often *loses* to it (a weak layer drags a strong one down). Anchored integration fixes that by design:
it returns the best single modality whenever the others add nothing, and improves on it only when a
modality carries orthogonal signal. Grounded in cooperative learning (Ding & Tibshirani, PNAS 2022)
and late-fusion benchmarks (the dominant modality is task-specific; weight by individual CV success).

## One call

```python
import numpy as np
from omniomics import multiomics as mo

# modalities: a dict of {name: X}, each X is (samples x features); y is binary labels (aligned rows)
res = mo.auto_integrate({"RNA": Xrna, "methylation": Xmeth, "CNV": Xcnv}, y, cv=5)

res["anchor"]          # modality chosen as the backbone (data-driven, per task)
res["auroc_anchor"]    # anchor-only out-of-fold AUROC (the floor)
res["auroc_combined"]  # after gated forward integration (>= floor)
res["delta"]           # combined - anchor  (>= ~0; > 0 only where a modality earns it)
res["added"]           # {modality: (n_folds_added, mean_beta)} -- which modalities entered
res["ranking"]         # [(name, mean_auroc, std, composite), ...] from select_anchor
```

## The pieces (use directly if you want control)

```python
mo.select_anchor({"RNA": Xrna, "methylation": Xmeth}, y)          # rank modalities; pick the anchor
mo.anchored_integrate(Xanchor, Xsecondary, y)                    # 2-modality gated residual
mo.forward_integrate({...}, y)                                   # N-modality greedy forward gating
mo.anchored_gate(anchor_prob, secondary_score, y, margin=0.01)   # the beta>=0 gate primitive
```

Key knobs: `gate_margin` (a secondary must beat the anchor by more than this to be added — guards
small-n selection noise), `inner_repeats` (averages the inner-CV gain estimate), `cv`.

## Reproduce the TCGA-BRCA demo

```bash
BRCA_DIR=/path/to/tcga_brca python run_auto_integrate.py    # writes auto_integrate_results.csv
```

| endpoint | n | anchor | anchor AUROC | combined | added |
|---|---|---|---|---|---|
| LumA/B (RNA-defined) | 417 | RNA | 0.940 | 0.940 | none |
| methylation-defined (positive control) | 417 | methylation | 0.980 | 0.980 | none |

`auto_integrate` picks the dominant modality per task and never falls below it; on a synthetic case
where a second modality carries orthogonal signal it engages and gains (see the method write-up §8).

## External validation — real subtype labels (not RNA-biased)

Run blind on nine expert-defined TCGA-BRCA labels (n=491), `auto_integrate` routes to the modality
that defines each task — proof the anchor selection follows the biology, not a built-in RNA preference:

```bash
for ep in methC1 methC2 methC3 methC4 methC5 PAM50Basal PAM50LumA PAM50LumB PAM50Her2; do
  python reports/dmoi_external_subtype.py $ep   # appends to external_subtype_results.csv
done
```

| label family | endpoints | anchor chosen | leader vs other |
|---|---|---|---|
| methylation clusters 1–5 (methylation-defined) | 5/5 | **methylation** | meth 0.70–0.97 vs RNA 0.63–0.91 |
| PAM50 Basal/LumA/LumB/Her2 (expression-defined) | 4/4 | **RNA** | RNA 0.92–0.99 vs meth 0.82–0.99 |

Methylation wins and is selected on every methylation-defined cluster; RNA wins and is selected on
every PAM50 call; the gate never falls below the leader on any of the nine (guarded in CI).

## Guarantees & honest limits

- **Never below the anchor** (up to small-sample noise, which `gate_margin` + `inner_repeats` control).
- **Adds a modality only where it earns a margin gain** — interpretable via `res["added"]`.
- It **cannot manufacture signal the assay did not measure**: across four real BRCA endpoints (incl.
  the epigenetic-clock age task) RNA dominated and nothing was added — genuine multi-omics gains are
  rarer than commonly claimed. Use it to *avoid false multi-omics claims* as much as to capture real ones.

Tested by `tests/test_dmoi_v2.py` (unit) and `tests/test_golden.py` (recorded-metric CI guards).
