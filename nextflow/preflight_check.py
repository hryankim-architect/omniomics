#!/usr/bin/env python3
"""Preflight consistency check for the three GSE57577 nf-core arms.

Validates — without touching the network or running anything — that each arm's ids.csv, sample map,
and (for chipseq) antibody/control wiring are internally consistent, so a multi-hour run does not
fail on a trivial typo. Returns a list of (check, ok, detail); exits non-zero on any failure.
"""
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRR = lambda s: s.startswith("SRR") and s[3:].isdigit()


def _read_ids(path):
    return [l.strip() for l in open(path) if l.strip()]


def check_arm(name, ids_csv, map_csv, chipseq=False):
    out = []
    ids = _read_ids(ids_csv)
    rows = list(csv.DictReader(open(map_csv)))
    map_runs = [r["run"] for r in rows]

    out.append((f"{name}: all ids look like SRR accessions", all(SRR(x) for x in ids),
                [x for x in ids if not SRR(x)]))
    out.append((f"{name}: ids.csv has no duplicates", len(ids) == len(set(ids)),
                sorted({x for x in ids if ids.count(x) > 1})))
    out.append((f"{name}: map runs unique", len(map_runs) == len(set(map_runs)),
                sorted({x for x in map_runs if map_runs.count(x) > 1})))
    out.append((f"{name}: ids set == map runs set", set(ids) == set(map_runs),
                {"ids_only": sorted(set(ids) - set(map_runs)), "map_only": sorted(set(map_runs) - set(ids))}))

    if chipseq:
        names = {r["sample"] for r in rows}
        controls = {r["control"] for r in rows if r.get("control")}
        out.append((f"{name}: every control references an existing sample", controls <= names,
                    sorted(controls - names)))
        chips = [r for r in rows if r.get("antibody")]
        inputs = [r for r in rows if not r.get("antibody")]
        out.append((f"{name}: has >=1 input (blank antibody)", len(inputs) >= 1, len(inputs)))
        out.append((f"{name}: every ChIP has a control + control_replicate",
                    all(r.get("control") and r.get("control_replicate") for r in chips),
                    [r["sample"] for r in chips if not (r.get("control") and r.get("control_replicate"))]))
    return out


def run():
    checks = []
    checks += check_arm("rnaseq", os.path.join(HERE, "ids.csv"), os.path.join(HERE, "srr_sample_map.csv"))
    checks += check_arm("methylseq", os.path.join(HERE, "methylseq", "ids.csv"),
                        os.path.join(HERE, "methylseq", "srr_sample_map.csv"))
    checks += check_arm("chipseq", os.path.join(HERE, "chipseq", "ids.csv"),
                        os.path.join(HERE, "chipseq", "srr_sample_map.csv"), chipseq=True)
    return checks


if __name__ == "__main__":
    checks = run()
    bad = 0
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + ("" if ok else f"   -> {detail}"))
        bad += 0 if ok else 1
    print(f"\n{'ALL PASS ✅' if not bad else f'{bad} FAILURE(S) ❌'}")
    sys.exit(1 if bad else 0)
