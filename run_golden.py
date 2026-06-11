#!/usr/bin/env python3
"""GSE57577 golden task (thin wrapper; logic lives in omniomics.golden_check).
Run:  python run_golden.py   or   omniomics-golden"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from omniomics.golden_check import run

if __name__ == "__main__":
    ok, _, _ = run()
    sys.exit(0 if ok else 1)
