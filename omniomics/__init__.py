"""omniomics — a manifest-driven, reusable multi-omics analysis engine.

Phase 0/1 prototype distilled from the GSE57577 reproduction. Turns the one-off
analysis into reusable components that run on any dataset described by a manifest,
and harmonizes many studies onto a common feature space for cross-study analysis.
"""
__all__ = ["geo", "loaders", "expression", "harmonize", "multiomics"]
__version__ = "0.2.0"
