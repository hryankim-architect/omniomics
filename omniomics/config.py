"""Path configuration so scripts run on any machine (sandbox, chi-mac-p, …).

Resolution order:
  - dmoi-brca-poc data : env DMOI_BRCA_DATA, else auto-detect common locations
  - methylation cache  : env OMNIOMICS_CACHE, else <repo>/data/brca_cache
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)               # omniomics-prototype/

def _first_existing(paths):
    for p in paths:
        if p and os.path.isdir(p):
            return p
    return None

def brca_data_dir():
    """Path to dmoi-brca-poc/data (contains tcga_brca/, metabric/, msigdb/)."""
    cand = _first_existing([
        os.environ.get("DMOI_BRCA_DATA"),
        os.path.expanduser("~/Downloads/AI/dmoi-brca-poc/data"),
        os.path.expanduser("~/AI/dmoi-brca-poc/data"),
        os.path.expanduser("~/dmoi-brca-poc/data"),
        "/sessions/sleepy-blissful-allen/mnt/AI/dmoi-brca-poc/data",  # sandbox
    ])
    if cand is None:
        raise FileNotFoundError(
            "Could not find dmoi-brca-poc/data. Set DMOI_BRCA_DATA, e.g.\n"
            "    export DMOI_BRCA_DATA=~/path/to/dmoi-brca-poc/data")
    return cand

def brca_tcga_dir(): return os.path.join(brca_data_dir(), "tcga_brca")
def brca_metabric_dir(): return os.path.join(brca_data_dir(), "metabric")
def hallmark_gmt(): return os.path.join(brca_data_dir(), "msigdb", "h.all.v2024.1.Hs.symbols.gmt")

def cache_dir():
    """Where prepare_brca.py writes methylation-context artifacts."""
    d = os.environ.get("OMNIOMICS_CACHE", os.path.join(_REPO, "data", "brca_cache"))
    os.makedirs(d, exist_ok=True)
    return d

def repo_dir(): return _REPO
