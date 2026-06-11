"""Console entry points (see pyproject [project.scripts])."""
import sys

def golden():
    """omniomics-golden — run the GSE57577 reproduction golden task; exit non-zero on failure."""
    from .golden_check import run
    ok, _, _ = run()
    sys.exit(0 if ok else 1)

def prepare_brca():
    from .prepare import main
    main()
