"""Portable path resolution for the v2 reproducibility bundle.

Set RV_REPO_ROOT to override the repository root if the bundle is relocated.
"""
import os
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
BUNDLE = SCRIPTS.parent                    # <repo>/repro/v2
REPO = Path(os.environ.get("RV_REPO_ROOT", BUNDLE.parents[1]))
SCI = BUNDLE / "results" / "science"
AUDIT = BUNDLE / "results" / "audit"
FIGS = BUNDLE / "figures"
PROV = BUNDLE / "provenance"
DATA = REPO / "data"
