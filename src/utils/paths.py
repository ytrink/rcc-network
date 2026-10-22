"""
src/utils/paths.py

Shared path helpers. All inputs live under <repo>/data (see README for the expected layout and
the Zenodo download); all outputs are written under <repo>/results.
"""

from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT_DIR / "data"
RESULTS_DIR = ROOT_DIR / "results"


def data_path(*parts) -> Path:
    """Path to an input file, e.g. data_path("scplus", "eRegulon_direct.tsv")."""
    return DATA_DIR.joinpath(*parts)


def results_dir(*parts) -> Path:
    """Output directory (created if needed), e.g. results_dir("figure2", "perturbation")."""
    d = RESULTS_DIR.joinpath(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d
