#!/usr/bin/env python3
"""
figures/figure4/python/perturbation/02_run_perturbation.py

Runs iterative TF perturbation simulations for normal-to-tumor
endothelial cell transition analysis.
Requires regressors from 01_get_regressor.py.

Note: Computationally intensive. Run after 01_get_regressor.py.
"""

import os
import sys
import gzip
import pickle
from datetime import datetime

import pandas as pd

# -- Paths ----------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))   # src/utils
from paths import data_path, results_dir

INPUT_DIR   = results_dir("figure4")                    # written by figure4/R/01_umap.R and 02_regulatory_influence.R
MODEL_DIR   = results_dir("figure4", "perturbation")                  # regressors.pkl.gz
PERTURB_DIR = results_dir("figure4", "perturbation", "over_iter")     # <TF>_perturbation_over_iter.pkl.gz
PLOTS_DIR   = results_dir("figure4", "perturbation", "arrow_plots")
ERROR_LOG   = MODEL_DIR / "perturbation_errors.txt"

# ── Parameters ─────────────────────────────────────────────────────────────────
CELL_TYPE = "endo"

TFS = [
    "FLI1", "SOX7", "EBF1", "RARB", "ELF1", "GATA2", "SOX17",
    "NFIA", "ZBTB20", "PBX1", "ETS1", "NFIB", "KLF9", "ID2",
    "HES1", "ERG", "ELK3", "MEF2C", "ETS2", "JUN", "EGR1",
    "KLF4", "FOS", "ATF3", "CEBPD", "CEBPB", "REL"
]

# ── Helpers ────────────────────────────────────────────────────────────────────
def timestamp(msg: str):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

# ── Load regressors ────────────────────────────────────────────────────────────
timestamp("Loading regressors...")
with gzip.open(MODEL_DIR / "regressors.pkl.gz", "rb") as f:
    regressors = pickle.load(f)
timestamp("Regressors loaded.")

# ── Load expression matrix ─────────────────────────────────────────────────────
timestamp("Loading expression matrix...")
X = pd.read_csv(
    INPUT_DIR / "X_counts.csv",
    index_col=0
).T
timestamp("Expression matrix loaded.")

# ── Run perturbation simulations ───────────────────────────────────────────────
timestamp(f"Running perturbation for {len(TFS)} TFs...")
for tf in TFS:
    pf.run_perturbation_sim(
        exprMatrix = X,
        tf         = tf,
        regressors = regressors,
        out_folder = str(PERTURB_DIR)
    )

timestamp("Done.")