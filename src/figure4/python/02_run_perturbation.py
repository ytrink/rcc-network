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

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT_DIR  = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DATA_DIR  = os.path.join(ROOT_DIR, "data")
UTILS_DIR = os.path.join(ROOT_DIR, "utils")
OUT_DIR   = os.path.join(ROOT_DIR, "figures", "figure4", "python", "perturbation", "results")

os.makedirs(OUT_DIR, exist_ok=True)
sys.path.append(UTILS_DIR)
import perturbation_functions as pf

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
with gzip.open(os.path.join(OUT_DIR, "..", "regressors.pkl.gz"), "rb") as f:
    regressors = pickle.load(f)
timestamp("Regressors loaded.")

# ── Load expression matrix ─────────────────────────────────────────────────────
timestamp("Loading expression matrix...")
X = pd.read_csv(
    os.path.join(DATA_DIR, CELL_TYPE, "X_counts.csv"),
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
        out_folder = OUT_DIR
    )

timestamp("Done.")