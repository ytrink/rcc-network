#!/usr/bin/env python3
"""
figures/figure3/python/perturbation/02_run_perturbation.py

Runs iterative TF perturbation simulations for CD8+ T cell exhaustion analysis.
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
OUT_DIR   = os.path.join(ROOT_DIR, "figures", "figure3", "python", "perturbation", "results")

os.makedirs(OUT_DIR, exist_ok=True)
sys.path.append(UTILS_DIR)
import perturbation_functions as pf

# ── Parameters ─────────────────────────────────────────────────────────────────
CELL_TYPE = "CD8"

TFS = [
    "FOSB", "BHLHE40", "RUNX1", "JUNB", "PRDM1", "RAD21",
    "IRF2", "ID2", "GATA3", "NFAT5", "FOXP1", "NR3C1",
    "ILF2", "MAZ", "EOMES", "BATF", "BCL11B", "IKZF3",
    "IKZF1", "RUNX3", "TBX21"
]

# ── Helpers ────────────────────────────────────────────────────────────────────
def timestamp(msg: str):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

# ── Filter to expressed TFs ────────────────────────────────────────────────────
expressed_genes = set(
    pd.read_csv(
        os.path.join(DATA_DIR, CELL_TYPE, "CD8_expressed_genes_10perc.csv"),
        index_col=0
    )["x"]
)
tfs = list(set(TFS).intersection(expressed_genes))

# ── Load regressors ────────────────────────────────────────────────────────────
timestamp("Loading regressors...")
with gzip.open(os.path.join(OUT_DIR, "..", "regressors.pkl.gz"), "rb") as f:
    regressors = pickle.load(f)
timestamp("Regressors loaded.")

# ── Load expression matrix ─────────────────────────────────────────────────────
timestamp("Loading expression matrix...")
X = pd.read_csv(
    os.path.join(DATA_DIR, CELL_TYPE, "cds_normalized_tcells_sampled_7500.csv"),
    index_col=0
).T
timestamp("Expression matrix loaded.")

# ── Run perturbation simulations ───────────────────────────────────────────────
timestamp(f"Running perturbation for {len(tfs)} TFs...")
for tf in tfs:
    pf.run_perturbation_sim(
        exprMatrix = X,
        tf         = tf,
        regressors = regressors,
        out_folder = OUT_DIR
    )

timestamp("Done.")