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

# -- Paths ----------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))   # src/utils
from paths import data_path, results_dir

INPUT_DIR   = data_path("perturbation", "CD8")        # raw-count matrices exported by preprocessing/monocle (see README)
MODEL_DIR   = results_dir("figure3", "perturbation")                  # regressors.pkl.gz
PERTURB_DIR = results_dir("figure3", "perturbation", "over_iter")     # <TF>_perturbation_over_iter.pkl.gz
PLOTS_DIR   = results_dir("figure3", "perturbation", "arrow_plots")
ERROR_LOG   = MODEL_DIR / "perturbation_errors.txt"

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
        data_path("li", "CD8_expressed_genes_10perc.csv"),
        index_col=0
    )["x"]
)
tfs = list(set(TFS).intersection(expressed_genes))

# ── Load regressors ────────────────────────────────────────────────────────────
timestamp("Loading regressors...")
with gzip.open(MODEL_DIR / "regressors.pkl.gz", "rb") as f:
    regressors = pickle.load(f)
timestamp("Regressors loaded.")

# ── Load expression matrix ─────────────────────────────────────────────────────
timestamp("Loading expression matrix...")
X = pd.read_csv(
    INPUT_DIR / "cds_normalized_tcells_sampled_7500.csv",
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
        out_folder = str(PERTURB_DIR)
    )

timestamp("Done.")