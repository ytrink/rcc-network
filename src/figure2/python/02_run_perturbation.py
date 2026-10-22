#!/usr/bin/env python3
"""
figures/figure2/python/perturbation/02_run_perturbation.py

Runs iterative TF perturbation simulations using pre-trained regressors.
Outputs per-TF perturbation results to disk for plotting in 03_plot_perturbation_iterations.py.

Note: This script is computationally intensive. Run after 01_get_regressor.py.
"""

import os
import gzip
import pickle
from datetime import datetime

import pandas as pd
import sys

# -- Paths ----------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))   # src/utils
from paths import data_path, results_dir

INPUT_DIR   = data_path("perturbation", "mac_mono")   # raw-count matrices exported by preprocessing/monocle (see README)
MODEL_DIR   = results_dir("figure2", "perturbation")                  # regressors.pkl.gz
PERTURB_DIR = results_dir("figure2", "perturbation", "over_iter")     # <TF>_perturbation_over_iter.pkl.gz
PLOTS_DIR   = results_dir("figure2", "perturbation", "arrow_plots")
ERROR_LOG   = MODEL_DIR / "perturbation_errors.txt"

# ── Parameters ─────────────────────────────────────────────────────────────────
CELL_TYPE = "mac_mono"

TFS = [
    "TCF7L2", "CEBPD", "MITF", "ATF3", "FOS", "CREB5", "HIF1A",
    "MEF2A", "ZFHX3", "USF2", "ETV6", "JDP2", "YBX1", "KLF4",
    "IRF8", "IRF5", "SPI1", "CEBPB", "CEBPA", "MAFB", "POU2F2",
    "REL", "HMGA1", "NR4A3", "MAF"
]


# ── Helpers ────────────────────────────────────────────────────────────────────
def timestamp(msg: str):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

# ── Load data ──────────────────────────────────────────────────────────────────
timestamp("Loading regressors...")
with gzip.open(INPUT_DIR / "regressors.pkl.gz", "rb") as f:
    regressors = pickle.load(f)
timestamp("Regressors loaded.")

timestamp("Loading expression matrix...")
X = pd.read_csv(
    INPUT_DIR / f"cds_{CELL_TYPE}_sampled10k.csv",
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