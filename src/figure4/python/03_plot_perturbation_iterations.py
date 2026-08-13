#!/usr/bin/env python3
"""
figures/figure4/python/perturbation/03_plot_perturbation_iterations.py

Plots predicted gene expression changes over perturbation iterations
for MEF2C perturbation in endothelial cells.
Requires outputs from 02_run_perturbation.py.
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT_DIR    = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
UTILS_DIR   = os.path.join(ROOT_DIR, "utils")
FIG_DIR     = os.path.join(ROOT_DIR, "figures", "figure4", "python", "perturbation")
PERTURB_DIR = os.path.join(FIG_DIR, "results")
PLOTS_DIR   = os.path.join(FIG_DIR, "iteration_plots")

os.makedirs(PLOTS_DIR, exist_ok=True)
sys.path.append(UTILS_DIR)
from iteration_func import plot_program_ylim

# ── Gene sets ──────────────────────────────────────────────────────────────────
ANGIOGENESIS_GENES = ["ANGPT2", "PLVAP", "MCAM", "THY1"]
VASCULAR_ECM_GENES = ["COL4A1", "COL4A2", "LAMC1", "NID1"]

# ── TF configurations ──────────────────────────────────────────────────────────
TF_CONFIGS = [
    {
        "tf":      "MEF2C",
        "targets": ANGIOGENESIS_GENES,
        "program": "angiogenesis"
    },
    {
        "tf":      "MEF2C",
        "targets": VASCULAR_ECM_GENES,
        "program": "vascular_ecm"
    }
]

# ── Plot ───────────────────────────────────────────────────────────────────────
for config in TF_CONFIGS:
    tf       = config["tf"]
    pkl_path = os.path.join(PERTURB_DIR, f"{tf}_perturbation_over_iter.pkl.gz")
    print(f"Plotting {tf} - {config['program']}...")

    plot_program_ylim(
        tf                 = tf,
        targets            = config["targets"],
        pkl_path           = pkl_path,
        program_given_name = config["program"],
        base_plot_dir      = PLOTS_DIR,
        ylim               = (-1.2, 1.2)
    )