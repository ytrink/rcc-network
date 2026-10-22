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

# -- Paths ----------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))   # src/utils
from paths import data_path, results_dir

INPUT_DIR   = results_dir("figure4")                    # written by figure4/R/01_umap.R and 02_regulatory_influence.R
MODEL_DIR   = results_dir("figure4", "perturbation")                  # regressors.pkl.gz
PERTURB_DIR = results_dir("figure4", "perturbation", "over_iter")     # <TF>_perturbation_over_iter.pkl.gz
PLOTS_DIR   = results_dir("figure4", "perturbation", "iteration_plots")
ERROR_LOG   = MODEL_DIR / "perturbation_errors.txt"

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
    pkl_path = str(PERTURB_DIR / f"{tf}_perturbation_over_iter.pkl.gz")
    print(f"Plotting {tf} - {config['program']}...")

    plot_program_ylim(
        tf                 = tf,
        targets            = config["targets"],
        pkl_path           = pkl_path,
        program_given_name = config["program"],
        base_plot_dir      = str(PLOTS_DIR),
        ylim               = (-1.2, 1.2)
    )