#!/usr/bin/env python3
"""
figures/figure3/python/perturbation/03_plot_perturbation_iterations.py

Plots predicted gene expression changes over perturbation iterations
for EOMES, BATF, IKZF3, and combined EOMES+BATF+IKZF3 perturbations.
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

INPUT_DIR   = data_path("perturbation", "CD8")        # raw-count matrices exported by preprocessing/monocle (see README)
MODEL_DIR   = results_dir("figure3", "perturbation")                  # regressors.pkl.gz
PERTURB_DIR = results_dir("figure3", "perturbation", "over_iter")     # <TF>_perturbation_over_iter.pkl.gz
PLOTS_DIR   = results_dir("figure3", "perturbation", "iteration_plots")
ERROR_LOG   = MODEL_DIR / "perturbation_errors.txt"

# ── TF configurations ──────────────────────────────────────────────────────────
TF_CONFIGS = [
    {
        "tf":    "EOMES",
        "genes": ["LAG3", "TIGIT", "TOX", "BATF"],
        "gene_colors": {
            "LAG3":  "#1f77b4",
            "TIGIT": "#ff7f0e",
            "TOX":   "#2ca02c",
            "BATF":  "#d62728"
        }
    },
    {
        "tf":    "BATF",
        "genes": ["LAG3", "TIGIT", "TOX", "EOMES"],
        "gene_colors": {
            "LAG3":  "#1f77b4",
            "TIGIT": "#ff7f0e",
            "TOX":   "#2ca02c",
            "EOMES": "#d62728"
        }
    },
    {
        "tf":    "IKZF3",
        "genes": ["LAG3", "TIGIT", "TOX"],
        "gene_colors": {
            "LAG3":  "#1f77b4",
            "TIGIT": "#ff7f0e",
            "TOX":   "#2ca02c"
        }
    },
    {
        "tf":    "EOMES_BATF_IKZF3",
        "genes": ["LAG3", "TIGIT", "TOX"],
        "gene_colors": {
            "LAG3":  "#1f77b4",
            "TIGIT": "#ff7f0e",
            "TOX":   "#2ca02c"
        }
    }
]

# ── Plot ───────────────────────────────────────────────────────────────────────
for config in TF_CONFIGS:
    tf       = config["tf"]
    pkl_path = str(PERTURB_DIR / f"{tf}_perturbation_over_iter.pkl.gz")
    print(f"Plotting {tf}...")

    plot_program_ylim(
        tf                 = tf,
        targets            = config["genes"],
        pkl_path           = pkl_path,
        program_given_name = tf,
        base_plot_dir      = str(PLOTS_DIR),
        ylim               = (-0.8, 0.8),
        gene_colors        = config["gene_colors"]
    )