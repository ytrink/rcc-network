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

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT_DIR    = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
UTILS_DIR   = os.path.join(ROOT_DIR, "utils")
FIG_DIR     = os.path.join(ROOT_DIR, "figures", "figure3", "python", "perturbation")
PERTURB_DIR = os.path.join(FIG_DIR, "results")
PLOTS_DIR   = os.path.join(FIG_DIR, "iteration_plots")

os.makedirs(PLOTS_DIR, exist_ok=True)
sys.path.append(UTILS_DIR)
from iteration_func import plot_program_ylim

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
    pkl_path = os.path.join(PERTURB_DIR, f"{tf}_perturbation_over_iter.pkl.gz")
    print(f"Plotting {tf}...")

    plot_program_ylim(
        tf                 = tf,
        targets            = config["genes"],
        pkl_path           = pkl_path,
        program_given_name = tf,
        base_plot_dir      = PLOTS_DIR,
        ylim               = (-0.8, 0.8),
        gene_colors        = config["gene_colors"]
    )