#!/usr/bin/env python3
"""
figures/figure2/python/perturbation/03_plot_perturbation_iterations.py

Plots predicted gene expression changes over perturbation iterations for MAFB and SPI1.
Requires outputs from 02_run_perturbation.py.
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT_DIR    = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
UTILS_DIR   = os.path.join(ROOT_DIR, "utils")
FIG_DIR     = os.path.join(ROOT_DIR, "figures", "figure2", "python", "perturbation")
PERTURB_DIR = os.path.join(FIG_DIR, "results")
PLOTS_DIR   = os.path.join(FIG_DIR, "iteration_plots")

os.makedirs(PLOTS_DIR, exist_ok=True)
sys.path.append(UTILS_DIR)
from iteration_func import plot_program_ylim

# ── TF configurations ──────────────────────────────────────────────────────────
TF_CONFIGS = [
    {
        "tf":    "MAFB",
        "genes": ["ARHGAP18", "CD74", "APOC1", "APOE", "TREM2", "LGMN", "CD163"],
        "gene_colors": {
            "ARHGAP18": "#1f77b4",
            "CD74":     "#ff7f0e",
            "APOC1":    "#2ca02c",
            "APOE":     "#d62728",
            "TREM2":    "#9467bd",
            "LGMN":     "#8c564b",
            "CD163":    "#17becf"
        }
    },
    {
        "tf":    "SPI1",
        "genes": ["ARHGAP18", "CD74", "APOE", "TREM2", "MSR1", "LGMN", "CD68"],
        "gene_colors": {
            "ARHGAP18": "#1f77b4",
            "CD74":     "#ff7f0e",
            "APOE":     "#d62728",
            "TREM2":    "#9467bd",
            "MSR1":     "#2ca02c",
            "LGMN":     "#8c564b",
            "CD68":     "#17becf"
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
        ylim               = (-1.2, 1.2),
        gene_colors        = config["gene_colors"]
    )