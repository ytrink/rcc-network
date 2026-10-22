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

# -- Paths ----------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))   # src/utils
from paths import data_path, results_dir

INPUT_DIR   = data_path("perturbation", "mac_mono")   # raw-count matrices exported by preprocessing/monocle (see README)
MODEL_DIR   = results_dir("figure2", "perturbation")                  # regressors.pkl.gz
PERTURB_DIR = results_dir("figure2", "perturbation", "over_iter")     # <TF>_perturbation_over_iter.pkl.gz
PLOTS_DIR   = results_dir("figure2", "perturbation", "iteration_plots")
ERROR_LOG   = MODEL_DIR / "perturbation_errors.txt"

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
    pkl_path = str(PERTURB_DIR / f"{tf}_perturbation_over_iter.pkl.gz")
    print(f"Plotting {tf}...")

    plot_program_ylim(
        tf                 = tf,
        targets            = config["genes"],
        pkl_path           = pkl_path,
        program_given_name = tf,
        base_plot_dir      = str(PLOTS_DIR),
        ylim               = (-1.2, 1.2),
        gene_colors        = config["gene_colors"]
    )