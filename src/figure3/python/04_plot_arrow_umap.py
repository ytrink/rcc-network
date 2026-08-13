#!/usr/bin/env python3
"""
figures/figure3/python/perturbation/04_plot_arrow_umap.py

Generates arrow UMAP plots showing predicted cell state shifts after TF perturbation
for CD8+ T cell exhaustion analysis.
Requires outputs from 02_run_perturbation.py.
"""

import os
import sys
from datetime import datetime

import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT_DIR  = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DATA_DIR  = os.path.join(ROOT_DIR, "data")
UTILS_DIR = os.path.join(ROOT_DIR, "utils")
FIG_DIR   = os.path.join(ROOT_DIR, "figures", "figure3", "python", "perturbation")

PERTURB_DIR = os.path.join(FIG_DIR, "results")
PLOTS_DIR   = os.path.join(FIG_DIR, "arrow_plots")
ERROR_LOG   = os.path.join(FIG_DIR, "perturbation_errors.txt")

os.makedirs(PLOTS_DIR, exist_ok=True)
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

# ── Load data ──────────────────────────────────────────────────────────────────
umap_df            = pd.read_csv(os.path.join(DATA_DIR, CELL_TYPE, "cds_tcells_umap_sampled_7500.csv"), index_col=0)
direct_eregulons   = pd.read_table(os.path.join(DATA_DIR, "scplus", "eRegulon_direct.tsv"))
extended_eregulons = pd.read_table(os.path.join(DATA_DIR, "scplus", "eRegulon_extended.tsv"))

# ── Skip TFs already processed ────────────────────────────────────────────────
tfs_done   = [name for name in os.listdir(PLOTS_DIR) if os.path.isdir(os.path.join(PLOTS_DIR, name))]
tfs_to_run = list(set(tfs).difference(set(tfs_done)))
timestamp(f"Running arrow plots for {len(tfs_to_run)} TFs ({len(tfs_done)} already done)")

# ── Run arrow plots ────────────────────────────────────────────────────────────
for tf in tfs_to_run:
    try:
        pf.run_arrow_analysis(
            umap_df            = umap_df,
            tf                 = tf,
            direct_eregulons   = direct_eregulons,
            extended_eregulons = extended_eregulons,
            plots_folder       = os.path.join(PLOTS_DIR, tf),
            out_folder         = PERTURB_DIR
        )
    except Exception as e:
        error_msg = (
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
            f"Failed for TF: {tf} — {type(e).__name__}: {e}"
        )
        print(error_msg)
        with open(ERROR_LOG, "a") as log_file:
            log_file.write(error_msg + "\n")