#!/usr/bin/env python3
"""
figures/figure4/python/perturbation/04_plot_arrow_umap.py

Generates arrow UMAP plots showing predicted cell state shifts after TF perturbation
for normal-to-tumor endothelial cell transition analysis.
Requires outputs from 02_run_perturbation.py.
"""

import os
import sys
from datetime import datetime

import pandas as pd

# -- Paths ----------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))   # src/utils
from paths import data_path, results_dir

INPUT_DIR   = results_dir("figure4")                    # written by figure4/R/01_umap.R and 02_regulatory_influence.R
MODEL_DIR   = results_dir("figure4", "perturbation")                  # regressors.pkl.gz
PERTURB_DIR = results_dir("figure4", "perturbation", "over_iter")     # <TF>_perturbation_over_iter.pkl.gz
PLOTS_DIR   = results_dir("figure4", "perturbation", "arrow_plots")
ERROR_LOG   = MODEL_DIR / "perturbation_errors.txt"

# ── Parameters ─────────────────────────────────────────────────────────────────
CELL_TYPE = "endo"

TFS = [
    "FLI1", "SOX7", "EBF1", "RARB", "ELF1", "GATA2", "SOX17",
    "NFIA", "ZBTB20", "PBX1", "ETS1", "NFIB", "KLF9", "ID2",
    "HES1", "ERG", "ELK3", "MEF2C", "ETS2", "JUN", "EGR1",
    "KLF4", "FOS", "ATF3", "CEBPD", "CEBPB", "REL"
]

# ── Helpers ────────────────────────────────────────────────────────────────────
def timestamp(msg: str):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

# ── Load data ──────────────────────────────────────────────────────────────────
umap_df = pd.read_csv(INPUT_DIR / "X_umap.csv", index_col=0)
umap_df.index = "X" + umap_df.index.astype(str)

direct_eregulons   = pd.read_table(data_path("scplus", "eRegulon_direct.tsv"))
extended_eregulons = pd.read_table(data_path("scplus", "eRegulons_extended.tsv"))

# ── Skip TFs already processed ────────────────────────────────────────────────
tfs_done   = [name for name in os.listdir(PLOTS_DIR) if os.path.isdir(str(PLOTS_DIR / name))]
tfs_to_run = list(set(TFS).difference(set(tfs_done)))
timestamp(f"Running arrow plots for {len(tfs_to_run)} TFs ({len(tfs_done)} already done)")

# ── Run arrow plots ────────────────────────────────────────────────────────────
for tf in tfs_to_run:
    try:
        pf.run_arrow_analysis(
            umap_df            = umap_df,
            tf                 = tf,
            direct_eregulons   = direct_eregulons,
            extended_eregulons = extended_eregulons,
            plots_folder       = str(PLOTS_DIR / tf),
            out_folder         = str(PERTURB_DIR)
        )
    except Exception as e:
        error_msg = (
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
            f"Failed for TF: {tf} — {type(e).__name__}: {e}"
        )
        print(error_msg)
        with open(ERROR_LOG, "a") as log_file:
            log_file.write(error_msg + "\n")