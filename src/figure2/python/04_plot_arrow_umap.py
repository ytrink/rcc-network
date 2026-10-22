#!/usr/bin/env python3
"""
figures/figure2/python/perturbation/04_plot_arrow_umap.py

Generates arrow UMAP plots showing predicted cell state shifts after TF perturbation.
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
umap_df            = pd.read_csv(INPUT_DIR / f"cds_{CELL_TYPE}_umap_sampled10k.csv", index_col=0)
direct_eregulons   = pd.read_table(data_path("scplus", "eRegulon_direct.tsv"))
extended_eregulons = pd.read_table(data_path("scplus", "eRegulons_extended.tsv"))

# ── Skip TFs already processed ────────────────────────────────────────────────
tfs_done = [
    name for name in os.listdir(PLOTS_DIR)
    if os.path.isdir(str(PLOTS_DIR / name))
]
tfs_to_run = list(set(TFS).difference(set(tfs_done)))
timestamp(f"Running arrow plots for {len(tfs_to_run)} TFs ({len(tfs_done)} already done)")

# ── Run arrow plots ────────────────────────────────────────────────────────────
for tf in tfs_to_run:
    tf_plots_dir = str(PLOTS_DIR / tf)
    try:
        pf.run_arrow_analysis(
            umap_df            = umap_df,
            tf                 = tf,
            direct_eregulons   = direct_eregulons,
            extended_eregulons = extended_eregulons,
            plots_folder       = tf_plots_dir,
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