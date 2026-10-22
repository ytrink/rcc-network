#!/usr/bin/env python3
"""
figures/figure3/python/perturbation/01_get_regressor.py

Trains gene expression regression models for TF perturbation analysis
of CD8+ T cell exhaustion transition.
Models are saved to disk for use in 02_run_perturbation.py.

Note: This script is computationally intensive and only needs to be run once.
"""

import os
import gzip
import pickle

import pandas as pd
from scenicplus.simulation import train_gene_expression_models

# -- Paths ----------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))   # src/utils
from paths import data_path, results_dir

INPUT_DIR   = data_path("perturbation", "CD8")        # raw-count matrices exported by preprocessing/monocle (see README)
MODEL_DIR   = results_dir("figure3", "perturbation")                  # regressors.pkl.gz
PERTURB_DIR = results_dir("figure3", "perturbation", "over_iter")     # <TF>_perturbation_over_iter.pkl.gz
PLOTS_DIR   = results_dir("figure3", "perturbation", "arrow_plots")
ERROR_LOG   = MODEL_DIR / "perturbation_errors.txt"

# ── Parameters ─────────────────────────────────────────────────────────────────
CELL_TYPE = "CD8"

TFS = [
    "FOSB", "BHLHE40", "RUNX1", "JUNB", "PRDM1", "RAD21",
    "IRF2", "ID2", "GATA3", "NFAT5", "FOXP1", "NR3C1",
    "ILF2", "MAZ", "EOMES", "BATF", "BCL11B", "IKZF3",
    "IKZF1", "RUNX3", "TBX21"
]

# ── Load data ──────────────────────────────────────────────────────────────────
# DEGs
degs_df = pd.read_csv(data_path("li", "degs_per_cell_type.csv"))
degs    = list(degs_df.loc[degs_df["cluster"] == "T_NK", "gene"])

# Filter TFs to those expressed in >=10% of cells
expressed_genes = set(
    pd.read_csv(
        data_path("li", "CD8_expressed_genes_10perc.csv"),
        index_col=0
    )["x"]
)
tfs = list(set(TFS).intersection(expressed_genes))

# Expression matrix
X = pd.read_csv(
    INPUT_DIR / "cds_normalized_tcells_sampled_7500.csv",
    index_col=0
).T

# eRegulon gene-TF mappings
direct_eregulons   = pd.read_table(data_path("scplus", "eRegulon_direct.tsv"))
extended_eregulons = pd.read_table(data_path("scplus", "eRegulons_extended.tsv"))

gene_tf = pd.concat([
    direct_eregulons[["Gene", "TF"]].drop_duplicates(),
    extended_eregulons[["Gene", "TF"]].drop_duplicates()
]).drop_duplicates()

gene_to_TF_raw = gene_tf.groupby("Gene")["TF"].apply(list).to_dict()

# ── Filter genes and TFs to those present in expression matrix ─────────────────
degs_filtered = set(degs).intersection(set(gene_to_TF_raw.keys()))

gene_to_TF_filtered = {
    k: [tf for tf in v if tf in X.columns]
    for k, v in gene_to_TF_raw.items()
    if k in X.columns and any(tf in X.columns for tf in v)
}

# ── Train regressors ───────────────────────────────────────────────────────────
regressors = train_gene_expression_models(
    df_EXP     = X,
    gene_to_TF = gene_to_TF_filtered,
    genes      = degs_filtered
)

# ── Save ───────────────────────────────────────────────────────────────────────
out_path = MODEL_DIR / "regressors.pkl.gz"
with gzip.open(out_path, "wb") as f:
    pickle.dump(regressors, f)

print(f"Saved regressors -> {out_path}")