#!/usr/bin/env python3
"""
figures/figure4/python/perturbation/01_get_regressor.py

Trains gene expression regression models for TF perturbation analysis
of normal-to-tumor endothelial cell transition.
Models are saved to disk for use in 02_run_perturbation.py.

Note: This script is computationally intensive and only needs to be run once.
"""

import os
import gzip
import pickle

import pandas as pd
from scenicplus.simulation import train_gene_expression_models

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DATA_DIR = os.path.join(ROOT_DIR, "data")
OUT_DIR  = os.path.join(ROOT_DIR, "figures", "figure4", "python", "perturbation")

os.makedirs(OUT_DIR, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────────────────────
CELL_TYPE = "endo"

TFS = [
    "FLI1", "SOX7", "EBF1", "RARB", "ELF1", "GATA2", "SOX17",
    "NFIA", "ZBTB20", "PBX1", "ETS1", "NFIB", "KLF9", "ID2",
    "HES1", "ERG", "ELK3", "MEF2C", "ETS2", "JUN", "EGR1",
    "KLF4", "FOS", "ATF3", "CEBPD", "CEBPB", "REL"
]

# ── Load DEGs ──────────────────────────────────────────────────────────────────
degs1 = list(pd.read_csv(os.path.join(DATA_DIR, CELL_TYPE, "degs_hec.csv"), index_col=0)["x"])
degs2 = list(pd.read_csv(os.path.join(DATA_DIR, CELL_TYPE, "degs_tec.csv"), index_col=0)["x"])
degs  = degs1 + degs2

# ── Load expression matrix ─────────────────────────────────────────────────────
X = pd.read_csv(
    os.path.join(DATA_DIR, CELL_TYPE, "X_counts.csv"),
    index_col=0
).T

# ── Load eRegulon gene-TF mappings ─────────────────────────────────────────────
direct_eregulons   = pd.read_table(os.path.join(DATA_DIR, "scplus", "eRegulon_direct.tsv"))
extended_eregulons = pd.read_table(os.path.join(DATA_DIR, "scplus", "eRegulon_extended.tsv"))

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
out_path = os.path.join(OUT_DIR, "regressors.pkl.gz")
with gzip.open(out_path, "wb") as f:
    pickle.dump(regressors, f)

print(f"Saved regressors -> {out_path}")