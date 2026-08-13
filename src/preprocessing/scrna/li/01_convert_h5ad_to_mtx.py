#!/usr/bin/env python3
"""
preprocessing/scrna/li/01_convert_h5ad_to_mtx.py

Converts Li et al. h5ad object to Matrix Market format for loading in R.
Adds patient-prefixed barcodes and exports metadata, row, and column names.

Input:  data/li/RCC_upload_final_raw_counts.h5ad
Output: data/li/matrix/adata_X.mtx
        data/li/matrix/metadata.csv
        data/li/matrix/barcodes.csv
        data/li/matrix/genes.csv
"""

import os
import scipy.io
import scanpy as sc
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DATA_DIR = os.path.join(ROOT_DIR, "data", "li")
OUT_DIR  = os.path.join(DATA_DIR, "matrix")

os.makedirs(OUT_DIR, exist_ok=True)

# ── Load ───────────────────────────────────────────────────────────────────────
adata = sc.read_h5ad(os.path.join(DATA_DIR, "RCC_upload_final_raw_counts.h5ad"))

# ── Add patient-prefixed barcodes ──────────────────────────────────────────────
adata.obs_names = [
    "LI_" + adata.obs.loc[bc, "patient"] + "_" + bc
    for bc in adata.obs_names
]
adata.obs["barcode"] = adata.obs_names

# ── Export ─────────────────────────────────────────────────────────────────────
scipy.io.mmwrite(os.path.join(OUT_DIR, "adata_X.mtx"), adata.X.tocsr())

adata.obs.to_csv(os.path.join(OUT_DIR, "metadata.csv"), index=False)
adata.obs_names.to_series().to_csv(os.path.join(OUT_DIR, "barcodes.csv"), index=False, header=False)
adata.var_names.to_series().to_csv(os.path.join(OUT_DIR, "genes.csv"),    index=False, header=False)

print(f"Exported matrix and metadata -> {OUT_DIR}")
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")