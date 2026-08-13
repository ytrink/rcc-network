"""
utils/perturbation_functions.py

Utility functions for TF perturbation simulations and visualization using SCENIC+.
Called by figures/figure2,3,4/python/perturbation/ scripts.
"""

import os
import gc
import gzip
import pickle

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scenicplus.simulation import (
    simulate_perturbation,
    plot_perturbation_effect_in_embedding
)

# ── Simulation ─────────────────────────────────────────────────────────────────

def run_perturbation_sim(exprMatrix, tf, regressors, out_folder, n_iter=5):
    """
    Run iterative TF perturbation simulation and save results as gzipped pickle.
    Skips simulation if output file already exists.

    Args:
        exprMatrix : pd.DataFrame — cells x genes expression matrix
        tf         : str — TF to perturb (set to 0)
        regressors : dict — pre-trained regressors from train_gene_expression_models()
        out_folder : str — directory to save output
        n_iter     : int — number of simulation iterations (default: 5)
    """
    os.makedirs(out_folder, exist_ok=True)
    out_path = os.path.join(out_folder, f"{tf}_perturbation_over_iter.pkl.gz")

    if os.path.exists(out_path):
        print(f"Output already exists for {tf}, skipping.")
        return

    print(f"Simulating perturbation for {tf}...")
    perturbation_over_iter = simulate_perturbation(
        df_EXP           = exprMatrix,
        perturbation      = {tf: 0},
        keep_intermediate = True,
        n_iter            = n_iter,
        regressors        = regressors
    )

    with gzip.open(out_path, "wb") as f:
        pickle.dump(perturbation_over_iter, f)

    print(f"Saved -> {out_path}")


# ── Arrow UMAP plot ────────────────────────────────────────────────────────────

def run_arrow_analysis(umap_df, tf, direct_eregulons, extended_eregulons,
                       plots_folder, out_folder):
    """
    Generate arrow UMAP plot showing predicted cell state shifts after TF perturbation.

    Args:
        umap_df            : pd.DataFrame — UMAP coordinates (cells x 2)
        tf                 : str — TF name
        direct_eregulons   : pd.DataFrame — direct eRegulon table from SCENIC+
        extended_eregulons : pd.DataFrame — extended eRegulon table from SCENIC+
        plots_folder       : str — directory to save plots
        out_folder         : str — directory containing perturbation pkl.gz files
    """
    os.makedirs(plots_folder, exist_ok=True)

    print(f"Loading perturbation results for {tf}...")
    with gzip.open(os.path.join(out_folder, f"{tf}_perturbation_over_iter.pkl.gz"), "rb") as f:
        perturbation_over_iter = pickle.load(f)

    perturbation_over_iter[5] = perturbation_over_iter[5].astype("float64")
    perturbation_over_iter[0] = perturbation_over_iter[0].astype("float64")

    embedding = umap_df.iloc[:, :2].to_numpy()

    eRegulon_data = pd.concat([
        direct_eregulons.drop_duplicates(),
        extended_eregulons.drop_duplicates()
    ]).drop_duplicates()

    print(f"Plotting arrow UMAP for {tf}...")
    fig, ax = plt.subplots()
    ax.scatter(umap_df.iloc[:, 0], umap_df.iloc[:, 1], alpha=0.3)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    plt.tight_layout()

    plot_perturbation_effect_in_embedding(
        perturbed_matrix = perturbation_over_iter[5],
        original_matrix  = perturbation_over_iter[0],
        embedding        = embedding,
        AUC_kwargs       = {},
        ax               = ax,
        eRegulons        = eRegulon_data,
        n_cpu            = 1
    )

    for ext in ("pdf", "svg"):
        plt.savefig(os.path.join(plots_folder, f"{tf}_perturbation_arrowplot.{ext}"))

    plt.close()
    gc.collect()
    print(f"Saved arrow plot for {tf}")