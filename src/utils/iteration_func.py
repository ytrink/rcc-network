"""
utils/iteration_func.py

Utility functions for plotting TF perturbation results over iterations.
Called by figures/figure2,3,4/python/perturbation/03_plot_perturbation_iterations.py
"""

import os
import gzip
import pickle

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_program_ylim(tf, targets, pkl_path, program_given_name, base_plot_dir,
                      ylim=None, gene_colors=None):
    """
    Load perturbation results and plot predicted expression changes over iterations.

    Args:
        tf                 : str — TF name
        targets            : list[str] — target genes to plot
        pkl_path           : str — path to perturbation pkl.gz file
        program_given_name : str — label for the biological program
        base_plot_dir      : str — output directory
        ylim               : tuple(float, float) or None — y-axis limits e.g. (-1.2, 1.2)
        gene_colors        : dict or None — {gene: color} mapping
    """
    print(f"Loading perturbation results for {tf}...")
    with gzip.open(pkl_path, "rb") as f:
        perturbation_over_iter = pickle.load(f)

    plots_folder = os.path.join(base_plot_dir, tf)
    _plot_iterations(
        perturbation_over_iter = perturbation_over_iter,
        plots_folder           = plots_folder,
        genes_to_show          = targets,
        tf                     = tf,
        program_name           = program_given_name,
        ylim                   = ylim,
        gene_colors            = gene_colors
    )


def _plot_iterations(perturbation_over_iter, plots_folder, genes_to_show, tf,
                     program_name, ylim=None, gene_colors=None):
    """
    Plot log2FC over perturbation iterations for a set of target genes.

    Args:
        perturbation_over_iter : list[pd.DataFrame] — simulated expression per iteration
        plots_folder           : str — output directory
        genes_to_show          : list[str] — genes to plot
        tf                     : str — TF name
        program_name           : str — biological program label
        ylim                   : tuple or None — y-axis limits
        gene_colors            : dict or None — {gene: color} mapping
    """
    os.makedirs(plots_folder, exist_ok=True)

    baseline = perturbation_over_iter[0][genes_to_show].mean()
    max_iter = min(6, len(perturbation_over_iter))

    fig, ax = plt.subplots(figsize=(5, 5))

    for gene in genes_to_show:
        log2fc = [
            np.log2(perturbation_over_iter[i][gene].mean() / baseline[gene])
            for i in range(1, max_iter)
        ]
        color = gene_colors.get(gene) if gene_colors else None
        ax.plot(range(1, max_iter), log2fc, linewidth=2, label=gene, color=color)

    ax.set_ylabel("Predicted $\\log_2$FC")
    ax.set_xlabel("Iteration")
    ax.axhline(y=0, color="black", linestyle="--")
    ax.legend(frameon=False)
    ax.grid(True, linestyle=":", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_xticks(range(1, max_iter))
    if ylim is not None:
        ax.set_ylim(ylim)

    plt.title(f"Predicted Expression Change: {tf}")
    plt.tight_layout()

    out_path = os.path.join(plots_folder, f"{program_name}_{tf}_log2fc_over_iterations.svg")
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved -> {out_path}")