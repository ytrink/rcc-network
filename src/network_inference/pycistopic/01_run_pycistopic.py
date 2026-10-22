#!/usr/bin/env python3
"""
src/network_inference/pycistopic/01_run_pycistopic.py

Candidate regulatory regions for SCENIC+ with pycisTopic (v2.0a0).

IMPORTANT: the original pycisTopic run was performed interactively, following the official SCENIC+ /
pycisTopic tutorial (https://pycistopic.readthedocs.io); no script was saved. This file reconstructs the
tutorial steps with the settings that can be verified from the saved outputs on the analysis server
(/home/be/trinkya/projects/2024_09_18_ccRCC_YU/network/with_fragments_downsampling/outs/):

  - samples: six Yu et al. scATAC-seq samples, 2 per WHO/ISUP grade
      grade I: RCC86, RCC99 | grade II: RCC84, RCC100 | grade III: RCC87, RCC104
    (verified in cistopic_obj_2.pkl: 21,339 cells; QC logs outs/qc/<sample>.pycistopic_qc.log)
  - cell type labels: label-transferred Yu scATAC annotation (column "predicted.id"; 6 broad types:
    fibroblasts, macrophages/monocytes, tumor, T/NK, B cells, endothelial)
  - consensus peaks: MACS2 on per-cell-type pseudobulks (outs/consensus_peak_calling/consensus_regions.bed)
  - QC: "pycistopic qc" CLI per sample (outs/qc/)
  - topic modelling: LDA (Mallet), model selected from outs/models.pkl -> cistopic_obj_2.pkl
  - region sets for pycisTarget: top 3,000 regions per topic, Otsu-binarized topics, and DARs
    (02_run_dars.py; adj. p < 0.05, log2FC > log2(1.5)) -> outs/region_sets/

Parameters not listed above were the tutorial defaults. Fill in the inputs below to re-run.
"""

import os
import pickle
import subprocess

import pandas as pd
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models
from pycisTopic.topic_binarization import binarize_topics

# -- Inputs (edit) ---------------------------------------------------------------
SAMPLES        = ["RCC84", "RCC86", "RCC87", "RCC99", "RCC100", "RCC104"]
FRAGMENTS_DIR  = "data/yu/scatac/fragments"            # <sample>/fragments.tsv.gz (GEO GSE207493)
CELL_DATA_CSV  = "data/yu/scatac/cell_data.csv"        # barcode, sample_id, predicted.id (scATAC annotation)
CHROMSIZES     = "data/annotation/hg38.chrom.sizes"
BLACKLIST      = "data/annotation/hg38-blacklist.v2.bed"
TSS_BED        = "data/annotation/hg38_tss.bed"
MALLET_PATH    = "mallet"
OUT_DIR        = "results/network_inference/pycistopic"
N_TOPICS       = [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]   # tutorial grid; the selected model is in models.pkl

os.makedirs(OUT_DIR, exist_ok=True)
fragments = {s: os.path.join(FRAGMENTS_DIR, s, "fragments.tsv.gz") for s in SAMPLES}
cell_data = pd.read_csv(CELL_DATA_CSV, index_col=0)
cell_data = cell_data[cell_data["sample_id"].isin(SAMPLES)]
chromsizes = pd.read_table(CHROMSIZES, header=None, names=["Chromosome", "End"])
chromsizes.insert(1, "Start", 0)

# -- 1. Pseudobulk per cell type -> MACS2 -> consensus peaks ----------------------
bw_paths, bed_paths = export_pseudobulk(
    input_data=cell_data, variable="predicted.id", sample_id_col="sample_id",
    chromsizes=chromsizes, bed_path=f"{OUT_DIR}/pseudobulk_bed", bigwig_path=f"{OUT_DIR}/pseudobulk_bw",
    path_to_fragments=fragments, n_cpu=8, temp_dir="/tmp")
narrow_peaks = peak_calling(macs_path="macs2", bed_paths=bed_paths, outdir=f"{OUT_DIR}/macs",
                            genome_size="hs", n_cpu=8, temp_dir="/tmp")
consensus = get_consensus_peaks(narrow_peaks_dict=narrow_peaks, peak_half_width=250,
                                chromsizes=chromsizes, path_to_blacklist=BLACKLIST)
consensus_bed = f"{OUT_DIR}/consensus_regions.bed"
consensus.to_bed(path=consensus_bed, keep=True, compression="infer", chain=False)

# -- 2. QC (pycisTopic CLI, one call per sample) ----------------------------------
for s in SAMPLES:
    subprocess.run(["pycistopic", "qc", "--fragments", fragments[s], "--regions", consensus_bed,
                    "--tss", TSS_BED, "--output", f"{OUT_DIR}/qc/{s}"], check=True)

# -- 3. cisTopic object (annotated cells) -----------------------------------------
objs = []
for s in SAMPLES:
    bcs = list(cell_data.index[cell_data["sample_id"] == s])
    objs.append(create_cistopic_object_from_fragments(
        path_to_fragments=fragments[s], path_to_regions=consensus_bed, path_to_blacklist=BLACKLIST,
        valid_bc=bcs, n_cpu=8, project=s))
cistopic_obj = merge(objs)
cistopic_obj.add_cell_data(cell_data)

# -- 4. LDA topic modelling and model selection -----------------------------------
models = run_cgs_models_mallet(cistopic_obj, n_topics=N_TOPICS, n_cpu=12, n_iter=500,
                               random_state=555, alpha=50, alpha_by_topic=True, eta=0.1,
                               eta_by_topic=False, mallet_path=MALLET_PATH, tmp_path="/tmp")
model = evaluate_models(models, select_model=None, return_model=True)
cistopic_obj.add_LDA_model(model)
with open(f"{OUT_DIR}/cistopic_obj_2.pkl", "wb") as f:
    pickle.dump(cistopic_obj, f)

# -- 5. Topic region sets (top 3,000 and Otsu) ------------------------------------
region_bin_topics_top3k = binarize_topics(cistopic_obj, method="ntop", ntop=3000)
region_bin_topics_otsu = binarize_topics(cistopic_obj, method="otsu")
for name, topics in [("Topics_top_3k", region_bin_topics_top3k), ("Topics_otsu", region_bin_topics_otsu)]:
    os.makedirs(f"{OUT_DIR}/region_sets/{name}", exist_ok=True)
    for topic, regions in topics.items():
        pd.Series(regions.index).str.replace(":", "\t").str.replace("-", "\t") \
            .to_csv(f"{OUT_DIR}/region_sets/{name}/{topic}.bed", index=False, header=False)
