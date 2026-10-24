#!/usr/bin/env python3
"""
src/network_inference/pycistopic/01_run_pycistopic.py

Candidate regulatory regions for SCENIC+ with pycisTopic (v2.0a0), following the pycisTopic 2.0 tutorial.

The analysis was run interactively (June 2025, run folder network/with_fragments_downsampling/ on the analysis
server). This script consolidates that session:
  - steps 1, 2 and 5 are the commands from the session's Python history, with paths made relative;
  - steps 3-4 (building the cisTopic object and the LDA models) follow the pycisTopic 2.0 tutorial, with all
    parameters taken from the saved objects (models.pkl, cistopic_obj_2.pkl).
Checks against the saved outputs: 408,701 consensus regions (500 bp); 21,339 cells, identical per sample to the
cells passing the Otsu QC thresholds; LDA models with 2, 10, 25 and 50 topics (Mallet, 500 iterations,
random_state 555, alpha 50 by topic, eta 0.1), 25-topic model selected.

Samples: six Yu et al. scATAC-seq samples, 2 per WHO/ISUP grade
  (grade I: RCC86, RCC99 | grade II: RCC84, RCC100 | grade III: RCC87, RCC104).
DARs and the DAR region sets: 02_run_dars.py.

Inputs (edit paths below):
  data/yu/scatac/fragments/<sample>/fragments.tsv.gz  (GEO GSE207493)
  data/yu/scatac/cell_data.csv     label-transferred scATAC annotation (index "<sample>_<barcode>", "predicted.id")
  data/annotation/hg38-blacklist.v2.bed
  mallet binary (Mallet LDA)
Output: results/network_inference/pycistopic/
"""

import os
import pickle
import subprocess

import pandas as pd

OUT_DIR = "results/network_inference/pycistopic"
FRAGMENTS_DIR = "data/yu/scatac/fragments"
CELL_DATA_CSV = "data/yu/scatac/cell_data.csv"
BLACKLIST = "data/annotation/hg38-blacklist.v2.bed"
MALLET_PATH = "mallet"
TMP_DIR = "/tmp"

SAMPLES = ["RCC84", "RCC86", "RCC87", "RCC99", "RCC100", "RCC104"]
fragments_dict = {s: os.path.join(FRAGMENTS_DIR, s, "fragments.tsv.gz") for s in SAMPLES}
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Pseudobulks per cell type -> MACS2 -> consensus peaks (from the session history) ───────────
cell_data = pd.read_csv(CELL_DATA_CSV, index_col=0)
chromsizes = pd.read_table("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
                           header=None, names=["Chromosome", "End"])
chromsizes.insert(1, "Start", 0)

# barcodes "<sample>_<barcode>" -> "<barcode>-<sample>"
cell_data["barcode"] = list(cell_data.index)
cell_data["sample"] = cell_data["barcode"].str.extract(r"^(RCC\d{2,3})_")
cell_data["clean_barcode"] = cell_data["barcode"].str.replace(r"^RCC\d{2,3}_", "", regex=True)
cell_data.index = cell_data["clean_barcode"] + "-" + cell_data["sample"]
cell_data = cell_data.drop(columns=["barcode", "clean_barcode"])
cell_data = cell_data[cell_data["sample"].isin(fragments_dict.keys())]

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks

peak_dir = os.path.join(OUT_DIR, "consensus_peak_calling")
for d in ["pseudobulk_bed_files", "pseudobulk_bw_files", "MACS"]:
    os.makedirs(os.path.join(peak_dir, d), exist_ok=True)

bw_paths, bed_paths = export_pseudobulk(
    input_data=cell_data,
    variable="predicted.id",
    sample_id_col="sample",
    chromsizes=chromsizes,
    bed_path=os.path.join(peak_dir, "pseudobulk_bed_files"),
    bigwig_path=os.path.join(peak_dir, "pseudobulk_bw_files"),
    path_to_fragments=fragments_dict,
    n_cpu=2,
    normalize_bigwig=True,
    temp_dir=TMP_DIR,
    split_pattern="-")

narrow_peak_dict = peak_calling(
    macs_path="macs2",
    bed_paths=bed_paths,
    outdir=os.path.join(peak_dir, "MACS"),
    genome_size="hs",
    n_cpu=4,
    input_format="BEDPE",
    shift=73,
    ext_size=146,
    keep_dup="all",
    q_value=0.05,
    _temp_dir=TMP_DIR)

consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict=narrow_peak_dict,
    peak_half_width=250,
    chromsizes=chromsizes,
    path_to_blacklist=BLACKLIST)
regions_bed = os.path.join(peak_dir, "consensus_regions.bed")
consensus_peaks.to_bed(path=regions_bed, keep=True, compression="infer", chain=False)

# ── 2. QC per sample with the pycistopic CLI (commands as in pycistopic_qc_commands.txt) ──────────
qc_dir = os.path.join(OUT_DIR, "qc")
os.makedirs(qc_dir, exist_ok=True)
tss_bed = os.path.join(qc_dir, "tss.bed")
# TSS annotation as in the pycisTopic tutorial (Ensembl human genes, UCSC hg38 chromosome names)
subprocess.run(["pycistopic", "tss", "get_tss", "--output", tss_bed, "--name", "hsapiens_gene_ensembl",
                "--to-chrom-source", "ucsc", "--ucsc", "hg38"], check=True)
for sample, fragments in fragments_dict.items():
    subprocess.run(["pycistopic", "qc", "--fragments", fragments, "--regions", regions_bed, "--tss", tss_bed,
                    "--output", os.path.join(qc_dir, sample), "--threads", "1"], check=True)

# ── 3. cisTopic object of the cells passing the automatic (Otsu) QC thresholds ────────────────────
import polars as pl
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge

cistopic_obj_list = []
for sample in SAMPLES:
    barcodes, thresholds = get_barcodes_passing_qc_for_sample(
        sample_id=sample, pycistopic_qc_output_dir=qc_dir,
        unique_fragments_threshold=None, tss_enrichment_threshold=None, frip_threshold=0,
        use_automatic_thresholds=True)
    metrics = pl.read_parquet(os.path.join(qc_dir, f"{sample}.fragments_stats_per_cb.parquet")) \
        .to_pandas().set_index("CB").loc[barcodes]
    cistopic_obj_list.append(create_cistopic_object_from_fragments(
        path_to_fragments=fragments_dict[sample], path_to_regions=regions_bed, path_to_blacklist=BLACKLIST,
        metrics=metrics, valid_bc=barcodes, n_cpu=1, project=sample, split_pattern="-"))
cistopic_obj = merge(cistopic_obj_list, project="merged_test")
cistopic_obj.add_cell_data(cell_data, split_pattern="-")
with open(os.path.join(OUT_DIR, "cistopic_obj_1.pkl"), "wb") as f:
    pickle.dump(cistopic_obj, f)

# ── 4. LDA topic models (Mallet) and model selection ─────────────────────────────────────────────
from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models

models = run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 10, 25, 50],
    n_cpu=12,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=TMP_DIR,
    save_path=TMP_DIR,
    mallet_path=MALLET_PATH)
with open(os.path.join(OUT_DIR, "models.pkl"), "wb") as f:
    pickle.dump(models, f)

model = evaluate_models(models, select_model=25, return_model=True)
cistopic_obj.add_LDA_model(model)
with open(os.path.join(OUT_DIR, "cistopic_obj_2.pkl"), "wb") as f:
    pickle.dump(cistopic_obj, f)

# ── 5. Topic region sets: top 3,000 regions and Otsu-binarized topics (from the session history) ─
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates

region_bin_topics_top_3k = binarize_topics(cistopic_obj, method="ntop", ntop=3000, plot=False, num_columns=5)
region_bin_topics_otsu = binarize_topics(cistopic_obj, method="otsu", plot=False, num_columns=5)

for name, topics in [("Topics_otsu", region_bin_topics_otsu), ("Topics_top_3k", region_bin_topics_top_3k)]:
    os.makedirs(os.path.join(OUT_DIR, "region_sets", name), exist_ok=True)
    for topic in topics:
        region_names_to_coordinates(topics[topic].index) \
            .sort_values(["Chromosome", "Start", "End"]) \
            .to_csv(os.path.join(OUT_DIR, "region_sets", name, f"{topic}.bed"),
                    sep="\t", header=False, index=False)
