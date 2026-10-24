#!/usr/bin/env python3
"""
src/network_inference/pycistopic/02_run_dars.py

Differentially accessible regions (DARs) per cell type from the pycisTopic object.
Same commands as /home/be/trinkya/projects/2024_09_18_ccRCC_YU/network/with_fragments_downsampling/run_dars_polars.py
(run in a pycisTopic "polars" environment, see https://github.com/aertslab/pycisTopic/issues/195).

  - highly variable regions on imputed accessibility (min_disp 0.05, min_mean 0.0125, max_mean 3, 20 bins)
  - DARs: find_diff_features on "predicted.id", adjpval_thr = 0.05, log2fc_thr = log2(1.5)
  - region sets written to <out>/region_sets/DARs_cell_type/<cell_type>.bed

Input:  results/network_inference/pycistopic/cistopic_obj_2.pkl (01_run_pycistopic.py)
Output: results/network_inference/pycistopic/{var_features.csv, DARs.pkl, region_sets/DARs_cell_type/}
"""

import os
import pickle

import numpy as np
import pandas as pd
from pycisTopic.diff_features import (
    calculate_per_region_mean_and_dispersion_on_normalized_imputed_acc,
    find_highly_variable_features, impute_accessibility, find_diff_features)
from pycisTopic.utils import region_names_to_coordinates

OUT_DIR = "results/network_inference/pycistopic"

with open(os.path.join(OUT_DIR, "cistopic_obj_2.pkl"), "rb") as f:
    cistopic_obj = pickle.load(f)

topic_region = cistopic_obj.selected_model.topic_region.to_numpy(dtype=np.float32)
cell_topic = cistopic_obj.selected_model.cell_topic.to_numpy(dtype=np.float32)
regions = cistopic_obj.selected_model.topic_region.index.tolist()

(region_names_to_keep, per_region_means, per_region_dispersions) = \
    calculate_per_region_mean_and_dispersion_on_normalized_imputed_acc(
        region_topic=topic_region, cell_topic=cell_topic, region_names=regions,
        scale_factor1=10**6, scale_factor2=10**4, regions_chunk_size=20000)

var_features = find_highly_variable_features(
    features=region_names_to_keep,
    per_region_means_on_normalized_imputed_acc=per_region_means,
    per_region_dispersions_on_normalized_imputed_acc=per_region_dispersions,
    min_disp=0.05, min_mean=0.0125, max_disp=float("inf"), max_mean=3, n_bins=20,
    n_top_features=None, plot=False)
pd.DataFrame(var_features, columns=["regions"]).to_csv(os.path.join(OUT_DIR, "var_features.csv"), index=False)

imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=var_features,
                                       scale_factor=10**6)

markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable="predicted.id",
                                  var_features=var_features, contrasts=None,
                                  adjpval_thr=0.05, log2fc_thr=np.log2(1.5), n_cpu=1)
with open(os.path.join(OUT_DIR, "DARs.pkl"), "wb") as f:
    pickle.dump(markers_dict, f)

os.makedirs(os.path.join(OUT_DIR, "region_sets", "DARs_cell_type"), exist_ok=True)
for cell_type in markers_dict:
    region_names_to_coordinates(markers_dict[cell_type].index) \
        .sort_values(["Chromosome", "Start", "End"]) \
        .to_csv(os.path.join(OUT_DIR, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
                sep="\t", header=False, index=False)
