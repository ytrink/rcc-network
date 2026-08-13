# Gene regulatory network analysis identifies transcription factors governing major cell-state transitions in the ccRCC tumor microenvironment

This repository contains reproducible scripts to analyze single-cell RNA-seq, single-cell ATAC-seq, and SCENIC+ gene regulatory networks (GRNs) for clear cell renal cell carcinoma (ccRCC).

---

## 📌 Repository Structure & Workflow

```
.
├── preprocessing/
│   ├── scrna/               # scRNA-seq QC, SCTransform normalization, Harmony integration & inferCNV
│   │   ├── yu/              # Yu et al. cohort processing (01_qc, 02_integration, 03_annotation, 04_infercnv)
│   │   ├── li/              # Li et al. cohort processing (01_convert_h5ad, 02_processing)
│   │   └── zvirblyte/       # Zvirblyte et al. cohort processing (01_processing)
│   ├── scatac/
│   │   └── yu/              # scATAC-seq peak filtering, Harmony integration & Signac label transfer
│   └── monocle/             # Monocle3 trajectory inference & cell embedding (01_mac_mono, 02_cd8_tcells)
│
├── network_inference/
│   ├── config/config.yaml   # SCENIC+ Snakemake configuration file
│   └── workflow/Snakefile   # SCENIC+ GRN inference pipeline
│
├── figure1/                 # AUCell regulon scoring, ComplexHeatmap & Leiden GRN community clustering
├── figure2/                 # Macrophage/Monocyte regulatory influence & in-silico TF perturbations
├── figure3/                 # CD8+ T-cell exhaustion regulatory influence & in-silico TF perturbations
├── figure4/                 # Endothelial / TEC regulatory influence & in-silico TF perturbations
├── supplementary/           # Centrality analysis (Reverse PageRank & Undirected Betweenness)
└── utils/                   # Shared helper functions & eRegulon to BEDPE genome track converter
```

---

## 🚀 How to Run

### 1. Preprocessing & Trajectories
Run cell-lineage processing and Monocle3 trajectory embeddings:
```bash
# scRNA-seq & scATAC-seq processing
Rscript preprocessing/scrna/yu/01_qc_yu.R
Rscript preprocessing/scrna/yu/02_integration_yu.R
Rscript preprocessing/scrna/yu/03_annotation_yu.R
Rscript preprocessing/scrna/yu/04_infercnv_yu.R

# Monocle3 embeddings
Rscript preprocessing/monocle/01_run_monocle3_mac_mono.R
Rscript preprocessing/monocle/02_run_monocle3_cd8_tcells.R
```

### 2. Network Inference (pycisTopic & SCENIC+)
Generate candidate regulatory regions and infer the GRN:
```bash
# pycisTopic peak calling, LDA topic modeling, and DAR extraction
python network_inference/pycistopic/01_run_pseudobulk_macs2_peaks.py
python network_inference/pycistopic/02_run_pycistopic_lda.py
python network_inference/pycistopic/03_run_dars_polars.py

# Run SCENIC+ GRN Snakemake workflow
snakemake --snakefile network_inference/workflow/Snakefile --cores 12
```

### 3. Downstream Figures & Perturbations
Reproduce figure panels and in-silico TF perturbation simulations:
```bash
# Figure 1: Regulon activity AUCell scoring & Leiden GRN community clustering
Rscript figure1/R/01_aucell_scores.R
Rscript figure1/R/02_aucell_heatmap.R
python figure1/python/01_grn_community_clustering.py

# Figures 2, 3, & 4: Regulatory influence & in-silico TF perturbations
Rscript figure2/R/02_regulatory_influence.R
python figure2/python/02_run_perturbation.py
python figure2/python/04_plot_arrow_umap.py

# Supplementary: Reverse PageRank & Betweenness Centrality
Rscript supplementary/compute_pagerank_and_betweeness_centrality.R
UMAP and genome plots can be reproduced using main figure scripts - please contact the lead author with any questions.
```

---

## 💾 Data Availability Note

> [!NOTE]
> Large data files, such as the SCENIC+ output files, are available upon request from the corresponding author.

