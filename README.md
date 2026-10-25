# Gene regulatory network analysis identifies transcription factors governing major cell-state transitions in the ccRCC tumor microenvironment

Scripts to reproduce the analyses and figures of the paper: a SCENIC+ gene regulatory network (GRN) inferred from
single-cell RNA-seq and ATAC-seq of 19 ccRCC patients (Yu et al.), regulon activity (AUCell) and network communities
(Leiden), regulatory influence scores, network centrality, and in-silico TF knockdowns in macrophages/monocytes,
CD8⁺ T cells (Li et al.) and endothelial cells (Zvirblyte et al.).

## Data

Processed data are archived on Zenodo: **https://doi.org/10.5281/zenodo.23648072**
Raw data: Yu et al. GEO [GSE207493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207493),
Zvirblyte et al. GEO [GSE242299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242299),
Li et al. [Mendeley Data](https://data.mendeley.com/datasets/g67bkbnhhg/1).

All scripts find the repository root through the `.here` file (R: `here`, Python: `src/utils/paths.py`), read inputs
from `data/` and write outputs to `results/`.

**Getting the data.** Small inputs are in this repository. The Zenodo files are zips whose contents sit under
`data/...`, so they unpack into place when unzipped at the repository root. From the repository root run:

```bash
bash download_data.sh
```

or download the zips from the Zenodo page and unzip them at the repository root.

```
data/
├── scplus/      eRegulon_direct.tsv, eRegulons_extended.tsv, scplus_direct_extended.graphml   [Zenodo: scplus_grn.zip]
│                leiden_RB_03_pos.csv (Table S1)                                                 [repo]
├── li/          Li_mac_mono.rds                                                                 [Zenodo: li_mac_mono_object.zip]
│                metadata_merged.csv                                                             [Zenodo: li_metadata.zip]
│                degs_per_cell_type.csv, {mac_mono,CD8}_expressed_genes_10perc.csv,
│                Li_sampled_broad2_5000_barcodes.csv (cells of Fig. 1C)                          [repo]
│                Li_integrated.rds, Li_CD8T_normalized.rds, Li_T.rds   (not archived: built from the Li et al. data
│                                                                        with src/preprocessing/scrna/li/)
├── zvirblyte/   zvir_endo.rds                                                                   [Zenodo: zvirblyte_endothelial.zip]
│                X_umap_published.csv (Fig. 4A/S12/S15 coordinates)                              [repo]
├── monocle/     cds_mac_mono_umap_sampled30k.csv, cds_tcells_umap.csv                           [Zenodo: monocle_umaps.zip]
├── perturbation/mac_mono/  cds_mac_mono_sampled15k.csv, cds_mac_mono_sampled10k.csv,
│                           cds_mac_mono_umap_sampled10k.csv                                     [Zenodo: perturbation_mac_mono.zip]
├── perturbation/CD8/       cds_normalized_tcells_sampled_7500.csv, cds_tcells_umap_sampled_7500.csv
│                                                                                                [Zenodo: perturbation_CD8.zip]
└── yu/          scrna/, scatac/  (GEO downloads; for preprocessing and network inference only)
```

With these files the figures can be reproduced from step 3 below. Figure 1C (AUCell) and the CD8⁺ T-cell influence
scores (Fig. 3B) additionally need `Li_integrated.rds` / `Li_CD8T_normalized.rds` from step 1.

## Run order

Steps 1–2 regenerate the network from raw data and are computationally heavy. To reproduce the figures, start at
step 3 with the Zenodo files.

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| **1. Preprocessing** | `src/preprocessing/scrna/yu/01_qc_yu.R` … `04_infercnv_yu.R` | Yu scRNA-seq QC (200–6000 features, >1000 UMI, <10% mt), scDblFinder, SCTransform, Harmony, annotation, inferCNV | `data/yu/scrna/` |
| | `src/preprocessing/scatac/yu/01_qc_atac_yu.R` … `03_annotation_transfer_yu.R` | Yu scATAC-seq QC, TF-IDF/LSI, Harmony (components 2–30 downstream), label transfer (rPCA) | `data/yu/scatac/` |
| | `src/preprocessing/scrna/li/01_convert_h5ad_to_mtx.py`, `02_processing_li.R` | Li et al. object | `data/li/Li_integrated.rds` |
| | `src/preprocessing/scrna/li/03_expressed_genes_and_markers.R` | ≥10%-expressed gene lists; per-cell-type marker genes (FindAllMarkers) | `data/li/*.csv` |
| | `src/preprocessing/scrna/zvirblyte/01_processing_zvirblyte.R` | Zvirblyte et al. endothelial object | `data/zvirblyte/zvir_endo.rds` |
| | `src/preprocessing/monocle/01_run_monocle3_mac_mono.R`, `02_run_monocle3_cd8_tcells.R` | Monocle3 UMAPs and the cell subsets used for perturbation | `data/monocle/`, `data/perturbation/` |
| **2. Network inference** | `src/network_inference/pycistopic/01_run_pycistopic.py` | Candidate regions (6 samples, 2 per WHO/ISUP grade): pseudobulk MACS2 peaks, QC, LDA topics. *Consolidated from the interactive pycisTopic session; parameters verified against the saved objects* | `results/network_inference/pycistopic/` |
| | `src/network_inference/pycistopic/02_run_dars.py` | DARs (adj. p < 0.05, log2FC > log2(1.5)) | region sets |
| | `snakemake --snakefile src/network_inference/workflow/Snakefile --configfile src/network_inference/config/config.yaml` | SCENIC+ (v1.0a2, precomputed hg38 cisTarget databases). Edit the absolute paths in `config.yaml` first | eRegulon tables |
| **3. Figure 1** | `src/figure1/python/01_grn_community_clustering.py` | Leiden communities (RB, resolution 0.3, activating links) | `results/figure1/community_clustering/` |
| | `src/figure1/R/01_aucell_scores.R` → `02_aucell_heatmap.R` | AUCell on Li et al. ccRCC patients; clustered heatmap (Fig. 1C) | `results/figure1/` |
| **4. Figure 2** (mac/mono) | `src/figure2/R/01_umap.R`, `02_regulatory_influence.R` | UMAP; pseudobulk DESeq2 (~patient + condition, padj < 0.1, \|log2FC\| > 0.5); influence scores | `results/figure2/` |
| | `src/figure2/python/01_get_regressor.py` → `02` → `03` → `04` | Knockdown simulations (regressors, 5 iterations, iteration plots, UMAP arrows) | `results/figure2/perturbation/` |
| **5. Figure 3** (CD8) | `src/figure3/R/01_umap.R`, `02_regulatory_influence.R`; `src/figure3/python/01` → `04` | Exhaustion score (LAG3+TIGIT, median split), DESeq2 (padj < 0.05, \|log2FC\| > 0.5), knockdowns | `results/figure3/` |
| **6. Figure 4** (endothelial) | `src/figure4/R/01_umap.R` → `02_regulatory_influence.R`; `src/figure4/python/01` → `04` | Re-embedding (SCTransform → Harmony → UMAP) and UMAP panels (Fig. 4A, S12, S15); DESeq2 (padj < 0.05, \|log2FC\| > 0.5); MEF2C targets for ToppGene; knockdowns | `results/figure4/` |
| **7. Supplementary** | `src/supplementary/compute_pagerank_and_betweeness_centrality.R` (after steps 4–6) | Reverse PageRank and betweenness in cell-type subnetworks (Figs. S6, S10, S13; Table S3) | `results/supplementary/centrality/` |
| | `src/supplementary/centrality_edge_weight_sensitivity.R` | Same centralities with the maximum or the sum of RE weights per TF–target edge (robustness check, Methods) | `results/supplementary/centrality/` |
| | `src/supplementary/tnk_subnetwork_for_cytoscape.R` | T/NK subnetwork of Fig. S3 (Table S3 CD8⁺ TFs → targets expressed in ≥10% of CD8⁺ T cells), for Cytoscape | `results/supplementary/tnk_subnetwork/` |
| | `src/utils/scenic_plus_to_bedpe.py` | Region-to-gene links for genome browser views (pyGenomeTracks) | |

## Environments

- R 4.4.1: Seurat 5.2.1, SeuratObject, Signac 1.14.0, harmony 1.2.3, DESeq2, monocle3 1.4.15, AUCell, GSEABase,
  ComplexHeatmap 2.20.0, circlize, igraph 2.1.4, scDblFinder 1.18.0, infercnv 1.20.0, here, ggrepel, data.table.
- Python 3.11: SCENIC+ 1.0a2 (`scenicplus.simulation` for the perturbations), pycisTopic 2.0a0, pandas, numpy,
  networkx 3.3, python-igraph, leidenalg, scanpy, matplotlib.

## Reproducibility notes

- **Stochastic steps.** Leiden clustering (the original run had no seed), UMAP and random cell sampling can differ
  between package versions. The published outputs are provided so the figures can be reproduced exactly:
  `data/scplus/leiden_RB_03_pos.csv` (Table S1), the Monocle UMAP/sample files, `Li_sampled_broad2_5000_barcodes.csv`
  (Fig. 1C cells) and, optionally, `zvirblyte/X_umap_published.csv` (Fig. 4A/S12/S15 coordinates).
- **Perturbation regressors** are trained, per cell type, for marker genes of that cell type (Li et al. FindAllMarkers,
  `degs_per_cell_type.csv`) or, for endothelial cells, TEC-vs-EC DEGs, on raw-count matrices of 7,500–15,000 cells.
  The macrophage/monocyte models of the paper were fit on a separate random sample (14,000 cells) of the same
  population.
- **Centrality** uses the GRN graphml (one edge per TF–target pair). The main results are robust to the choice of edge weight
  (single region-to-gene link, maximum, or sum over regions).
- `degs_per_cell_type.csv`: the original script was not saved; `03_expressed_genes_and_markers.R` uses parameters
  inferred from the published file. Use the published file to reproduce the figures exactly.
