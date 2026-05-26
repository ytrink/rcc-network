# Systems-level modeling of cellular states in the ccRCC tumor microenvironment through gene regulatory networks


## Overview
This study maps the gene regulatory networks governing cell states across the tumor 
microenvironment (TME) of clear cell renal cell carcinoma (ccRCC). Using single-cell 
multiomic data, we constructed a global gene regulatory network and identified key 
transcription factors predicted to drive clinically relevant cell-state transitions — 
including monocyte-to-TAM differentiation, CD8⁺ T cell exhaustion, and normal-to-tumor 
endothelial cell transitions. This repository contains all code to reproduce the 
analyses and figures presented in the paper.

---

## Data Availability

| Dataset | Source | Accession |
|---------|--------|-----------|
| Yu et al. | GEO | [GSE207493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207493) |
| Li et al. | Mendeley Data | [g67bkbnhhg](https://data.mendeley.com/datasets/g67bkbnhhg/1) |
| Zvirblyte et al. | GEO | [GSE242299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242299) |

Download and place files in the `data/` directory before running any scripts.


---

## Requirements

### R (v4.4)
| Package | Version |
|---------|---------|
| Seurat | 5.2.1 |
| Signac | 1.14.0 |
| igraph | 2.1.4 |

### Python (v3.x)
| Package | Version |
|---------|---------|
| SCENIC+ | 1.0a2 |
| igraph | 2.1.4 |


---
