#!/usr/bin/env python3
"""
figures/figure1/python/01_grn_community_clustering.py

Performs Leiden community clustering on the SCENIC+ gene regulatory network (GRN)
using positive regulation links (regulation == 1) and interaction weights (rho_TF2G).

Input:  data/eRegulon_direct.tsv
        data/eRegulons_extended.tsv
Output: data/community_clustering/leiden_communities_pos_res0.3.csv
        data/community_clustering/tf_community_mapping.csv
        data/community_clustering/grn_positive_communities.graphml
"""

import os
import pandas as pd
import networkx as nx
import igraph as ig
import leidenalg as la

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR = os.path.join(ROOT_DIR, "data")
OUT_DIR  = os.path.join(ROOT_DIR, "data", "community_clustering")

os.makedirs(OUT_DIR, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────────────────────
RESOLUTION_PARAM = 0.3  # Resolution parameter matching manuscript Methods text
SEED             = 101  # Reproducibility seed

# ── Load eRegulon data ─────────────────────────────────────────────────────────
print("Loading SCENIC+ eRegulon tables...")
direct_path   = os.path.join(DATA_DIR, "eRegulon_direct.tsv")
extended_path = os.path.join(DATA_DIR, "eRegulons_extended.tsv")

if not os.path.exists(extended_path):
    extended_path = os.path.join(DATA_DIR, "eRegulon_extended.tsv")

direct_eregulons   = pd.read_table(direct_path)
extended_eregulons = pd.read_table(extended_path)

# Filter for positive regulation links only (regulation == 1)
direct_pos   = direct_eregulons[direct_eregulons["regulation"] == 1]
extended_pos = extended_eregulons[extended_eregulons["regulation"] == 1]

# Combine unique TF -> Gene edge pairs with correlation weights (rho_TF2G)
edge_df = pd.concat([
    direct_pos[["TF", "Gene", "rho_TF2G"]].drop_duplicates(),
    extended_pos[["TF", "Gene", "rho_TF2G"]].drop_duplicates()
]).drop_duplicates()

print(f"Total positive regulatory edges: {len(edge_df)}")

# ── Build NetworkX Graph ───────────────────────────────────────────────────────
G_pos = nx.from_pandas_edgelist(
    edge_df,
    source="TF",
    target="Gene",
    edge_attr="rho_TF2G",
    create_using=nx.DiGraph()
)

# Extract list of all TFs present in the network
all_tfs = set(edge_df["TF"].unique())
print(f"Total nodes: {G_pos.number_of_nodes()} (TFs: {len(all_tfs)})")

# ── Convert to igraph for Leiden Clustering ────────────────────────────────────
print("Converting to igraph for community detection...")
igraph_G = ig.Graph(directed=True)
igraph_G.add_vertices(list(G_pos.nodes))
igraph_G.add_edges(list(G_pos.edges))

# Edge weights: absolute value of SCENIC+ rho_TF2G
edge_weights = [abs(G_pos[u][v]["rho_TF2G"]) for u, v in G_pos.edges]
igraph_G.es["rho_TF2G"] = edge_weights
igraph_G.vs["name"]     = list(G_pos.nodes)

# ── Run Leiden Community Detection ─────────────────────────────────────────────
print(f"Running Leiden algorithm (RBConfigurationVertexPartition, resolution = {RESOLUTION_PARAM})...")

partition = la.find_partition(
    igraph_G,
    la.RBConfigurationVertexPartition,
    weights="rho_TF2G",
    resolution_parameter=RESOLUTION_PARAM,
    seed=SEED
)

print(f"Identified {len(partition)} network communities.")

# ── Extract Community Assignments ──────────────────────────────────────────────
community_records = []
for comm_id, community_indices in enumerate(partition):
    community_nodes = [igraph_G.vs[idx]["name"] for idx in community_indices]
    community_tfs   = [node for node in community_nodes if node in all_tfs]
    
    for node in community_nodes:
        is_tf = node in all_tfs
        community_records.append({
            "node": node,
            "is_tf": is_tf,
            "community_id": f"Community_{comm_id}",
            "community_size": len(community_nodes),
            "num_community_tfs": len(community_tfs),
            "community_tfs": ", ".join(community_tfs)
        })

comm_df = pd.DataFrame(community_records)

# ── Save Output Table ──────────────────────────────────────────────────────────
out_csv = os.path.join(OUT_DIR, "leiden_communities_pos_res0.3.csv")
comm_df.to_csv(out_csv, index=False)
print(f"Saved community table -> {out_csv}")

# Export TF-only community mapping
tf_comm_df = comm_df[comm_df["is_tf"]][["node", "community_id"]].rename(columns={"node": "TF"})
tf_csv = os.path.join(OUT_DIR, "tf_community_mapping.csv")
tf_comm_df.to_csv(tf_csv, index=False)
print(f"Saved TF community mapping -> {tf_csv}")

# ── Attach Communities to NetworkX Graph & Save GraphML ───────────────────────
comm_map = comm_df.set_index("node")["community_id"].to_dict()
nx.set_node_attributes(G_pos, comm_map, "community_id")

out_graphml = os.path.join(OUT_DIR, "grn_positive_communities.graphml")
nx.write_graphml(G_pos, out_graphml)
print(f"Saved GraphML network -> {out_graphml}")

print("Community clustering complete!")
