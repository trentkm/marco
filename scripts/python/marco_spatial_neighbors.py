"""
Spatial neighborhood analysis of Marco+ spots in SCP2771 (Visium DSS colitis).
Finds nearest neighbors of Marco+ spots and characterizes their gene expression.
"""

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from pathlib import Path

BASE = Path("/Volumes/repos/uc")
DATA = BASE / "data" / "single-cell-portal" / "SCP2771"
RESULTS = BASE / "analysis"

# ---- Load data ----
print("Loading SCP2771 data...")

# Expression
expr_dir = DATA / "expression" / "66e1e0765e4f1f173a036b83"
with open(expr_dir / "genes.tsv") as f:
    genes = [line.strip().split("\t")[0] for line in f]
with open(expr_dir / "barcodes.tsv") as f:
    barcodes = [line.strip().split("\t")[0] for line in f]
mat = mmread(expr_dir / "matrix.mtx").tocsc()

# Metadata
meta = pd.read_csv(DATA / "metadata" / "dss_spatial_meta.tsv", sep="\t", index_col=0)
if meta.index[0] == "TYPE":
    meta = meta.iloc[1:]

# Spatial coordinates
spatial = pd.read_csv(DATA / "cluster" / "dss_spatial_cluster.tsv", sep="\t", index_col=0)
if spatial.index[0] == "TYPE":
    spatial = spatial.iloc[1:]
spatial["X"] = spatial["X"].astype(float)
spatial["Y"] = spatial["Y"].astype(float)

# Build expression DataFrame (sparse → dense only for Marco and key markers)
marco_idx = genes.index("Marco")
marco_vals = np.asarray(mat.getrow(marco_idx).toarray().flatten())
expr_df = pd.DataFrame({"Marco": marco_vals}, index=barcodes)

# Merge everything
common = expr_df.index.intersection(meta.index).intersection(spatial.index)
df = pd.DataFrame(index=common)
df["Marco"] = expr_df.loc[common, "Marco"].values
df["X"] = spatial.loc[common, "X"].values
df["Y"] = spatial.loc[common, "Y"].values
df["Category"] = spatial.loc[common, "Category"].values
df["dss_time_point"] = meta.loc[common, "dss_time_point"].values

print(f"Total spots with spatial + expression + metadata: {len(df)}")
print(f"Marco+ spots: {(df['Marco'] > 0).sum()}")

# ---- Focus on D12 (peak inflammation) ----
d12 = df[df["dss_time_point"] == "D12"].copy()
print(f"\nD12 spots: {len(d12)}")
print(f"D12 Marco+ spots: {(d12['Marco'] > 0).sum()}")

# ---- Build spatial KD-tree for D12 ----
coords = d12[["X", "Y"]].values
tree = cKDTree(coords)

# Find nearest neighbors for each Marco+ spot
marco_pos_mask = d12["Marco"] > 0
marco_pos_idx = np.where(marco_pos_mask.values)[0]
marco_neg_idx = np.where(~marco_pos_mask.values)[0]

# Query k nearest neighbors (excluding self)
K = 6  # 6 neighbors (Visium hexagonal grid has 6 immediate neighbors)
distances, neighbor_indices = tree.query(coords[marco_pos_idx], k=K+1)
# Remove self (first neighbor)
neighbor_indices = neighbor_indices[:, 1:]
neighbor_distances = distances[:, 1:]

print(f"\nNearest neighbor analysis (K={K}):")
print(f"  Median neighbor distance: {np.median(neighbor_distances[:, 0]):.3f}")
print(f"  Mean neighbor distance: {np.mean(neighbor_distances[:, 0]):.3f}")

# ---- Analyze cluster categories of neighbors ----
print(f"\n{'='*60}")
print("Spatial cluster composition around Marco+ spots")
print(f"{'='*60}")

d12_barcodes = d12.index.tolist()

# Categories of Marco+ spots
marco_cats = d12.loc[marco_pos_mask, "Category"].value_counts().sort_index()
print(f"\nMarco+ spots by cluster:")
for cat, n in marco_cats.items():
    pct = n / marco_cats.sum() * 100
    print(f"  Cluster {cat}: {n} ({pct:.1f}%)")

# Categories of neighbors of Marco+ spots
neighbor_cats = []
for row in neighbor_indices:
    for idx in row:
        if idx < len(d12_barcodes):
            neighbor_cats.append(d12.iloc[idx]["Category"])

neighbor_cat_counts = pd.Series(neighbor_cats).value_counts().sort_index()
print(f"\nNeighbors of Marco+ spots by cluster:")
for cat, n in neighbor_cat_counts.items():
    pct = n / len(neighbor_cats) * 100
    print(f"  Cluster {cat}: {n} ({pct:.1f}%)")

# Compare to overall D12 cluster distribution
overall_cats = d12["Category"].value_counts().sort_index()
print(f"\nOverall D12 cluster distribution:")
for cat, n in overall_cats.items():
    pct = n / len(d12) * 100
    print(f"  Cluster {cat}: {n} ({pct:.1f}%)")

# Enrichment: which clusters are over/under-represented near Marco+ spots
print(f"\nCluster enrichment near Marco+ spots (observed/expected):")
for cat in sorted(d12["Category"].unique()):
    obs = neighbor_cat_counts.get(cat, 0) / len(neighbor_cats) if len(neighbor_cats) > 0 else 0
    exp = overall_cats.get(cat, 0) / len(d12)
    ratio = obs / exp if exp > 0 else float("inf")
    print(f"  Cluster {cat}: {ratio:.2f}x (obs={obs*100:.1f}%, exp={exp*100:.1f}%)")

# ---- Gene expression in neighbors vs non-neighbors ----
print(f"\n{'='*60}")
print("Gene expression: Marco+ spot neighbors vs rest of D12")
print(f"{'='*60}")

# Get all unique neighbor spot indices (excluding Marco+ spots themselves)
all_neighbor_idx = set()
for row in neighbor_indices:
    for idx in row:
        if idx < len(d12) and idx not in marco_pos_idx:
            all_neighbor_idx.add(idx)

# Also get "distant" spots (not Marco+ and not neighbors)
all_marco_and_neighbor_idx = set(marco_pos_idx) | all_neighbor_idx
distant_idx = [i for i in range(len(d12)) if i not in all_marco_and_neighbor_idx]

print(f"  Marco+ spots: {len(marco_pos_idx)}")
print(f"  Neighbor spots (non-Marco+): {len(all_neighbor_idx)}")
print(f"  Distant spots: {len(distant_idx)}")

# Extract expression for all genes in these groups
# Use sparse matrix efficiently
d12_barcode_positions = [barcodes.index(b) for b in d12.index]

neighbor_positions = [d12_barcode_positions[i] for i in all_neighbor_idx]
distant_positions = [d12_barcode_positions[i] for i in distant_idx]
marco_positions = [d12_barcode_positions[i] for i in marco_pos_idx]

# For each gene, compare neighbor expression vs distant expression
print("\nComputing DE: neighbors of Marco+ vs distant spots...")
results = []
for g_idx, gene in enumerate(genes):
    row = mat.getrow(g_idx).toarray().flatten()

    neighbor_vals = row[neighbor_positions]
    distant_vals = row[distant_positions]
    marco_vals_g = row[marco_positions]

    mean_neighbor = np.mean(neighbor_vals)
    mean_distant = np.mean(distant_vals)
    mean_marco = np.mean(marco_vals_g)

    pct_neighbor = (neighbor_vals > 0).mean() * 100
    pct_distant = (distant_vals > 0).mean() * 100
    pct_marco = (marco_vals_g > 0).mean() * 100

    # Only test genes expressed in >1% of either group
    if pct_neighbor > 1 or pct_distant > 1:
        try:
            stat, pval = mannwhitneyu(neighbor_vals, distant_vals, alternative="two-sided")
        except:
            pval = 1.0

        if mean_distant > 0:
            log2fc = np.log2((mean_neighbor + 0.01) / (mean_distant + 0.01))
        else:
            log2fc = np.log2(mean_neighbor + 0.01) - np.log2(0.01)

        results.append({
            "gene": gene,
            "mean_marco_spots": round(mean_marco, 4),
            "pct_marco_spots": round(pct_marco, 1),
            "mean_neighbors": round(mean_neighbor, 4),
            "pct_neighbors": round(pct_neighbor, 1),
            "mean_distant": round(mean_distant, 4),
            "pct_distant": round(pct_distant, 1),
            "log2fc_neighbor_vs_distant": round(log2fc, 3),
            "pval": pval,
        })

results_df = pd.DataFrame(results)

# Multiple testing correction (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
reject, padj, _, _ = multipletests(results_df["pval"], method="fdr_bh")
results_df["padj"] = padj
results_df = results_df.sort_values("log2fc_neighbor_vs_distant", ascending=False)

# Save full results
results_df.to_csv(RESULTS / "scp2771_spatial_neighbor_de.csv", index=False)

# Print top upregulated in neighbors
sig_up = results_df[(results_df["padj"] < 0.05) & (results_df["log2fc_neighbor_vs_distant"] > 0.5)]
sig_down = results_df[(results_df["padj"] < 0.05) & (results_df["log2fc_neighbor_vs_distant"] < -0.5)]

print(f"\nSignificant genes enriched near Marco+ spots: {len(sig_up)}")
print(f"Significant genes depleted near Marco+ spots: {len(sig_down)}")

print(f"\nTop 30 genes ENRICHED in Marco+ spot neighbors (vs distant):")
print(sig_up.head(30)[["gene", "mean_marco_spots", "pct_marco_spots", "mean_neighbors", "pct_neighbors", "mean_distant", "pct_distant", "log2fc_neighbor_vs_distant", "padj"]].to_string(index=False))

print(f"\nTop 30 genes DEPLETED near Marco+ spots:")
print(sig_down.sort_values("log2fc_neighbor_vs_distant").head(30)[["gene", "mean_marco_spots", "pct_marco_spots", "mean_neighbors", "pct_neighbors", "mean_distant", "pct_distant", "log2fc_neighbor_vs_distant", "padj"]].to_string(index=False))

# ---- Key marker categories in the spatial neighborhood ----
print(f"\n{'='*60}")
print("Specific marker genes in Marco+ spots vs neighbors vs distant")
print(f"{'='*60}")

marker_categories = {
    "Immune / Myeloid": ["Cd68", "Cd14", "Csf1r", "Itgam", "Adgre1", "Cx3cr1", "Ccr2", "Ly6c2", "Lyz2"],
    "T cells": ["Cd3e", "Cd3d", "Cd4", "Cd8a", "Foxp3", "Tbx21", "Gata3", "Rorc", "Il17a"],
    "B cells": ["Cd19", "Cd79a", "Ms4a1", "Ighm"],
    "Neutrophils": ["S100a8", "S100a9", "Ly6g", "Cxcr2"],
    "Fibroblasts / Stroma": ["Col1a1", "Col1a2", "Col3a1", "Vim", "Acta2", "Des", "Pdgfra"],
    "Epithelial": ["Epcam", "Krt20", "Muc2", "Lgr5", "Cdx2"],
    "Endothelial": ["Pecam1", "Cdh5", "Vwf"],
    "Tissue remodeling": ["Mmp2", "Mmp9", "Mmp10", "Mmp13", "Timp1", "Tgfb1"],
    "Inflammatory cytokines": ["Il1b", "Il6", "Tnf", "Ifng", "Il10", "Il33"],
    "Scavenger receptors": ["Marco", "Msr1", "Cd163", "Mrc1", "Cd36", "Scarb1", "Colec12"],
    "Yolk sac / tissue-resident": ["Lyve1", "Folr2", "Timd4", "Vsig4", "Mertk", "Cx3cr1"],
}

for category, markers in marker_categories.items():
    print(f"\n--- {category} ---")
    print(f"  {'Gene':12s} {'Marco+ spots':>14s} {'Neighbors':>14s} {'Distant':>14s} {'Neighbor/Distant':>16s}")
    for gene in markers:
        if gene in genes:
            g_idx = genes.index(gene)
            row = mat.getrow(g_idx).toarray().flatten()

            m_vals = row[marco_positions]
            n_vals = row[neighbor_positions]
            d_vals = row[distant_positions]

            m_pct = f"{(m_vals>0).mean()*100:.1f}%"
            n_pct = f"{(n_vals>0).mean()*100:.1f}%"
            d_pct = f"{(d_vals>0).mean()*100:.1f}%"

            n_mean = np.mean(n_vals)
            d_mean = np.mean(d_vals)
            ratio = n_mean / d_mean if d_mean > 0 else float("inf")

            print(f"  {gene:12s} {m_pct:>14s} {n_pct:>14s} {d_pct:>14s} {ratio:>14.2f}x")
        else:
            print(f"  {gene:12s} {'not found':>14s}")

# ---- Also analyze all timepoints for spatial context ----
print(f"\n{'='*60}")
print("Marco+ spatial neighborhood by timepoint")
print(f"{'='*60}")

for tp in ["D0", "D12", "D30", "D73"]:
    tp_data = df[df["dss_time_point"] == tp]
    n_marco = (tp_data["Marco"] > 0).sum()
    if n_marco == 0:
        print(f"\n{tp}: No Marco+ spots")
        continue

    tp_coords = tp_data[["X", "Y"]].values
    tp_tree = cKDTree(tp_coords)
    tp_marco_idx = np.where((tp_data["Marco"] > 0).values)[0]

    # Get neighbor clusters
    tp_distances, tp_neighbors = tp_tree.query(tp_coords[tp_marco_idx], k=K+1)
    tp_neighbors = tp_neighbors[:, 1:]

    tp_neighbor_cats = []
    for row in tp_neighbors:
        for idx in row:
            if idx < len(tp_data):
                tp_neighbor_cats.append(tp_data.iloc[idx]["Category"])

    print(f"\n{tp}: {n_marco} Marco+ spots")
    print(f"  Marco+ cluster distribution: {dict(tp_data.loc[tp_data['Marco']>0, 'Category'].value_counts().sort_index())}")
    if tp_neighbor_cats:
        nc = pd.Series(tp_neighbor_cats).value_counts().sort_index()
        print(f"  Neighbor cluster distribution: {dict(nc)}")

print(f"\n{'='*60}")
print("Analysis complete.")
print(f"Full results saved to: {RESULTS / 'scp2771_spatial_neighbor_de.csv'}")
