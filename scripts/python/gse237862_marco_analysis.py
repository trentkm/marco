"""
GSE237862: Colonic muscularis macrophage analysis for MARCO expression.
Stavely et al., Inflamm Bowel Dis, 2025.

Dataset: 4 samples of colonic muscularis propria (10x Chromium scRNA-seq)
  - Naive male, Naive female, DSS male, DSS female

Goal: Load raw Cell Ranger output, QC, normalize, cluster, annotate cell types
using known markers, then extract Marco expression by cell type and condition.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import gzip
from pathlib import Path
from scipy.stats import mannwhitneyu

sc.settings.verbosity = 1

BASE = Path("/Volumes/repos/uc")
DATA = BASE / "data" / "gse237862"
ANALYSIS = BASE / "analysis"

# Sample metadata
SAMPLES = {
    "GSM7655535_RS2": {"sex": "male", "treatment": "naive"},
    "GSM7655536_RS3": {"sex": "male", "treatment": "DSS"},
    "GSM7655537_RS5": {"sex": "female", "treatment": "naive"},
    "GSM7655538_RS6": {"sex": "female", "treatment": "DSS"},
}

# ---- 1. Load and merge all samples ----
print("=" * 60)
print("Loading samples...")
print("=" * 60)

adatas = []
for sample_id, info in SAMPLES.items():
    prefix = DATA / sample_id
    mtx = sio.mmread(str(prefix) + "_matrix.mtx.gz")

    with gzip.open(str(prefix) + "_barcodes.tsv.gz", "rt") as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(str(prefix) + "_features.tsv.gz", "rt") as f:
        features = [line.strip().split("\t") for line in f]

    gene_ids = [feat[0] for feat in features]
    gene_names = [feat[1] for feat in features]

    adata = sc.AnnData(
        X=mtx.T.tocsr(),  # transpose: cells x genes
        obs=pd.DataFrame(index=[f"{sample_id}_{bc}" for bc in barcodes]),
        var=pd.DataFrame({"gene_ids": gene_ids, "gene_names": gene_names},
                         index=gene_names),
    )
    # Handle duplicate gene names
    adata.var_names_make_unique()

    adata.obs["sample"] = sample_id
    adata.obs["sex"] = info["sex"]
    adata.obs["treatment"] = info["treatment"]

    print(f"  {sample_id} ({info['treatment']}, {info['sex']}): "
          f"{adata.n_obs} cells x {adata.n_vars} genes")
    adatas.append(adata)

# Merge
adata = sc.concat(adatas, join="inner")
print(f"\nMerged: {adata.n_obs} cells x {adata.n_vars} genes")

# ---- 2. QC and filtering ----
print("\n" + "=" * 60)
print("QC filtering...")
print("=" * 60)

# Mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

print(f"  Before filtering: {adata.n_obs} cells")

# Standard QC thresholds
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=6000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
sc.pp.filter_genes(adata, min_cells=3)

print(f"  After filtering: {adata.n_obs} cells x {adata.n_vars} genes")

# Save raw counts for DE later
adata.layers["counts"] = adata.X.copy()

# ---- 3. Normalize, HVG, PCA, neighbors, UMAP, clustering ----
print("\n" + "=" * 60)
print("Normalization and clustering...")
print("=" * 60)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Store normalized data
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50)

# Transfer PCA to main object
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)

print(f"  Found {adata.obs['leiden'].nunique()} clusters")
print(f"  Cluster sizes:")
print(adata.obs["leiden"].value_counts().to_string())

# ---- 4. Cell type annotation using marker genes ----
print("\n" + "=" * 60)
print("Cell type annotation...")
print("=" * 60)

# Key markers from the paper and general knowledge:
# Macrophages: Adgre1 (F4/80), Cd68, Csf1r, Cx3cr1
# Lyve1+ macrophages: Lyve1, Cd163, Mrc1 (CD206), Folr2
# MARCO+ macrophages: Marco, Msr1
# Monocytes: Ly6c2, Ccr2, S100a8, S100a9
# Neutrophils: S100a8, S100a9, Cxcr2, Ly6g
# T cells: Cd3d, Cd3e, Cd4, Cd8a
# B cells: Cd19, Cd79a, Ms4a1
# NK cells: Nkg7, Klrb1c, Ncr1
# DCs: Itgax (CD11c), H2-Ab1 (MHC-II), Flt3
# Mast cells: Kit, Cpa3, Mcpt4
# Enteric neurons: Tubb3, Elavl4, Ret, Phox2b
# Enteric glia: Sox10, Plp1, S100b, Gfap
# Smooth muscle: Acta2, Myh11, Des
# ICC (interstitial cells of Cajal): Kit, Ano1
# Fibroblasts: Col1a1, Col3a1, Pdgfra
# Endothelial: Pecam1 (CD31), Cdh5, Kdr

MARKERS = {
    "Macrophage": ["Adgre1", "Cd68", "Csf1r", "Cx3cr1"],
    "Lyve1_Mac": ["Lyve1", "Cd163", "Mrc1", "Folr2"],
    "Marco": ["Marco"],
    "Monocyte": ["Ly6c2", "Ccr2"],
    "Neutrophil": ["S100a8", "S100a9", "Cxcr2"],
    "T_cell": ["Cd3d", "Cd3e"],
    "B_cell": ["Cd79a", "Ms4a1"],
    "NK": ["Nkg7", "Klrb1c"],
    "DC": ["Itgax", "Flt3"],
    "Mast": ["Kit", "Cpa3"],
    "Neuron": ["Tubb3", "Elavl4", "Ret"],
    "Glia": ["Sox10", "Plp1", "S100b"],
    "Smooth_muscle": ["Acta2", "Myh11"],
    "Fibroblast": ["Col1a1", "Pdgfra"],
    "Endothelial": ["Pecam1", "Cdh5"],
}

# Check which markers are present
print("\nMarker gene availability:")
for cat, genes in MARKERS.items():
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    print(f"  {cat}: {len(present)}/{len(genes)} present"
          + (f" (missing: {missing})" if missing else ""))

# Score each cluster for each cell type
# Use mean expression of marker genes per cluster
cluster_scores = {}
for cluster in adata.obs["leiden"].cat.categories:
    mask = adata.obs["leiden"] == cluster
    cluster_data = adata.raw[mask]
    scores = {}
    for cat, genes in MARKERS.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            vals = np.asarray(cluster_data[:, present].X.mean(axis=0)).flatten()
            scores[cat] = vals.mean()
        else:
            scores[cat] = 0.0
    cluster_scores[cluster] = scores

score_df = pd.DataFrame(cluster_scores).T
score_df.index.name = "cluster"
print("\nCluster marker scores (mean normalized expression):")
print(score_df.round(3).to_string())

# Assign cell types based on highest scoring category
# But use a decision tree for more nuanced annotation
def annotate_cluster(scores, cluster_id, adata, mask):
    """Annotate a cluster based on marker scores and additional logic."""
    # Check for non-immune cell types first (often have very distinct markers)
    if scores.get("Neuron", 0) > 0.5:
        return "Enteric neuron"
    if scores.get("Glia", 0) > 0.5:
        return "Enteric glia"
    if scores.get("Smooth_muscle", 0) > 1.0:
        return "Smooth muscle"
    if scores.get("Fibroblast", 0) > 0.5 and scores.get("Macrophage", 0) < 0.3:
        return "Fibroblast"
    if scores.get("Endothelial", 0) > 0.5 and scores.get("Macrophage", 0) < 0.3:
        return "Endothelial"

    # Immune cells
    if scores.get("Macrophage", 0) > 0.3:
        # Is it Lyve1+ subtype?
        if scores.get("Lyve1_Mac", 0) > 0.3:
            return "Lyve1+ macrophage"
        return "Cx3cr1+ macrophage"
    if scores.get("Monocyte", 0) > 0.3:
        return "Monocyte"
    if scores.get("T_cell", 0) > 0.3:
        return "T cell"
    if scores.get("NK", 0) > 0.3:
        return "NK cell"
    if scores.get("B_cell", 0) > 0.3:
        return "B cell"
    if scores.get("DC", 0) > 0.2:
        return "DC"
    if scores.get("Mast", 0) > 0.3:
        return "Mast cell"
    if scores.get("Neutrophil", 0) > 0.5 and scores.get("Macrophage", 0) < 0.2:
        return "Neutrophil"

    # Fallback: highest score
    top = max(scores, key=scores.get)
    return f"Unknown ({top})"


cell_type_map = {}
for cluster in adata.obs["leiden"].cat.categories:
    mask = adata.obs["leiden"] == cluster
    scores = cluster_scores[cluster]
    ct = annotate_cluster(scores, cluster, adata, mask)
    cell_type_map[cluster] = ct
    n = mask.sum()
    print(f"  Cluster {cluster} (n={n}): {ct}")

adata.obs["cell_type"] = adata.obs["leiden"].map(cell_type_map)

print("\nCell type summary:")
print(adata.obs["cell_type"].value_counts().to_string())

# ---- 5. Marco expression analysis ----
print("\n" + "=" * 60)
print("Marco expression analysis")
print("=" * 60)

if "Marco" in adata.var_names:
    marco_vals = np.asarray(adata.raw[:, "Marco"].X.todense()).flatten()
    adata.obs["marco_expr"] = marco_vals
    adata.obs["marco_positive"] = marco_vals > 0

    n_marco_pos = (marco_vals > 0).sum()
    print(f"  Total cells: {adata.n_obs}")
    print(f"  Marco+ cells: {n_marco_pos} ({n_marco_pos/adata.n_obs*100:.1f}%)")
    print(f"  Mean expression: {marco_vals.mean():.4f}")
    print(f"  Max expression: {marco_vals.max():.2f}")

    # By cell type
    print("\n--- Marco by cell type ---")
    ct_summary = adata.obs.groupby("cell_type").agg(
        n_cells=("marco_expr", "size"),
        mean_expr=("marco_expr", "mean"),
        pct_positive=("marco_positive", lambda x: x.mean() * 100),
        max_expr=("marco_expr", "max"),
    ).sort_values("pct_positive", ascending=False).round(3)
    print(ct_summary.to_string())

    # By treatment
    print("\n--- Marco by treatment (naive vs DSS) ---")
    tx_summary = adata.obs.groupby("treatment").agg(
        n_cells=("marco_expr", "size"),
        mean_expr=("marco_expr", "mean"),
        pct_positive=("marco_positive", lambda x: x.mean() * 100),
        max_expr=("marco_expr", "max"),
    ).round(3)
    print(tx_summary.to_string())

    # By treatment AND cell type
    print("\n--- Marco by treatment x cell type ---")
    txct_summary = adata.obs.groupby(["treatment", "cell_type"]).agg(
        n_cells=("marco_expr", "size"),
        mean_expr=("marco_expr", "mean"),
        pct_positive=("marco_positive", lambda x: x.mean() * 100),
        max_expr=("marco_expr", "max"),
    ).round(3)
    # Only show rows with any Marco+ cells or macrophage types
    show = txct_summary[
        (txct_summary["pct_positive"] > 0) |
        txct_summary.index.get_level_values("cell_type").str.contains("macro|mono|Mac", case=False)
    ].sort_values(["treatment", "pct_positive"], ascending=[True, False])
    print(show.to_string())

    # By sex
    print("\n--- Marco by sex ---")
    sex_summary = adata.obs.groupby("sex").agg(
        n_cells=("marco_expr", "size"),
        mean_expr=("marco_expr", "mean"),
        pct_positive=("marco_positive", lambda x: x.mean() * 100),
        max_expr=("marco_expr", "max"),
    ).round(3)
    print(sex_summary.to_string())

    # By sample
    print("\n--- Marco by sample ---")
    sample_summary = adata.obs.groupby(["sample", "treatment", "sex"]).agg(
        n_cells=("marco_expr", "size"),
        mean_expr=("marco_expr", "mean"),
        pct_positive=("marco_positive", lambda x: x.mean() * 100),
        max_expr=("marco_expr", "max"),
    ).round(3)
    print(sample_summary.to_string())

    # Statistical test: naive vs DSS for Marco expression
    print("\n--- Statistical test: naive vs DSS ---")
    naive_marco = adata.obs.loc[adata.obs["treatment"] == "naive", "marco_expr"]
    dss_marco = adata.obs.loc[adata.obs["treatment"] == "DSS", "marco_expr"]
    stat, pval = mannwhitneyu(naive_marco, dss_marco, alternative="two-sided")
    print(f"  Mann-Whitney U: stat={stat:.0f}, p={pval:.2e}")
    print(f"  Naive mean: {naive_marco.mean():.4f}, DSS mean: {dss_marco.mean():.4f}")

    # Within macrophages only
    mac_mask = adata.obs["cell_type"].str.contains("macrophage|Macrophage", case=False)
    if mac_mask.sum() > 0:
        print("\n--- Statistical test: naive vs DSS (macrophages only) ---")
        naive_mac = adata.obs.loc[mac_mask & (adata.obs["treatment"] == "naive"), "marco_expr"]
        dss_mac = adata.obs.loc[mac_mask & (adata.obs["treatment"] == "DSS"), "marco_expr"]
        if len(naive_mac) > 0 and len(dss_mac) > 0:
            stat, pval = mannwhitneyu(naive_mac, dss_mac, alternative="two-sided")
            print(f"  Mann-Whitney U: stat={stat:.0f}, p={pval:.2e}")
            print(f"  Naive macs: n={len(naive_mac)}, mean={naive_mac.mean():.4f}, "
                  f"pct+={( naive_mac > 0).mean()*100:.1f}%")
            print(f"  DSS macs:   n={len(dss_mac)}, mean={dss_mac.mean():.4f}, "
                  f"pct+={(dss_mac > 0).mean()*100:.1f}%")

    # ---- 6. Co-expression with key markers ----
    print("\n" + "=" * 60)
    print("Co-expression analysis: Marco+ vs Marco- macrophages")
    print("=" * 60)

    mac_cells = adata.raw[mac_mask] if mac_mask.sum() > 0 else None

    if mac_cells is not None and mac_cells.n_obs > 10:
        coexpr_genes = [
            "Marco", "Lyve1", "Cd163", "Mrc1", "Folr2", "Cx3cr1", "Csf1r",
            "Cd68", "Adgre1", "Msr1", "Colec12", "Timd4",
            "Ccr2", "Ly6c2", "S100a8", "S100a9", "Fcn1",
            "Il1b", "Tnf", "Il6", "Cxcl10", "Nos2",
            "Il10", "Tgfb1", "Mertk", "Axl",
            "H2-Ab1", "H2-Aa", "H2-Eb1",
            "C1qa", "C1qb", "C1qc",
        ]
        present_genes = [g for g in coexpr_genes if g in adata.var_names]

        marco_mac = adata.obs.loc[mac_mask, "marco_positive"]
        mac_pos_mask = marco_mac.values
        mac_neg_mask = ~marco_mac.values

        n_pos = mac_pos_mask.sum()
        n_neg = mac_neg_mask.sum()
        print(f"  Marco+ macrophages: {n_pos}")
        print(f"  Marco- macrophages: {n_neg}")

        if n_pos >= 5:
            print(f"\n  {'Gene':<12} {'Marco+ mean':>12} {'Marco- mean':>12} {'logFC':>8} {'p-value':>12}")
            print("  " + "-" * 58)

            for gene in present_genes:
                vals = np.asarray(mac_cells[:, gene].X.todense()).flatten()
                pos_vals = vals[mac_pos_mask]
                neg_vals = vals[mac_neg_mask]

                pos_mean = pos_vals.mean()
                neg_mean = neg_vals.mean()

                if neg_mean > 0:
                    logfc = np.log2(max(pos_mean, 0.001) / neg_mean)
                elif pos_mean > 0:
                    logfc = float("inf")
                else:
                    logfc = 0.0

                if pos_vals.std() > 0 or neg_vals.std() > 0:
                    _, pval = mannwhitneyu(pos_vals, neg_vals, alternative="two-sided")
                else:
                    pval = 1.0

                sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
                print(f"  {gene:<12} {pos_mean:>12.3f} {neg_mean:>12.3f} {logfc:>8.2f} {pval:>10.2e} {sig}")

    # ---- 7. Lyve1+ macrophage fate during DSS ----
    print("\n" + "=" * 60)
    print("Lyve1+ macrophage fate during DSS colitis")
    print("=" * 60)

    lyve1_mask = adata.obs["cell_type"] == "Lyve1+ macrophage"
    if lyve1_mask.sum() > 0:
        lyve1_by_tx = adata.obs[lyve1_mask].groupby("treatment").agg(
            n_cells=("marco_expr", "size"),
            mean_marco=("marco_expr", "mean"),
            pct_marco_pos=("marco_positive", lambda x: x.mean() * 100),
        ).round(3)
        print(lyve1_by_tx.to_string())

        # Check Lyve1 expression across all macrophages by treatment
        if "Lyve1" in adata.var_names:
            print("\n  Lyve1 expression in all macrophages by treatment:")
            lyve1_vals = np.asarray(adata.raw[:, "Lyve1"].X.todense()).flatten()
            adata.obs["lyve1_expr"] = lyve1_vals
            mac_lyve1 = adata.obs[mac_mask].groupby("treatment").agg(
                n_cells=("lyve1_expr", "size"),
                mean_lyve1=("lyve1_expr", "mean"),
                pct_lyve1_pos=("lyve1_expr", lambda x: (x > 0).mean() * 100),
            ).round(3)
            print(mac_lyve1.to_string())
    else:
        print("  No Lyve1+ macrophage cluster identified")
        # Still check Lyve1 expression in macrophages
        if "Lyve1" in adata.var_names and mac_mask.sum() > 0:
            print("\n  Lyve1 expression in all macrophages by treatment:")
            lyve1_vals = np.asarray(adata.raw[:, "Lyve1"].X.todense()).flatten()
            adata.obs["lyve1_expr"] = lyve1_vals
            mac_lyve1 = adata.obs[mac_mask].groupby("treatment").agg(
                n_cells=("lyve1_expr", "size"),
                mean_lyve1=("lyve1_expr", "mean"),
                pct_lyve1_pos=("lyve1_expr", lambda x: (x > 0).mean() * 100),
            ).round(3)
            print(mac_lyve1.to_string())

    # ---- 8. DE: Marco+ vs Marco- in muscularis macrophages ----
    print("\n" + "=" * 60)
    print("DE analysis: Marco+ vs Marco- muscularis macrophages")
    print("=" * 60)

    if mac_mask.sum() > 0 and (adata.obs.loc[mac_mask, "marco_positive"]).sum() >= 10:
        adata_mac = adata[mac_mask].copy()
        adata_mac.obs["marco_group"] = adata_mac.obs["marco_positive"].map(
            {True: "Marco+", False: "Marco-"}
        )

        sc.tl.rank_genes_groups(adata_mac, groupby="marco_group", groups=["Marco+"],
                                reference="Marco-", method="wilcoxon",
                                layer="counts", use_raw=False)

        de_results = sc.get.rank_genes_groups_df(adata_mac, group="Marco+")
        de_sig = de_results[
            (de_results["pvals_adj"] < 0.05) &
            (de_results["logfoldchanges"].abs() > 0.25)
        ].copy()

        n_up = (de_sig["logfoldchanges"] > 0).sum()
        n_down = (de_sig["logfoldchanges"] < 0).sum()
        print(f"  Significant DE genes: {len(de_sig)} ({n_up} up, {n_down} down)")

        print(f"\n  Top 30 upregulated in Marco+:")
        top_up = de_sig[de_sig["logfoldchanges"] > 0].head(30)
        print(top_up[["names", "logfoldchanges", "pvals_adj", "scores"]].to_string(index=False))

        print(f"\n  Top 30 downregulated in Marco+:")
        top_down = de_sig[de_sig["logfoldchanges"] < 0].sort_values("logfoldchanges").head(30)
        print(top_down[["names", "logfoldchanges", "pvals_adj", "scores"]].to_string(index=False))

        # Save full DE results
        de_results.to_csv(ANALYSIS / "gse237862_marco_de_results.csv", index=False)
        print(f"\n  Full DE results saved to analysis/gse237862_marco_de_results.csv")

    # ---- 9. DE: Naive vs DSS within macrophages ----
    print("\n" + "=" * 60)
    print("DE analysis: Naive vs DSS in muscularis macrophages")
    print("=" * 60)

    if mac_mask.sum() > 0:
        adata_mac2 = adata[mac_mask].copy()

        n_naive = (adata_mac2.obs["treatment"] == "naive").sum()
        n_dss = (adata_mac2.obs["treatment"] == "DSS").sum()
        print(f"  Naive macrophages: {n_naive}")
        print(f"  DSS macrophages: {n_dss}")

        if n_naive >= 10 and n_dss >= 10:
            sc.tl.rank_genes_groups(adata_mac2, groupby="treatment", groups=["DSS"],
                                    reference="naive", method="wilcoxon",
                                    layer="counts", use_raw=False)

            de_tx = sc.get.rank_genes_groups_df(adata_mac2, group="DSS")
            de_tx_sig = de_tx[
                (de_tx["pvals_adj"] < 0.05) &
                (de_tx["logfoldchanges"].abs() > 0.25)
            ]

            n_up = (de_tx_sig["logfoldchanges"] > 0).sum()
            n_down = (de_tx_sig["logfoldchanges"] < 0).sum()
            print(f"  Significant DE genes: {len(de_tx_sig)} ({n_up} up in DSS, {n_down} down)")

            # Check if Marco itself changes
            marco_de = de_tx[de_tx["names"] == "Marco"]
            if len(marco_de) > 0:
                print(f"\n  Marco in DSS vs naive:")
                print(f"    logFC = {marco_de.iloc[0]['logfoldchanges']:.3f}, "
                      f"padj = {marco_de.iloc[0]['pvals_adj']:.2e}")

            print(f"\n  Top 20 upregulated in DSS macrophages:")
            top_up_tx = de_tx_sig[de_tx_sig["logfoldchanges"] > 0].head(20)
            print(top_up_tx[["names", "logfoldchanges", "pvals_adj"]].to_string(index=False))

            print(f"\n  Top 20 downregulated in DSS macrophages:")
            top_down_tx = de_tx_sig[de_tx_sig["logfoldchanges"] < 0].sort_values(
                "logfoldchanges").head(20)
            print(top_down_tx[["names", "logfoldchanges", "pvals_adj"]].to_string(index=False))

            # Save
            de_tx.to_csv(ANALYSIS / "gse237862_dss_vs_naive_mac_de.csv", index=False)
            print(f"\n  Saved to analysis/gse237862_dss_vs_naive_mac_de.csv")

else:
    print("  Marco not found in gene list!")

# ---- 10. Summary output ----
print("\n" + "=" * 60)
print("Saving summary CSV")
print("=" * 60)

summary_rows = []
for (tx, ct), grp in adata.obs.groupby(["treatment", "cell_type"]):
    summary_rows.append({
        "dataset": "GSE237862",
        "treatment": tx,
        "cell_type": ct,
        "n_cells": len(grp),
        "mean_marco": grp["marco_expr"].mean(),
        "pct_marco_positive": grp["marco_positive"].mean() * 100,
        "max_marco": grp["marco_expr"].max(),
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(ANALYSIS / "gse237862_marco_summary.csv", index=False)
print(f"  Saved to analysis/gse237862_marco_summary.csv")
print(summary_df.sort_values("pct_marco_positive", ascending=False).round(3).to_string(index=False))

print("\n" + "=" * 60)
print("DONE")
print("=" * 60)
