"""
MARCO expression analysis across single-cell datasets.
Checks MARCO/Marco expression by cell type, species, disease state, and spatial location.
"""

import os
import gzip
import numpy as np
import pandas as pd
import scipy.io as sio
import scanpy as sc
from pathlib import Path

BASE = Path("/Volumes/repos/uc")
DATA = BASE / "data"
SCP = DATA / "single-cell-portal"
RESULTS = BASE / "analysis" / "results"
RESULTS.mkdir(exist_ok=True)


def load_genes(path):
    """Load gene names from a tsv file (plain or gzipped)."""
    p = str(path)
    if p.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            genes = [line.strip().split("\t")[0] for line in f]
    else:
        with open(p) as f:
            genes = [line.strip().split("\t")[0] for line in f]
    return genes


def load_barcodes(path):
    """Load cell barcodes from a tsv file (plain or gzipped)."""
    p = str(path)
    if p.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            barcodes = [line.strip().split("\t")[0] for line in f]
    else:
        with open(p) as f:
            barcodes = [line.strip().split("\t")[0] for line in f]
    return barcodes


def load_metadata(path, sep="\t"):
    """Load metadata, skipping the TYPE row if present."""
    df = pd.read_csv(path, sep=sep, index_col=0, compression="infer")
    if df.index[0] == "TYPE":
        df = df.iloc[1:]
    return df


def extract_marco_from_mtx(mtx_path, genes, barcodes, gene_name="MARCO"):
    """Extract MARCO expression vector from an MTX matrix (genes x cells)."""
    idx = None
    for i, g in enumerate(genes):
        if g.upper() == gene_name.upper():
            idx = i
            break
    if idx is None:
        return None

    mat = sio.mmread(mtx_path)
    # MTX is genes x cells — extract the MARCO row
    marco_row = mat.getrow(idx).toarray().flatten()
    return pd.Series(marco_row, index=barcodes, name=gene_name)


def summarize_expression(expr_series, metadata, group_col, dataset_name):
    """Summarize MARCO expression by a grouping variable."""
    df = metadata.copy()
    df["marco_expr"] = expr_series.reindex(df.index).fillna(0).values
    df["marco_positive"] = df["marco_expr"] > 0

    summary = df.groupby(group_col).agg(
        n_cells=("marco_expr", "size"),
        mean_expr=("marco_expr", "mean"),
        pct_positive=("marco_positive", lambda x: x.mean() * 100),
        max_expr=("marco_expr", "max"),
    ).round(3)

    summary["dataset"] = dataset_name
    return summary


# =============================================================================
# Phase 2 & 3: Process each dataset
# =============================================================================

all_summaries = []

# ---- SCP259: Human UC colon (Smillie et al. 2019) ----
print("=" * 60)
print("SCP259: Human UC Colon (Smillie et al. 2019)")
print("=" * 60)

scp259 = SCP / "SCP259"
meta259 = load_metadata(scp259 / "metadata" / "all.meta2.txt")

for compartment in ["Epi", "Imm"]:
    expr_dir = list((scp259 / "expression").iterdir())
    # Find the right expression dir by checking for the genes file
    for d in expr_dir:
        genes_file = d / f"{compartment}.genes.tsv"
        if genes_file.exists():
            genes = load_genes(genes_file)
            barcodes = load_barcodes(d / f"{compartment}.barcodes2.tsv")
            mtx_file = d / f"gene_sorted-{compartment}.matrix.mtx"
            marco = extract_marco_from_mtx(mtx_file, genes, barcodes, "MARCO")
            if marco is not None:
                # Join with metadata
                common = meta259.index.intersection(marco.index)
                marco_common = marco.loc[common]
                meta_common = meta259.loc[common]

                print(f"\n--- {compartment} compartment ---")
                print(f"  Cells with MARCO data: {len(common)}")
                print(f"  MARCO+ cells: {(marco_common > 0).sum()}")

                # By cell type (Cluster)
                by_cluster = summarize_expression(
                    marco_common, meta_common, "Cluster",
                    f"SCP259_{compartment}"
                )
                print(f"\n  By cell type:")
                pos = by_cluster[by_cluster["pct_positive"] > 0].sort_values("pct_positive", ascending=False)
                if len(pos) > 0:
                    print(pos.to_string(index=True))
                else:
                    print("  No MARCO+ cells found")

                # By health status
                by_health = summarize_expression(
                    marco_common, meta_common, "Health",
                    f"SCP259_{compartment}"
                )
                print(f"\n  By health status:")
                print(by_health.to_string(index=True))

                all_summaries.append(by_cluster.reset_index())
                all_summaries.append(by_health.reset_index())
            break

# ---- SCP1845: Human normal gut (immune cells) ----
print("\n" + "=" * 60)
print("SCP1845: Human Normal Gut Immune Cells")
print("=" * 60)

scp1845 = SCP / "SCP1845"
norm_dir = scp1845 / "expression" / "627509d0696e78b69ff49cb8"
genes_file = norm_dir / "global_features_normalized.tsv.gz"
barcodes_file = norm_dir / "global_barcodes_normalized.tsv.gz"
mtx_file = norm_dir / "global_normalized_matrix.mtx.gz"

genes = load_genes(genes_file)
barcodes = load_barcodes(barcodes_file)
meta1845 = load_metadata(scp1845 / "metadata" / "global_meta.tsv.gz")

marco = extract_marco_from_mtx(mtx_file, genes, barcodes, "MARCO")
if marco is not None:
    common = meta1845.index.intersection(marco.index)
    marco_common = marco.loc[common]
    meta_common = meta1845.loc[common]

    print(f"  Cells with MARCO data: {len(common)}")
    print(f"  MARCO+ cells: {(marco_common > 0).sum()}")

    # By cell type
    ct_col = "Manually_curated_celltype" if "Manually_curated_celltype" in meta_common.columns else "Cell_category"
    by_ct = summarize_expression(marco_common, meta_common, ct_col, "SCP1845")
    print(f"\n  By cell type ({ct_col}):")
    pos = by_ct[by_ct["pct_positive"] > 0].sort_values("pct_positive", ascending=False)
    if len(pos) > 0:
        print(pos.to_string(index=True))
    else:
        print("  No MARCO+ cells found")

    # By organ
    if "organ__ontology_label" in meta_common.columns:
        by_organ = summarize_expression(marco_common, meta_common, "organ__ontology_label", "SCP1845")
        print(f"\n  By organ:")
        print(by_organ.to_string(index=True))
        all_summaries.append(by_organ.reset_index())

    all_summaries.append(by_ct.reset_index())

# ---- SCP2038: Mouse colon scRNA-seq ----
print("\n" + "=" * 60)
print("SCP2038: Mouse Colon scRNA-seq (Parigi et al.)")
print("=" * 60)

scp2038 = SCP / "SCP2038"
# Use the h5ad file which is more convenient
h5ad_path = scp2038 / "anndata" / "scRNAseq.h5ad"
if h5ad_path.exists():
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Shape: {adata.shape}")

    if "Marco" in adata.var_names:
        marco_vals = np.asarray(adata[:, "Marco"].X.todense()).flatten()
        marco_series = pd.Series(marco_vals, index=adata.obs_names, name="Marco")

        print(f"  Marco+ cells: {(marco_vals > 0).sum()}")

        # Check available metadata columns
        print(f"  Metadata columns: {list(adata.obs.columns)}")

        # By cell type
        ct_cols = [c for c in adata.obs.columns if "cell" in c.lower() or "cluster" in c.lower() or "type" in c.lower() or "annot" in c.lower()]
        if ct_cols:
            ct_col = ct_cols[0]
        else:
            ct_col = adata.obs.columns[0]

        by_ct = summarize_expression(marco_series, adata.obs, ct_col, "SCP2038")
        print(f"\n  By cell type ({ct_col}):")
        pos = by_ct[by_ct["pct_positive"] > 0].sort_values("pct_positive", ascending=False)
        if len(pos) > 0:
            print(pos.to_string(index=True))
        else:
            print("  No Marco+ cells found")
        all_summaries.append(by_ct.reset_index())

        # By disease/condition if available
        disease_cols = [c for c in adata.obs.columns if "disease" in c.lower() or "condition" in c.lower() or "treatment" in c.lower()]
        for dc in disease_cols:
            by_disease = summarize_expression(marco_series, adata.obs, dc, "SCP2038")
            print(f"\n  By {dc}:")
            print(by_disease.to_string(index=True))
            all_summaries.append(by_disease.reset_index())
    else:
        print("  Marco not found in var_names")
else:
    print("  h5ad file not found, skipping")

# ---- SCP2760: Mouse colon scRNA-seq ----
print("\n" + "=" * 60)
print("SCP2760: Mouse Colon scRNA-seq")
print("=" * 60)

scp2760 = SCP / "SCP2760"
expr_dir = scp2760 / "expression" / "66cf96e51a832ffe7a16da09"
genes = load_genes(expr_dir / "features.tsv")
barcodes = load_barcodes(expr_dir / "barcodes.tsv")
meta2760 = load_metadata(scp2760 / "metadata" / "meta.tsv")

marco = extract_marco_from_mtx(expr_dir / "counts.mtx.gz", genes, barcodes, "Marco")
if marco is not None:
    common = meta2760.index.intersection(marco.index)
    marco_common = marco.loc[common]
    meta_common = meta2760.loc[common]

    print(f"  Cells with Marco data: {len(common)}")
    print(f"  Marco+ cells: {(marco_common > 0).sum()}")

    # By cell annotation
    ct_col = "cell_annotation" if "cell_annotation" in meta_common.columns else meta_common.columns[0]
    by_ct = summarize_expression(marco_common, meta_common, ct_col, "SCP2760")
    print(f"\n  By cell type ({ct_col}):")
    pos = by_ct[by_ct["pct_positive"] > 0].sort_values("pct_positive", ascending=False)
    if len(pos) > 0:
        print(pos.to_string(index=True))
    else:
        print("  No Marco+ cells found")
    all_summaries.append(by_ct.reset_index())

    # By disease/condition
    for dc in ["disease__ontology_label", "disease", "microbe"]:
        if dc in meta_common.columns:
            by_d = summarize_expression(marco_common, meta_common, dc, "SCP2760")
            print(f"\n  By {dc}:")
            print(by_d.to_string(index=True))
            all_summaries.append(by_d.reset_index())

# ---- SCP2771: Mouse colon Visium spatial (DSS colitis) ----
print("\n" + "=" * 60)
print("SCP2771: Mouse Colon Visium Spatial (DSS Colitis)")
print("=" * 60)

scp2771 = SCP / "SCP2771"
expr_dir1 = scp2771 / "expression" / "66e1e0765e4f1f173a036b83"
genes = load_genes(expr_dir1 / "genes.tsv")
barcodes = load_barcodes(expr_dir1 / "barcodes.tsv")
meta2771 = load_metadata(scp2771 / "metadata" / "dss_spatial_meta.tsv")

marco = extract_marco_from_mtx(expr_dir1 / "matrix.mtx", genes, barcodes, "Marco")
if marco is not None:
    common = meta2771.index.intersection(marco.index)
    marco_common = marco.loc[common]
    meta_common = meta2771.loc[common]

    print(f"  Cells with Marco data: {len(common)}")
    print(f"  Marco+ cells: {(marco_common > 0).sum()}")

    # By disease/timepoint
    if "dss_time_point" in meta_common.columns:
        by_tp = summarize_expression(marco_common, meta_common, "dss_time_point", "SCP2771")
        print(f"\n  By DSS timepoint:")
        print(by_tp.to_string(index=True))
        all_summaries.append(by_tp.reset_index())

    if "disease__ontology_label" in meta_common.columns:
        by_disease = summarize_expression(marco_common, meta_common, "disease__ontology_label", "SCP2771")
        print(f"\n  By disease status:")
        print(by_disease.to_string(index=True))
        all_summaries.append(by_disease.reset_index())

    # Spatial location
    cluster_file = scp2771 / "cluster" / "dss_spatial_cluster.tsv"
    if cluster_file.exists():
        spatial = load_metadata(cluster_file)
        spatial_common = spatial.loc[spatial.index.intersection(marco_common.index)]
        marco_spatial = marco_common.loc[spatial_common.index]
        spatial_common["marco_expr"] = marco_spatial.values
        spatial_common["marco_positive"] = marco_spatial.values > 0

        spatial_out = spatial_common[spatial_common["marco_positive"]].copy()
        if len(spatial_out) > 0:
            print(f"\n  Marco+ spots with spatial coordinates: {len(spatial_out)}")
            spatial_out.to_csv(RESULTS / "scp2771_marco_spatial.csv")
            print(f"  Saved spatial data to results/scp2771_marco_spatial.csv")

# ---- Gut Cell Atlas (already analyzed — include summary) ----
print("\n" + "=" * 60)
print("Gut Cell Atlas: Human Colon Immune (previously analyzed)")
print("=" * 60)

gut_atlas_h5ad = DATA / "gut-atlas" / "Colon_cell_atlas.h5ad"
if gut_atlas_h5ad.exists():
    adata_ga = sc.read_h5ad(gut_atlas_h5ad)
    if "MARCO" in adata_ga.var_names:
        marco_vals = np.asarray(adata_ga[:, "MARCO"].X.todense()).flatten()
        marco_series = pd.Series(marco_vals, index=adata_ga.obs_names, name="MARCO")

        print(f"  Shape: {adata_ga.shape}")
        print(f"  MARCO+ cells: {(marco_vals > 0).sum()}")

        ct_col = "cell_type" if "cell_type" in adata_ga.obs.columns else adata_ga.obs.columns[0]
        by_ct = summarize_expression(marco_series, adata_ga.obs, ct_col, "GutAtlas")
        pos = by_ct[by_ct["pct_positive"] > 0].sort_values("pct_positive", ascending=False)
        print(f"\n  By cell type:")
        if len(pos) > 0:
            print(pos.to_string(index=True))
        all_summaries.append(by_ct.reset_index())

        if "region" in adata_ga.obs.columns:
            by_region = summarize_expression(marco_series, adata_ga.obs, "region", "GutAtlas")
            print(f"\n  By region:")
            print(by_region.to_string(index=True))
            all_summaries.append(by_region.reset_index())

# =============================================================================
# Save combined results
# =============================================================================
if all_summaries:
    combined = pd.concat(all_summaries, ignore_index=True)
    combined.to_csv(RESULTS / "marco_all_datasets_summary.csv", index=False)
    print(f"\n{'=' * 60}")
    print(f"Combined summary saved to: {RESULTS / 'marco_all_datasets_summary.csv'}")
    print(f"Total summary rows: {len(combined)}")
