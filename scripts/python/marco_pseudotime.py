"""
Pseudotime analysis of myeloid cells to trace monocyte-to-macrophage
differentiation and locate MARCO+ cells along the trajectory.

Uses scanpy diffusion pseudotime (DPT) on:
1. Gut Atlas — colon-specific myeloid cells
2. SCP1845 — multi-organ myeloid cells (subset to gut)
3. SCP259 — UC colon myeloid cells
"""

import numpy as np
import pandas as pd
import scanpy as sc
import warnings
from pathlib import Path
from scipy.io import mmread
import gzip

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

BASE = Path("/Volumes/repos/uc")
DATA = BASE / "data"
RESULTS = BASE / "analysis" / "results"


def run_pseudotime(adata, marco_gene, dataset_name, root_type=None, ct_col="cell_type"):
    """
    Run diffusion pseudotime on myeloid cells.
    Root cell is chosen from the specified root_type (e.g., monocytes).
    """
    print(f"\n{'='*70}")
    print(f"Pseudotime: {dataset_name}")
    print(f"{'='*70}")
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    # Get MARCO expression
    marco_vals = np.asarray(adata[:, marco_gene].X.todense()).flatten()
    adata.obs["marco_expr"] = marco_vals
    adata.obs["marco_status"] = ["MARCO+" if v > 0 else "MARCO-" for v in marco_vals]

    # Standard preprocessing
    print("  Preprocessing...")
    adata_pp = adata.copy()

    # Check if data looks raw (integer counts) or normalized
    sample_vals = adata_pp.X[:100, :100].toarray() if hasattr(adata_pp.X, 'toarray') else adata_pp.X[:100, :100]
    is_raw = np.allclose(sample_vals, sample_vals.astype(int))

    if is_raw:
        print("  Data appears raw — normalizing and log-transforming")
        sc.pp.normalize_total(adata_pp, target_sum=1e4)
        sc.pp.log1p(adata_pp)
    else:
        print("  Data appears pre-normalized")

    # HVG, PCA, neighbors, UMAP
    sc.pp.highly_variable_genes(adata_pp, n_top_genes=2000, subset=True)
    sc.tl.pca(adata_pp, n_comps=30)
    sc.pp.neighbors(adata_pp, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_pp)
    sc.tl.diffmap(adata_pp, n_comps=10)

    # Find root cell — pick a monocyte (or specified root type) with lowest MARCO
    if root_type:
        root_mask = adata_pp.obs[ct_col].isin(root_type if isinstance(root_type, list) else [root_type])
        if root_mask.sum() > 0:
            root_candidates = adata_pp.obs[root_mask]
            # Pick the monocyte with lowest MARCO expression
            root_idx = root_candidates["marco_expr"].idxmin()
            root_pos = list(adata_pp.obs_names).index(root_idx)
            print(f"  Root cell: {root_idx} (type: {adata_pp.obs.loc[root_idx, ct_col]})")
        else:
            print(f"  Warning: root type {root_type} not found, using first cell")
            root_pos = 0
    else:
        root_pos = 0

    adata_pp.uns["iroot"] = root_pos
    sc.tl.dpt(adata_pp)

    # Transfer results back
    adata.obs["dpt_pseudotime"] = adata_pp.obs["dpt_pseudotime"].values
    adata.obs["umap_1"] = adata_pp.obsm["X_umap"][:, 0]
    adata.obs["umap_2"] = adata_pp.obsm["X_umap"][:, 1]

    # Analyze MARCO along pseudotime
    print(f"\n  --- MARCO along pseudotime ---")

    # Bin pseudotime into deciles
    valid_pt = adata.obs[adata.obs["dpt_pseudotime"] < np.inf].copy()
    valid_pt["pt_bin"] = pd.qcut(valid_pt["dpt_pseudotime"], q=10, labels=False, duplicates="drop")

    pt_summary = valid_pt.groupby("pt_bin").agg(
        n_cells=("marco_expr", "size"),
        mean_marco=("marco_expr", "mean"),
        pct_marco_pos=("marco_status", lambda x: (x == "MARCO+").mean() * 100),
        mean_pseudotime=("dpt_pseudotime", "mean"),
    ).round(3)
    print(pt_summary.to_string())

    # MARCO+ vs MARCO- pseudotime distribution
    marco_pos_pt = valid_pt[valid_pt["marco_status"] == "MARCO+"]["dpt_pseudotime"]
    marco_neg_pt = valid_pt[valid_pt["marco_status"] == "MARCO-"]["dpt_pseudotime"]

    print(f"\n  MARCO+ pseudotime: mean={marco_pos_pt.mean():.3f}, median={marco_pos_pt.median():.3f}")
    print(f"  MARCO- pseudotime: mean={marco_neg_pt.mean():.3f}, median={marco_neg_pt.median():.3f}")

    from scipy.stats import mannwhitneyu
    if len(marco_pos_pt) > 0 and len(marco_neg_pt) > 0:
        stat, pval = mannwhitneyu(marco_pos_pt, marco_neg_pt, alternative="two-sided")
        print(f"  Mann-Whitney U test: U={stat:.0f}, p={pval:.2e}")

    # Cell type along pseudotime
    print(f"\n  --- Cell types along pseudotime ---")
    ct_pt = valid_pt.groupby(ct_col).agg(
        n_cells=("dpt_pseudotime", "size"),
        mean_pt=("dpt_pseudotime", "mean"),
        median_pt=("dpt_pseudotime", "median"),
        pct_marco_pos=("marco_status", lambda x: (x == "MARCO+").mean() * 100),
    ).round(3).sort_values("mean_pt")
    print(ct_pt.to_string())

    # Save results
    out_df = valid_pt[["dpt_pseudotime", "umap_1", "umap_2", "marco_expr", "marco_status", ct_col]].copy()
    out_path = RESULTS / f"{dataset_name}_pseudotime.csv"
    out_df.to_csv(out_path)
    print(f"\n  Results saved to: {out_path.name}")

    return adata, adata_pp


# =============================================================================
# 1. Gut Cell Atlas — colon myeloid
# =============================================================================
print("Loading Gut Cell Atlas...")
adata_ga = sc.read_h5ad(DATA / "gut-atlas" / "Colon_cell_atlas.h5ad")

myeloid_types = ["Monocyte", "LYVE1 Macrophage", "Macrophage", "cycling DCs",
                 "cDC1", "cDC2", "pDC"]
myeloid_ga = adata_ga[adata_ga.obs["cell_type"].isin(myeloid_types)].copy()

run_pseudotime(
    myeloid_ga, "MARCO", "gut_atlas",
    root_type="Monocyte", ct_col="cell_type"
)

# =============================================================================
# 2. SCP1845 — gut-only myeloid cells
# =============================================================================
print("\n\nLoading SCP1845...")
scp1845 = DATA / "single-cell-portal" / "SCP1845"
norm_dir = scp1845 / "expression" / "627509d0696e78b69ff49cb8"

with gzip.open(norm_dir / "global_features_normalized.tsv.gz", "rt") as f:
    genes = [line.strip().split("\t")[0] for line in f]
with gzip.open(norm_dir / "global_barcodes_normalized.tsv.gz", "rt") as f:
    barcodes = [line.strip().split("\t")[0] for line in f]

meta = pd.read_csv(scp1845 / "metadata" / "global_meta.tsv.gz", sep="\t",
                    index_col=0, compression="gzip")
if meta.index[0] == "TYPE":
    meta = meta.iloc[1:]

print("  Reading expression matrix...")
mat = mmread(norm_dir / "global_normalized_matrix.mtx.gz").tocsc()

adata_1845 = sc.AnnData(X=mat.T, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
common = adata_1845.obs.index.intersection(meta.index)
adata_1845 = adata_1845[common].copy()
for col in meta.columns:
    adata_1845.obs[col] = meta.loc[common, col].values

# Subset to gut myeloid cells only
gut_organs = ["jejunal epithelium", "lamina propria", "ileum", "caecum",
              "sigmoid colon", "transverse colon", "duodenum"]
myeloid_types_1845 = [
    "Alveolar macrophages", "Intermediate macrophages", "Classical monocytes",
    "Nonclassical monocytes", "Erythrophagocytic macrophages",
    "Intestinal macrophages", "Cycling", "DC1", "DC2", "migDC", "pDC",
]

gut_mask = adata_1845.obs["organ__ontology_label"].isin(gut_organs)
myeloid_mask = adata_1845.obs["Manually_curated_celltype"].isin(myeloid_types_1845)
gut_myeloid = adata_1845[gut_mask & myeloid_mask].copy()

print(f"  Gut myeloid cells: {gut_myeloid.n_obs}")

if gut_myeloid.n_obs >= 50:
    run_pseudotime(
        gut_myeloid, "MARCO", "scp1845_gut",
        root_type="Classical monocytes", ct_col="Manually_curated_celltype"
    )
else:
    print("  Too few gut myeloid cells, running all-organ myeloid instead...")
    all_myeloid = adata_1845[myeloid_mask].copy()
    run_pseudotime(
        all_myeloid, "MARCO", "scp1845_all",
        root_type="Classical monocytes", ct_col="Manually_curated_celltype"
    )

# =============================================================================
# 3. SCP259 — UC colon myeloid
# =============================================================================
print("\n\nLoading SCP259...")
scp259 = DATA / "single-cell-portal" / "SCP259"

for d in (scp259 / "expression").iterdir():
    if (d / "Imm.genes.tsv").exists():
        with open(d / "Imm.genes.tsv") as f:
            genes_259 = [line.strip() for line in f]
        with open(d / "Imm.barcodes2.tsv") as f:
            barcodes_259 = [line.strip() for line in f]
        print("  Reading expression matrix...")
        mat_259 = mmread(d / "gene_sorted-Imm.matrix.mtx").tocsc()
        break

meta_259 = pd.read_csv(scp259 / "metadata" / "all.meta2.txt", sep="\t",
                        index_col=0, low_memory=False)
if meta_259.index[0] == "TYPE":
    meta_259 = meta_259.iloc[1:]

adata_259 = sc.AnnData(X=mat_259.T, obs=pd.DataFrame(index=barcodes_259),
                        var=pd.DataFrame(index=genes_259))
common = adata_259.obs.index.intersection(meta_259.index)
adata_259 = adata_259[common].copy()
for col in meta_259.columns:
    adata_259.obs[col] = meta_259.loc[common, col].values

# Subset to myeloid
myeloid_259 = ["Macrophages", "Inflammatory Monocytes", "Cycling Monocytes", "DC1", "DC2"]
myeloid_mask = adata_259.obs["Cluster"].isin(myeloid_259)
myeloid_adata = adata_259[myeloid_mask].copy()

run_pseudotime(
    myeloid_adata, "MARCO", "scp259",
    root_type="Inflammatory Monocytes", ct_col="Cluster"
)

# Also check: does pseudotime differ by Health status among MARCO+ cells?
print("\n  --- MARCO+ pseudotime by Health status (SCP259) ---")
valid = myeloid_adata.obs[
    (myeloid_adata.obs["dpt_pseudotime"] < np.inf) &
    (myeloid_adata.obs["marco_status"] == "MARCO+")
].copy()
if len(valid) > 0:
    by_health = valid.groupby("Health").agg(
        n_marco_pos=("dpt_pseudotime", "size"),
        mean_pt=("dpt_pseudotime", "mean"),
        median_pt=("dpt_pseudotime", "median"),
    ).round(3)
    print(by_health.to_string())

print("\n" + "=" * 70)
print("All pseudotime analyses complete.")
