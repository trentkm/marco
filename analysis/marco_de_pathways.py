"""
Differential expression and pathway analysis: MARCO+ vs MARCO- myeloid cells.
Identifies co-upregulated/downregulated genes and enriched GO/KEGG pathways.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

BASE = Path("/Volumes/repos/uc")
DATA = BASE / "data"
RESULTS = BASE / "analysis" / "results"

# Gene set libraries to test
GENE_SET_LIBS = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2021_Human",
]

MOUSE_GENE_SET_LIBS = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2019_Mouse",
]


def run_de_and_enrichment(adata, marco_gene, group_col, myeloid_labels,
                          dataset_name, gene_set_libs, min_marco_cells=20):
    """
    Run DE (MARCO+ vs MARCO-) within myeloid cells, then pathway enrichment.
    """
    print(f"\n{'='*60}")
    print(f"{dataset_name}: DE & Pathway Analysis")
    print(f"{'='*60}")

    # Subset to myeloid cells
    mask = adata.obs[group_col].isin(myeloid_labels)
    myeloid = adata[mask].copy()
    print(f"  Myeloid cells: {myeloid.n_obs}")

    # Get MARCO expression
    marco_vals = np.asarray(myeloid[:, marco_gene].X.todense()).flatten()
    myeloid.obs["marco_status"] = ["MARCO+" if v > 0 else "MARCO-" for v in marco_vals]

    n_pos = (myeloid.obs["marco_status"] == "MARCO+").sum()
    n_neg = (myeloid.obs["marco_status"] == "MARCO-").sum()
    print(f"  MARCO+ cells: {n_pos}")
    print(f"  MARCO- cells: {n_neg}")

    if n_pos < min_marco_cells:
        print(f"  Skipping — too few MARCO+ cells (<{min_marco_cells})")
        return None, None

    # Run DE: Wilcoxon rank-sum
    sc.tl.rank_genes_groups(myeloid, groupby="marco_status", groups=["MARCO+"],
                            reference="MARCO-", method="wilcoxon", use_raw=False)

    # Extract results
    de_results = sc.get.rank_genes_groups_df(myeloid, group="MARCO+")
    de_results = de_results.sort_values("scores", ascending=False)

    # Filter significant genes
    sig = de_results[de_results["pvals_adj"] < 0.05].copy()
    up = sig[sig["logfoldchanges"] > 0.25].head(50)
    down = sig[sig["logfoldchanges"] < -0.25].head(50)

    print(f"\n  Significant DE genes (padj < 0.05): {len(sig)}")
    print(f"  Upregulated in MARCO+ (logFC > 0.25): {len(sig[sig['logfoldchanges'] > 0.25])}")
    print(f"  Downregulated in MARCO+ (logFC < -0.25): {len(sig[sig['logfoldchanges'] < -0.25])}")

    if len(up) > 0:
        print(f"\n  Top 25 upregulated in MARCO+ cells:")
        print(up[["names", "logfoldchanges", "pvals_adj", "scores"]].head(25).to_string(index=False))

    if len(down) > 0:
        print(f"\n  Top 25 downregulated in MARCO+ cells:")
        down_sorted = sig[sig["logfoldchanges"] < -0.25].sort_values("logfoldchanges").head(25)
        print(down_sorted[["names", "logfoldchanges", "pvals_adj", "scores"]].to_string(index=False))

    # Save full DE results
    de_outfile = RESULTS / f"{dataset_name}_marco_de_results.csv"
    de_results.to_csv(de_outfile, index=False)
    print(f"\n  Full DE results saved to: {de_outfile.name}")

    # Pathway enrichment on upregulated genes
    up_genes = list(sig[sig["logfoldchanges"] > 0.25]["names"])
    down_genes = list(sig[sig["logfoldchanges"] < -0.25]["names"])

    enrichment_results = {}
    for direction, gene_list in [("upregulated", up_genes), ("downregulated", down_genes)]:
        if len(gene_list) < 5:
            print(f"\n  Skipping {direction} enrichment — too few genes ({len(gene_list)})")
            continue

        print(f"\n  Running enrichment on {len(gene_list)} {direction} genes...")
        for lib in gene_set_libs:
            try:
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=lib,
                    organism="human" if marco_gene == "MARCO" else "mouse",
                    outdir=None,
                    no_plot=True,
                )
                sig_terms = enr.results[enr.results["Adjusted P-value"] < 0.05]
                if len(sig_terms) > 0:
                    print(f"\n  {lib} — {direction} ({len(sig_terms)} significant terms):")
                    top = sig_terms.head(15)[["Term", "Adjusted P-value", "Combined Score", "Genes"]].copy()
                    top["Adjusted P-value"] = top["Adjusted P-value"].map("{:.2e}".format)
                    top["Combined Score"] = top["Combined Score"].round(1)
                    print(top.to_string(index=False))
                    enrichment_results[f"{lib}_{direction}"] = sig_terms
            except Exception as e:
                print(f"  {lib} enrichment failed: {e}")

    # Save enrichment results
    if enrichment_results:
        all_enr = pd.concat(
            [df.assign(source=name) for name, df in enrichment_results.items()],
            ignore_index=True
        )
        enr_outfile = RESULTS / f"{dataset_name}_marco_enrichment.csv"
        all_enr.to_csv(enr_outfile, index=False)
        print(f"\n  Enrichment results saved to: {enr_outfile.name}")

    return de_results, enrichment_results


# =============================================================================
# Dataset 1: Gut Cell Atlas (Human colon immune)
# =============================================================================
print("Loading Gut Cell Atlas...")
adata_ga = sc.read_h5ad(DATA / "gut-atlas" / "Colon_cell_atlas.h5ad")

myeloid_types_ga = ["Monocyte", "LYVE1 Macrophage", "Macrophage", "cycling DCs",
                    "cDC1", "cDC2", "pDC"]
run_de_and_enrichment(
    adata_ga, "MARCO", "cell_type", myeloid_types_ga,
    "gut_atlas", GENE_SET_LIBS
)

# =============================================================================
# Dataset 2: SCP1845 (Human cross-tissue immune)
# =============================================================================
print("\nLoading SCP1845...")
# This is a large MTX dataset — use the normalized matrix
from scipy.io import mmread
import gzip

scp1845 = DATA / "single-cell-portal" / "SCP1845"
norm_dir = scp1845 / "expression" / "627509d0696e78b69ff49cb8"

with gzip.open(norm_dir / "global_features_normalized.tsv.gz", "rt") as f:
    genes_1845 = [line.strip().split("\t")[0] for line in f]
with gzip.open(norm_dir / "global_barcodes_normalized.tsv.gz", "rt") as f:
    barcodes_1845 = [line.strip().split("\t")[0] for line in f]

meta_1845 = pd.read_csv(scp1845 / "metadata" / "global_meta.tsv.gz", sep="\t",
                         index_col=0, compression="gzip")
if meta_1845.index[0] == "TYPE":
    meta_1845 = meta_1845.iloc[1:]

print("  Reading expression matrix...")
mat_1845 = mmread(norm_dir / "global_normalized_matrix.mtx.gz").tocsc()

# Build AnnData
adata_1845 = sc.AnnData(
    X=mat_1845.T,  # cells x genes
    obs=pd.DataFrame(index=barcodes_1845),
    var=pd.DataFrame(index=genes_1845),
)
# Add metadata
common = adata_1845.obs.index.intersection(meta_1845.index)
adata_1845 = adata_1845[common].copy()
for col in meta_1845.columns:
    adata_1845.obs[col] = meta_1845.loc[common, col].values

myeloid_types_1845 = [
    "Alveolar macrophages", "Intermediate macrophages", "Classical monocytes",
    "Nonclassical monocytes", "Erythrophagocytic macrophages",
    "Intestinal macrophages", "Cycling", "DC1", "DC2", "migDC", "pDC",
]
run_de_and_enrichment(
    adata_1845, "MARCO", "Manually_curated_celltype", myeloid_types_1845,
    "scp1845", GENE_SET_LIBS
)

# =============================================================================
# Dataset 3: SCP259 (Human UC — immune compartment)
# =============================================================================
print("\nLoading SCP259 immune compartment...")
scp259 = DATA / "single-cell-portal" / "SCP259"

# Find expression dir
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

adata_259 = sc.AnnData(
    X=mat_259.T,
    obs=pd.DataFrame(index=barcodes_259),
    var=pd.DataFrame(index=genes_259),
)
common = adata_259.obs.index.intersection(meta_259.index)
adata_259 = adata_259[common].copy()
for col in meta_259.columns:
    adata_259.obs[col] = meta_259.loc[common, col].values

myeloid_types_259 = ["Macrophages", "Inflammatory Monocytes", "Cycling Monocytes",
                     "DC1", "DC2"]
run_de_and_enrichment(
    adata_259, "MARCO", "Cluster", myeloid_types_259,
    "scp259", GENE_SET_LIBS
)

# =============================================================================
# Dataset 4: SCP2771 (Mouse DSS colitis — Visium spatial, D12 only)
# =============================================================================
print("\nLoading SCP2771 (D12 DSS spots only)...")
scp2771 = DATA / "single-cell-portal" / "SCP2771"
expr_dir = scp2771 / "expression" / "66e1e0765e4f1f173a036b83"

with open(expr_dir / "genes.tsv") as f:
    genes_2771 = [line.strip().split("\t")[0] for line in f]
with open(expr_dir / "barcodes.tsv") as f:
    barcodes_2771 = [line.strip().split("\t")[0] for line in f]

mat_2771 = mmread(expr_dir / "matrix.mtx").tocsc()

meta_2771 = pd.read_csv(scp2771 / "metadata" / "dss_spatial_meta.tsv", sep="\t",
                         index_col=0)
if meta_2771.index[0] == "TYPE":
    meta_2771 = meta_2771.iloc[1:]

adata_2771 = sc.AnnData(
    X=mat_2771.T,
    obs=pd.DataFrame(index=barcodes_2771),
    var=pd.DataFrame(index=genes_2771),
)
common = adata_2771.obs.index.intersection(meta_2771.index)
adata_2771 = adata_2771[common].copy()
for col in meta_2771.columns:
    adata_2771.obs[col] = meta_2771.loc[common, col].values

# For spatial data, compare D12 spots (peak inflammation) Marco+ vs Marco-
d12 = adata_2771[adata_2771.obs["dss_time_point"] == "D12"].copy()
print(f"  D12 spots: {d12.n_obs}")

marco_vals = np.asarray(d12[:, "Marco"].X.todense()).flatten()
d12.obs["marco_status"] = ["MARCO+" if v > 0 else "MARCO-" for v in marco_vals]
n_pos = (d12.obs["marco_status"] == "MARCO+").sum()
print(f"  Marco+ spots at D12: {n_pos}")

if n_pos >= 20:
    sc.tl.rank_genes_groups(d12, groupby="marco_status", groups=["MARCO+"],
                            reference="MARCO-", method="wilcoxon", use_raw=False)
    de_2771 = sc.get.rank_genes_groups_df(d12, group="MARCO+")
    sig = de_2771[de_2771["pvals_adj"] < 0.05]
    up = sig[sig["logfoldchanges"] > 0.25]
    down = sig[sig["logfoldchanges"] < -0.25]

    print(f"\n  Significant DE genes: {len(sig)}")
    print(f"  Upregulated in Marco+ spots: {len(up)}")
    print(f"  Downregulated in Marco+ spots: {len(down)}")

    if len(up) > 0:
        print(f"\n  Top 25 upregulated in Marco+ D12 spots:")
        print(up.head(25)[["names", "logfoldchanges", "pvals_adj", "scores"]].to_string(index=False))

    if len(down) > 0:
        print(f"\n  Top 25 downregulated in Marco+ D12 spots:")
        print(down.sort_values("logfoldchanges").head(25)[["names", "logfoldchanges", "pvals_adj", "scores"]].to_string(index=False))

    de_2771.to_csv(RESULTS / "scp2771_d12_marco_de_results.csv", index=False)

    # Enrichment — convert mouse genes to uppercase for human gene set libs
    up_genes_mouse = list(up["names"])
    down_genes_mouse = list(down.sort_values("logfoldchanges")["names"])

    for direction, gene_list in [("upregulated", up_genes_mouse), ("downregulated", down_genes_mouse)]:
        if len(gene_list) < 5:
            continue
        # Mouse gene names: capitalize first letter only -> uppercase for enrichr
        gene_list_upper = [g.upper() for g in gene_list]
        print(f"\n  Running enrichment on {len(gene_list)} {direction} genes (D12 Marco+ spots)...")
        for lib in GENE_SET_LIBS:  # Use human libs with uppercased mouse genes
            try:
                enr = gp.enrichr(
                    gene_list=gene_list_upper,
                    gene_sets=lib,
                    organism="human",
                    outdir=None,
                    no_plot=True,
                )
                sig_terms = enr.results[enr.results["Adjusted P-value"] < 0.05]
                if len(sig_terms) > 0:
                    print(f"\n  {lib} — {direction} ({len(sig_terms)} significant terms):")
                    top = sig_terms.head(15)[["Term", "Adjusted P-value", "Combined Score", "Genes"]].copy()
                    top["Adjusted P-value"] = top["Adjusted P-value"].map("{:.2e}".format)
                    top["Combined Score"] = top["Combined Score"].round(1)
                    print(top.to_string(index=False))
            except Exception as e:
                print(f"  {lib} failed: {e}")

print("\n" + "=" * 60)
print("All analyses complete.")
