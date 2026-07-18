"""Step 4: Integrate — PCA, Harmony batch correction, neighbors/UMAP, and a
sample-driven-clustering sanity check.

Order (per scanpy/harmonypy docs: Harmony adjusts principal components, so
it must run after PCA but before the neighbor graph): PCA on the scaled HVG
matrix -> Harmony (batch key = "sample") -> neighbors on the Harmony-adjusted
PCs -> UMAP -> Leiden clustering, purely as a per-cluster sample-composition
diagnostic (not a cell-type call — that's Step 5).

n_pcs=30: standard practice for feeding Harmony (Korsunsky et al. 2019;
consistent with recent integration-benchmarking guidance, e.g. Nature
Methods 2025 feature-selection/integration benchmark), not tuned against
this dataset specifically.

Integration-QC is a visual UMAP-by-sample check plus a per-cluster sample-
composition table (max-sample-fraction per cluster) — a kBET/LISI-style
stat without pulling in scib-metrics, per docs/agents/codebase-agent.md's
"don't over-engineer integration QC for this scope" guidance.

MK447 watch: Step 3's logic-agent review flagged that MK447 (low-depth,
QC-flagged — see qc.py LOW_CONFIDENCE_SAMPLES) carries stress-response/
dissociation-artifact signal that could produce an MK447-isolated cluster
post-integration for technical reasons, not novel biology. This module
does not resolve that call — it only reports whether any cluster is
MK447-dominated, for a human/logic-agent to interpret.
"""

from pathlib import Path

import harmonypy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parent.parent
IN_PATH = REPO_ROOT / "data" / "processed" / "normalized.h5ad"
OUT_PATH = REPO_ROOT / "data" / "processed" / "integrated.h5ad"
PLOT_DIR = REPO_ROOT / "data" / "processed" / "integrate_plots"

N_PCS = 30
LEIDEN_RESOLUTION = 1.0
DOMINANCE_FLAG_THRESHOLD = 0.8  # a cluster >=80% one sample is worth flagging


def run_pca(adata: sc.AnnData) -> sc.AnnData:
    sc.pp.pca(adata, n_comps=N_PCS)
    return adata


def run_harmony(adata: sc.AnnData) -> sc.AnnData:
    """Call harmonypy directly rather than sc.external.pp.harmony_integrate.

    The installed harmonypy (2.0.0) returns Z_corr already as cells x PCs,
    but scanpy's wrapper unconditionally transposes it (`Z_corr.T`),
    assuming the pre-2.0 PCs x cells orientation — a real version
    incompatibility, not a data issue (reproduces on a synthetic array with
    the installed harmonypy). Orientation is checked against n_obs here so
    this keeps working regardless of which convention a given harmonypy
    version uses.
    """
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["sample"])
    z_corr = np.asarray(ho.Z_corr)
    if z_corr.shape[0] != adata.n_obs:
        z_corr = z_corr.T
    assert z_corr.shape == (adata.n_obs, N_PCS), f"unexpected Z_corr shape {z_corr.shape}"
    adata.obsm["X_pca_harmony"] = z_corr
    return adata


def run_neighbors_umap(adata: sc.AnnData) -> sc.AnnData:
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    sc.tl.umap(adata)
    return adata


def run_leiden(adata: sc.AnnData, resolution: float = LEIDEN_RESOLUTION) -> sc.AnnData:
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2)
    return adata


def plot_umap_by_sample_and_cluster(adata: sc.AnnData) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(adata, color="sample", ax=axes[0], show=False, title="UMAP by sample")
    sc.pl.umap(adata, color="leiden", ax=axes[1], show=False, title="UMAP by leiden cluster")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "umap_sample_vs_cluster.png", dpi=150)
    plt.close(fig)
    print(f"Saved {PLOT_DIR / 'umap_sample_vs_cluster.png'}")


def sample_driven_clustering_check(adata: sc.AnnData) -> pd.DataFrame:
    """Per-cluster sample composition: n cells, dominant sample, its
    fraction. Clusters at/above DOMINANCE_FLAG_THRESHOLD are printed as
    flagged — most relevant here is whether any flagged cluster is MK447.
    """
    comp = (
        pd.crosstab(adata.obs["leiden"], adata.obs["sample"])
        .div(adata.obs.groupby("leiden", observed=True).size(), axis=0)
    )
    dominant_sample = comp.idxmax(axis=1)
    dominant_frac = comp.max(axis=1)
    sizes = adata.obs.groupby("leiden", observed=True).size()

    summary = pd.DataFrame({
        "n_cells": sizes,
        "dominant_sample": dominant_sample,
        "dominant_frac": dominant_frac,
    }).sort_values("n_cells", ascending=False)

    print(f"\nPer-cluster sample composition ({len(summary)} leiden clusters, "
          f"resolution={LEIDEN_RESOLUTION}):")
    print(f"{'cluster':>7s} {'n_cells':>8s} {'dominant_sample':>16s} {'dominant_frac':>14s}  flag")
    for cluster, row in summary.iterrows():
        flagged = row["dominant_frac"] >= DOMINANCE_FLAG_THRESHOLD
        note = ""
        if flagged:
            note = f"FLAGGED — {'MK447' if row['dominant_sample'] == 'MK447' else 'sample'}-dominated"
        print(f"{cluster:>7s} {int(row['n_cells']):8d} {row['dominant_sample']:>16s} "
              f"{row['dominant_frac']:14.2f}  {note}")

    n_flagged = int((summary["dominant_frac"] >= DOMINANCE_FLAG_THRESHOLD).sum())
    mk447_flagged = summary[(summary["dominant_frac"] >= DOMINANCE_FLAG_THRESHOLD)
                             & (summary["dominant_sample"] == "MK447")]
    print(f"\n{n_flagged} / {len(summary)} clusters are >={DOMINANCE_FLAG_THRESHOLD:.0%} "
          "one sample.")
    if not mk447_flagged.empty:
        print("MK447-dominated cluster(s) present — flagging for logic-agent review "
              "per the Step 3 carry-forward note (technical dropout/stress-response "
              "artifact vs. novel biology is not resolved here).")
    else:
        print("No MK447-dominated cluster — MK447 cells integrated with other samples' "
              "clusters rather than forming an isolated group.")
    return summary


def integrate() -> sc.AnnData:
    adata = sc.read_h5ad(IN_PATH)
    n_cells_in, n_genes_in = adata.n_obs, adata.n_vars
    adata = run_pca(adata)
    adata = run_harmony(adata)
    adata = run_neighbors_umap(adata)
    adata = run_leiden(adata)
    plot_umap_by_sample_and_cluster(adata)
    sample_driven_clustering_check(adata)
    print(f"\nInput:  {n_cells_in} cells x {n_genes_in} genes ({IN_PATH.name})")
    print(f"Output: {adata.n_obs} cells, X_pca_harmony dim={adata.obsm['X_pca_harmony'].shape[1]}")
    return adata


if __name__ == "__main__":
    integrated = integrate()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    integrated.write_h5ad(OUT_PATH)
    print(f"\nSaved integrated data (Step 5 input) to {OUT_PATH}")
