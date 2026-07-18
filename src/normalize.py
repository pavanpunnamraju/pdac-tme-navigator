"""Step 3: Normalize — shifted-log normalization, HVG selection, scaling.

Shifted-logarithm (sc.pp.normalize_total + sc.pp.log1p) is used as the
default normalization, per docs/agents/codebase-agent.md best practices: it
is the most robust general-purpose choice for downstream dimensionality
reduction (PCA/Harmony in Step 4), unlike more complex model-based methods
(e.g. scran pooling, analytic Pearson residuals for normalization itself)
which need more tuning to pay off.

HVG selection footnote: analytic Pearson residuals (sc.experimental.pp.
highly_variable_genes, flavor="pearson_residuals") are reported in the
literature (Lause et al. 2021) to outperform log-normalization-based HVG
selection specifically for detecting rare cell populations, which matters
here given the project cares about tumor-infiltrating rare TME populations.
Not implemented — the raw counts are still present in adata.raw / a layer
if this needs revisiting, and per the same best-practices doc this should
only be added if the simple seurat-flavor dispersion approach demonstrably
falls short (e.g. known rare population markers not surviving HVG filtering
downstream), not preemptively.

HVG selection uses batch_key="sample": with Step 4 (Harmony) integrating
across all 8 patients next, an HVG set not dominated by inter-patient
variation is the safer default going into integration.

MK447 is excluded from the HVG *voting* only (its cells stay in the final
output). An empirical sensitivity check — HVG mask with all 8 samples
voting vs. with MK447 dropped from voting, same n_top_genes/batch_key —
found only 88.95% overlap (1779/2000 genes), below the 90% bar for calling
this a non-issue. MK447 is the smallest, depth-limited, QC-flagged sample
(see qc.py LOW_CONFIDENCE_SAMPLES); its inclusion in voting measurably
pulls in genes the other 7 samples don't select, consistent with it
contributing noise-floor variance rather than signal to the dispersion
ranking. See select_hvgs() for the mechanics.

Order: normalize_total -> log1p -> HVG selection -> scale -> subset to HVGs.
Raw counts are preserved in adata.layers["counts"] before the in-place
normalize_total/log1p step (NOT in adata.raw — .raw holds the full
log-normalized matrix, pre-HVG-subset, for reference/plotting only). PCA
itself is out of scope for this step (Step 4).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parent.parent
IN_PATH = REPO_ROOT / "data" / "processed" / "qc_filtered.h5ad"
OUT_PATH = REPO_ROOT / "data" / "processed" / "normalized.h5ad"
PLOT_DIR = REPO_ROOT / "data" / "processed" / "normalize_plots"

N_TOP_GENES = 2000
SCALE_MAX_VALUE = 10
HVG_VOTING_EXCLUDE = ["MK447"]  # noise floor pulls in spurious genes — see module docstring


def normalize_and_log(adata: sc.AnnData) -> sc.AnnData:
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata


def select_hvgs(adata: sc.AnnData, n_top_genes: int = N_TOP_GENES) -> sc.AnnData:
    """HVG mask computed on samples excluding HVG_VOTING_EXCLUDE, then
    applied back to all samples (whose cells are all kept in the output).
    See module docstring for the sensitivity check behind this exclusion.
    """
    voting = adata[~adata.obs["sample"].isin(HVG_VOTING_EXCLUDE)].copy()
    sc.pp.highly_variable_genes(
        voting, n_top_genes=n_top_genes, flavor="seurat", batch_key="sample"
    )
    for col in ["means", "dispersions", "dispersions_norm", "highly_variable"]:
        adata.var[col] = voting.var[col].values
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"HVG selection: {n_hvg} / {adata.n_vars} genes flagged highly variable "
          f"(target {n_top_genes}, batch_key='sample', voting excludes {HVG_VOTING_EXCLUDE})")
    return adata


def plot_hvg_mean_variance(adata: sc.AnnData) -> None:
    """Mean vs. normalized dispersion, HVGs highlighted — sanity check that
    the HVG flag is actually tracking dispersion rather than picking noise."""
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    means = adata.var["means"]
    disp = adata.var["dispersions_norm"]
    hvg = adata.var["highly_variable"]

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(means[~hvg], disp[~hvg], s=4, alpha=0.3, color="gray", label="not HVG")
    ax.scatter(means[hvg], disp[hvg], s=4, alpha=0.5, color="firebrick", label="HVG")
    ax.set_xscale("log")
    ax.set_xlabel("mean expression (log-normalized)")
    ax.set_ylabel("normalized dispersion")
    ax.set_title(f"HVG selection ({int(hvg.sum())} genes)")
    ax.legend(markerscale=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "hvg_mean_dispersion.png", dpi=150)
    plt.close(fig)
    print(f"Saved {PLOT_DIR / 'hvg_mean_dispersion.png'}")


def scale_and_subset(adata: sc.AnnData) -> sc.AnnData:
    """Stash full log-normalized matrix in .raw, subset to HVGs, then scale.

    Scaling (zero mean, unit variance, clipped at SCALE_MAX_VALUE) is done
    only on the HVG subset that Step 4/PCA will actually consume — no need
    to scale genes that won't enter the embedding.
    """
    adata.raw = adata
    hvg_adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg_adata, max_value=SCALE_MAX_VALUE)
    return hvg_adata


def normalize() -> sc.AnnData:
    adata = sc.read_h5ad(IN_PATH)
    n_cells_in, n_genes_in = adata.n_obs, adata.n_vars
    adata = normalize_and_log(adata)
    adata = select_hvgs(adata)
    plot_hvg_mean_variance(adata)
    scaled = scale_and_subset(adata)
    print(f"\nInput:  {n_cells_in} cells x {n_genes_in} genes ({IN_PATH.name})")
    print(f"Output: {scaled.n_obs} cells x {scaled.n_vars} genes (HVG-subset, scaled)")
    return scaled


if __name__ == "__main__":
    scaled = normalize()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    scaled.write_h5ad(OUT_PATH)
    print(f"\nSaved normalized/scaled data (Step 4 input) to {OUT_PATH}")
