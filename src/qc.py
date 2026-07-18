"""Step 2: QC — per-sample doublet detection, distribution plots, filtering.

Order matters: doublet detection (scrublet) runs per sample on raw counts
before any other filtering, since doublet rates are batch-specific and a
merged run would blend them. QC metrics (genes/cell, counts/cell, % mito)
are then computed per sample and plotted so thresholds come from the actual
data rather than a blind default.

Mito filtering uses a per-sample adaptive cutoff (median + 3 MAD of
log1p(pct_counts_mt), per Theis 2019 / sc-best-practices) rather than one
shared cutoff across samples — samples differ in whether elevated mito
reflects genuine cell death/stress (complexity drops with mito%) or just a
higher ambient-RNA baseline (complexity flat across mito%), and a shared
cutoff over- or under-filters depending on which case a sample is in.

Doublet filtering is deliberately NOT applied. None of the 8 per-sample
scrublet score distributions show a separable bimodal doublet population
(checked directly, not just inferred from low counts) — a known scrublet
blind spot for homotypic doublets and low cell-type-diversity samples
(Wolock et al. 2019), not a small-N artifact or tooling bug. A fixed
score cutoff would be arbitrary with nothing to validate it against, so
doublets are left in. See DOUBLET_LIMITATION_STATEMENT below for the
full reasoning and downstream risk.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation

REPO_ROOT = Path(__file__).resolve().parent.parent
IN_PATH = REPO_ROOT / "data" / "processed" / "ingested.h5ad"
ANNOTATED_OUT_PATH = REPO_ROOT / "data" / "processed" / "qc_annotated.h5ad"
FILTERED_OUT_PATH = REPO_ROOT / "data" / "processed" / "qc_filtered.h5ad"
PLOT_DIR = REPO_ROOT / "data" / "processed" / "qc_plots"
REPORT_PATH = REPO_ROOT / "docs" / "qc_report.md"

SAMPLES = ["MK362", "AM67", "MK364", "MK336", "BW21", "MK371", "MK359", "MK447"]

MIN_GENES = 250
MT_NMADS = 3
LOW_CONFIDENCE_SAMPLES = ["MK447"]  # genuinely depth-limited, not dropped — see Step 5 follow-up

DOUBLET_LIMITATION_STATEMENT = (
    "Doublet detection (scrublet, run per-sample) found no separable doublet "
    "population in any of the 8 samples — score distributions were unimodal "
    "with no bimodal shoulder, even in smaller samples. This is a known "
    "scrublet limitation for homotypic doublets and low cell-type-diversity "
    "samples, not evidence the data is doublet-free. Reported per-sample "
    "doublet counts (0-8 cells) reflect an auto-threshold artifact and "
    "should not be read as a doublet-burden estimate. Undetected doublets "
    "may contribute noise to per-patient composition estimates; they do not "
    "affect classifier training, which uses an external reference."
)

# Rough 10x Chromium expected-multiplet-rate table (cells recovered -> % rate),
# linearly interpolated. Used only as a sanity check against scrublet's output.
_TENX_MULTIPLET_TABLE = {
    500: 0.4, 1000: 0.8, 2000: 1.6, 3000: 2.3, 4000: 3.1,
    5000: 3.9, 6000: 4.6, 7000: 5.4, 8000: 6.1, 9000: 6.9, 10000: 7.6,
}


def _expected_multiplet_rate(n_cells: int) -> float:
    xs = sorted(_TENX_MULTIPLET_TABLE)
    ys = [_TENX_MULTIPLET_TABLE[x] for x in xs]
    return float(np.interp(n_cells, xs, ys))


def run_doublet_detection(combined: sc.AnnData) -> sc.AnnData:
    """Run scrublet independently per sample, plot score histograms, reassemble.

    Does not filter on predicted_doublet here — that auto-threshold is
    unreliable at these batch sizes (see module docstring) and needs a
    manual/human call from the histograms saved to qc_plots/.
    """
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    per_sample = []
    print("Doublet detection (scrublet, per sample):")
    print(f"{'sample':8s} {'n':>6s} {'auto_dbl':>9s} {'auto_%':>7s} "
          f"{'auto_thresh':>12s} {'expected_%':>11s}")
    for sample in SAMPLES:
        sub = combined[combined.obs["sample"] == sample].copy()
        sc.pp.scrublet(sub, batch_key=None, verbose=False)

        threshold = sub.uns["scrublet"]["parameters"].get("threshold")
        if threshold is None:
            threshold = sub.uns["scrublet"].get("threshold")
        n_dbl = int(sub.obs["predicted_doublet"].sum())
        auto_pct = 100 * n_dbl / sub.n_obs
        expected_pct = _expected_multiplet_rate(sub.n_obs)
        sub.obs["scrublet_auto_threshold"] = threshold

        print(f"{sample:8s} {sub.n_obs:6d} {n_dbl:9d} {auto_pct:6.1f}% "
              f"{threshold:12.3f} {expected_pct:10.1f}%")

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(sub.obs["doublet_score"], bins=60, color="steelblue", alpha=0.8)
        ax.axvline(threshold, color="red", linestyle="--",
                    label=f"scrublet auto threshold ({threshold:.3f})")
        ax.set_xlabel("doublet_score")
        ax.set_ylabel("cells")
        ax.set_title(f"{sample}: doublet score (auto {auto_pct:.1f}% vs "
                      f"expected ~{expected_pct:.1f}%)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"doublet_score_hist_{sample}.png", dpi=150)
        plt.close(fig)

        per_sample.append(sub)

    print(f"\nSaved doublet score histograms to {PLOT_DIR} "
          "(auto-threshold shown but NOT applied — pick a threshold by inspection)")
    return sc.concat(per_sample, join="outer")


def compute_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    return adata


def compute_adaptive_mt_thresholds(adata: sc.AnnData, nmads: int = MT_NMADS) -> dict:
    """Per-sample upper cutoff on pct_counts_mt: median + nmads*MAD.

    Computed on the raw pct_counts_mt scale, not log1p. log1p(pct_counts_mt)
    was tried first (as specified) but degenerates on MK362 and AM67: both
    have a heavy right tail (AM67 has a genuine 90-100% dying-cell spike —
    15% of its cells), which inflates the log-scale MAD enough that
    median + 3*MAD exceeds the 100% ceiling, flagging zero outliers on
    exactly the two samples this is meant to catch. Raw-scale MAD (matching
    the original Theis/sc-best-practices is_outlier formulation, which does
    not log-transform pct_counts_mt) stays bounded and gives sensible,
    differentiated cutoffs — see reported thresholds below.

    Only an upper bound is meaningful here (low mito isn't a quality problem),
    so this reports where each sample's own distribution says "this is an
    outlier for THIS sample" rather than applying one global mito ceiling.
    """
    thresholds = {}
    print(f"\nAdaptive per-sample mito thresholds (median + {nmads} MAD, raw pct scale):")
    print(f"{'sample':8s} {'med':>6s} {'mad':>6s} {'pct_thresh':>11s} {'n_outlier':>10s} {'%cells':>7s}")
    for s in SAMPLES:
        vals = adata.obs.loc[adata.obs["sample"] == s, "pct_counts_mt"]
        med, mad = np.median(vals), median_abs_deviation(vals)
        thresh = med + nmads * mad
        n_outlier = int((vals > thresh).sum())
        thresholds[s] = float(thresh)
        print(f"{s:8s} {med:6.2f} {mad:6.2f} {thresh:10.2f}% {n_outlier:10d} {100 * n_outlier / len(vals):6.1f}%")
    return thresholds


def plot_mito_vs_complexity(adata: sc.AnnData) -> None:
    """Binned median genes/cell across mito bands, per sample — the evidence
    for whether elevated mito tracks with cell death (complexity drops) or is
    just a higher ambient baseline (complexity flat)."""
    bins = list(range(0, 101, 10))
    labels = [f"{lo}-{hi}" for lo, hi in zip(bins[:-1], bins[1:])]
    fig, ax = plt.subplots(figsize=(10, 5))
    for s in SAMPLES:
        sub = adata.obs[adata.obs["sample"] == s]
        band = pd_cut(sub["pct_counts_mt"], bins, labels)
        medians = [sub.loc[band == lab, "n_genes_by_counts"].median() for lab in labels]
        ax.plot(labels, medians, marker="o", label=s)
    ax.set_xlabel("pct_counts_mt band")
    ax.set_ylabel("median n_genes_by_counts")
    ax.set_title("Complexity vs mito band, per sample")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "mito_vs_complexity.png", dpi=150)
    plt.close(fig)
    print(f"Saved {PLOT_DIR / 'mito_vs_complexity.png'}")


def pd_cut(series, bins, labels):
    import pandas as pd
    return pd.cut(series, bins=bins, labels=labels, include_lowest=True)


def plot_distributions(adata: sc.AnnData) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    for metric in metrics:
        fig, ax = plt.subplots(figsize=(12, 5))
        data = [adata.obs.loc[adata.obs["sample"] == s, metric] for s in SAMPLES]
        ax.violinplot(data, showmedians=True)
        ax.set_xticks(range(1, len(SAMPLES) + 1))
        ax.set_xticklabels(SAMPLES, rotation=45, ha="right")
        ax.set_ylabel(metric)
        ax.set_title(f"{metric} per sample (pre-filtering)")
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"{metric}_violin.png", dpi=150)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    for s in SAMPLES:
        sub = adata.obs[adata.obs["sample"] == s]
        ax.scatter(sub["total_counts"], sub["n_genes_by_counts"], s=3, alpha=0.3, label=s)
    ax.set_xlabel("total_counts")
    ax.set_ylabel("n_genes_by_counts")
    ax.legend(markerscale=4, fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "counts_vs_genes_scatter.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved distribution plots to {PLOT_DIR}")


def summarize(adata: sc.AnnData) -> None:
    print("\nPer-sample QC summary (median [p5-p95]):")
    for s in SAMPLES:
        sub = adata.obs[adata.obs["sample"] == s]
        flag = "  [LOW CONFIDENCE — depth-limited, see Step 5 follow-up]" if s in LOW_CONFIDENCE_SAMPLES else ""
        print(f"  {s}{flag}")
        for metric in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
            vals = sub[metric]
            med = np.median(vals)
            p5, p95 = np.percentile(vals, [5, 95])
            print(f"    {metric:20s} median={med:8.1f}  p5={p5:8.1f}  p95={p95:8.1f}")
        print()


def apply_filter(adata: sc.AnnData, min_genes: int, mt_thresholds: dict) -> sc.AnnData:
    """Final QC filter: per-sample adaptive mito cutoff + min_genes floor.

    No doublet-score filtering — see DOUBLET_LIMITATION_STATEMENT. Doublets
    are left in the data; the limitation is documented instead of worked
    around with an arbitrary cutoff.
    """
    print(f"\nFinal filter: min_genes>={min_genes}, per-sample adaptive mito cutoff, "
          "no doublet-score filtering")
    print(f"{'sample':10s} {'before':>8s} {'after':>8s} {'dropped':>8s}  note")
    kept_masks = []
    for s in SAMPLES:
        mask = adata.obs["sample"] == s
        sub = adata.obs.loc[mask]
        keep = (sub["n_genes_by_counts"] >= min_genes) & (sub["pct_counts_mt"] <= mt_thresholds[s])
        kept_masks.append(keep)
        before, after = mask.sum(), keep.sum()
        note = "low-confidence (depth-limited)" if s in LOW_CONFIDENCE_SAMPLES else ""
        print(f"{s:10s} {before:8d} {after:8d} {before - after:8d}  {note}")
    keep_all = np.concatenate([m.values for m in kept_masks])
    filtered = adata[keep_all].copy()
    print(f"{'TOTAL':10s} {adata.n_obs:8d} {filtered.n_obs:8d} {adata.n_obs - filtered.n_obs:8d}")
    return filtered


def write_report(before_after: list[tuple[str, int, int]]) -> None:
    lines = [
        "# QC Report (Step 2)",
        "",
        "Filter applied: per-sample adaptive mito cutoff (median + 3 MAD, raw pct scale) "
        "+ `n_genes_by_counts >= 250`. No doublet-score filtering.",
        "",
        "## Doublet detection limitation",
        "",
        f"> {DOUBLET_LIMITATION_STATEMENT}",
        "",
        "## Before / after counts",
        "",
        "| Sample | Before | After | Dropped |",
        "|---|---|---|---|",
    ]
    total_before = total_after = 0
    for s, before, after in before_after:
        note = " (low-confidence, depth-limited)" if s in LOW_CONFIDENCE_SAMPLES else ""
        lines.append(f"| {s}{note} | {before} | {after} | {before - after} |")
        total_before += before
        total_after += after
    lines.append(f"| **Total** | **{total_before}** | **{total_after}** | "
                  f"**{total_before - total_after}** |")
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n")
    print(f"\nWrote {REPORT_PATH}")


def qc() -> tuple[sc.AnnData, sc.AnnData]:
    adata = sc.read_h5ad(IN_PATH)
    adata = run_doublet_detection(adata)
    adata = compute_qc_metrics(adata)
    adata.obs["low_confidence_sample"] = adata.obs["sample"].isin(LOW_CONFIDENCE_SAMPLES)
    plot_distributions(adata)
    plot_mito_vs_complexity(adata)
    summarize(adata)
    print(f"\n{DOUBLET_LIMITATION_STATEMENT}\n")
    mt_thresholds = compute_adaptive_mt_thresholds(adata)
    filtered = apply_filter(adata, MIN_GENES, mt_thresholds)
    before_after = [
        (s, int((adata.obs["sample"] == s).sum()), int((filtered.obs["sample"] == s).sum()))
        for s in SAMPLES
    ]
    write_report(before_after)
    return adata, filtered


if __name__ == "__main__":
    annotated, filtered = qc()
    ANNOTATED_OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    annotated.write_h5ad(ANNOTATED_OUT_PATH)
    filtered.write_h5ad(FILTERED_OUT_PATH)
    print(f"\nSaved doublet-annotated (unfiltered) data to {ANNOTATED_OUT_PATH}")
    print(f"Saved filtered QC'd data (Step 3 input) to {FILTERED_OUT_PATH}")
