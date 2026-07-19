"""
Step 6: aggregate per-cell classifier predictions into a per-sample TME
composition breakdown (counts + proportions per cell type, one row per sample).
"""

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
PRED_CSV = REPO_ROOT / "data" / "processed" / "local_cell_type_predictions.csv"
OUT_CSV = REPO_ROOT / "data" / "processed" / "composition.csv"

# Fixed order matching the Zenodo 6024273 reference's 10-category coarse
# cell-type taxonomy (see classify.py / design doc reference-dataset note).
CELL_TYPES = [
    "Ductal cell type 1",
    "Ductal cell type 2",
    "Acinar cell",
    "Endocrine cell",
    "Fibroblast cell",
    "Stellate cell",
    "Endothelial cell",
    "Macrophage cell",
    "T cell",
    "B cell",
]

# Below this count, a cell type's proportion for a sample is too noisy to
# trust at face value (small-N denominator effect, not a modeling issue).
LOW_COUNT_FLAG_THRESHOLD = 20


def build_composition(pred_csv: Path = PRED_CSV) -> pd.DataFrame:
    preds = pd.read_csv(pred_csv)
    counts = (
        preds.groupby(["sample", "predicted_cell_type"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=CELL_TYPES, fill_value=0)
        .sort_index()
    )
    totals = counts.sum(axis=1)
    proportions = counts.div(totals, axis=0)

    out = pd.DataFrame(index=counts.index)
    out["total_cells"] = totals
    for ct in CELL_TYPES:
        out[f"count_{ct}"] = counts[ct]
    for ct in CELL_TYPES:
        out[f"prop_{ct}"] = proportions[ct]

    return out.reset_index().rename(columns={"sample": "sample"})


def flag_low_counts(counts_wide: pd.DataFrame) -> list[str]:
    """Return human-readable notes for (sample, cell_type) cells with counts
    below LOW_COUNT_FLAG_THRESHOLD, i.e. proportions worth treating as noisy."""
    notes = []
    for _, row in counts_wide.iterrows():
        sample = row["sample"]
        for ct in CELL_TYPES:
            n = row[f"count_{ct}"]
            if 0 < n < LOW_COUNT_FLAG_THRESHOLD:
                prop = row[f"prop_{ct}"]
                notes.append(
                    f"{sample}: {ct} n={n} ({prop:.1%} of {row['total_cells']} cells) — low count, noisy proportion"
                )
    return notes


def main():
    composition = build_composition()
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    composition.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV} ({len(composition)} samples x {len(CELL_TYPES)} cell types)")

    print("\nPer-sample composition (%):")
    pct_cols = ["sample"] + [f"prop_{ct}" for ct in CELL_TYPES]
    pct_view = composition[pct_cols].copy()
    for ct in CELL_TYPES:
        pct_view[f"prop_{ct}"] = (pct_view[f"prop_{ct}"] * 100).round(1)
    pct_view.columns = ["sample"] + CELL_TYPES
    print(pct_view.to_string(index=False))

    low_count_notes = flag_low_counts(composition)
    if low_count_notes:
        print(f"\nLow-count flags (n < {LOW_COUNT_FLAG_THRESHOLD}, proportion likely noisy):")
        for note in low_count_notes:
            print(f"  - {note}")


if __name__ == "__main__":
    main()
