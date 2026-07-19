"""Step 5: Cell-type classifier — the Week 1 go/no-go gate.

Train a supervised classifier on the Zenodo 6024273 reference's `Cell_type`
labels (label-transferred calls across 5 integrated public PDAC cohorts —
see src/reference_prep/), evaluate it against a naive marker-threshold
baseline on a held-out split of the *reference* (local data has no ground
truth), then apply the trained classifier to the local integrated cells.

Feature space: the intersection of reference genes and local genes (both
are log-normalized expression, comparable on that basis), restricted to the
top HVGs *within that intersection as computed on the reference* — the
reference's own feature-selection is what a "train on reference, apply
elsewhere" workflow should be keyed to, not the local Harmony-integration
HVG set (which was chosen for a different purpose in normalize.py).

Baseline: for each cell, z-score the 6 canonical markers (within the
held-out split) and assign the cell type whose marker is highest -- a
per-cell argmax over a fixed marker->cell-type mapping. It exists to be
beaten, not to be a strawman; results are reported even if the classifier
doesn't win.
"""

from pathlib import Path

import anndata as ad
import lightgbm as lgb
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import accuracy_score, f1_score, precision_recall_fscore_support
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

REPO_ROOT = Path(__file__).resolve().parent.parent
REF_PATH = REPO_ROOT / "data" / "reference" / "pk_all.h5ad"
LOCAL_PATH = REPO_ROOT / "data" / "processed" / "integrated.h5ad"
MODEL_DIR = REPO_ROOT / "data" / "models"
MODEL_PATH = MODEL_DIR / "cell_type_classifier.txt"
OUT_PATH = REPO_ROOT / "data" / "processed" / "classified.h5ad"
PRED_CSV = REPO_ROOT / "data" / "processed" / "local_cell_type_predictions.csv"

LABEL_COL = "Cell_type"
N_HVG_FEATURES = 2000
TEST_SIZE = 0.2
RANDOM_STATE = 0

# Canonical markers, per docs/agents/codebase-agent.md, extended to cover as
# many of the reference's 10 Cell_type classes as have a clean canonical
# marker. "Ductal cell type 1" vs "Ductal cell type 2" is deliberately left
# without a distinguishing marker here -- Peng et al. 2019 (Cell Res
# 29:725-738) report that ductal subtype 2 vs 1 is a malignancy/CNV-driven
# distinction, not a marker-gene one, so EPCAM (pan-ductal) is the only
# canonical marker available and it can't separate the two subtypes. This is
# a baseline limitation, not a fixable mapping gap.
MARKER_TO_CELL_TYPE = {
    "EPCAM": "Ductal cell type 1",
    "PTPRC": "T cell",
    "COL1A1": "Fibroblast cell",
    "PECAM1": "Endothelial cell",
    "CD3D": "T cell",
    "CD68": "Macrophage cell",
    "PRSS1": "Acinar cell",
    "INS": "Endocrine cell",
    "RGS5": "Stellate cell",
    "MS4A1": "B cell",
}


def load_reference() -> sc.AnnData:
    ref = sc.read_h5ad(REF_PATH)
    ref = ref[ref.obs[LABEL_COL].notna()].copy()
    print(f"Reference: {ref.n_obs} cells x {ref.n_vars} genes, "
          f"{ref.obs[LABEL_COL].nunique()} cell types")
    return ref


def select_feature_genes(ref_train: sc.AnnData, local: sc.AnnData) -> list[str]:
    """HVG selection sees only the train split of the reference (`ref_train`
    must already be subset to idx_train) so the held-out test cells never
    influence feature selection -- avoids leaking test-set variance into the
    feature space used to evaluate on that same test set.
    """
    local_genes = set(local.raw.var_names)
    common = [g for g in ref_train.var_names if g in local_genes]
    print(f"Genes in common between reference and local (raw, all genes): {len(common)}")

    ref_common = ref_train[:, common].copy()
    sc.pp.highly_variable_genes(ref_common, n_top_genes=N_HVG_FEATURES, flavor="seurat")
    feature_genes = ref_common.var_names[ref_common.var["highly_variable"]].tolist()

    for marker in MARKER_TO_CELL_TYPE:
        if marker in common and marker not in feature_genes:
            feature_genes.append(marker)
    print(f"Feature genes for classifier: {len(feature_genes)} "
          f"(reference HVGs within common set, markers force-included)")
    return feature_genes


def naive_baseline_predict(expr: pd.DataFrame) -> np.ndarray:
    """Per-cell argmax over z-scored canonical marker expression, mapped to
    the coarse cell type each marker denotes. Cells default to whichever
    marker has the highest z-score even at low absolute expression -- this
    is the "naive" part; it is not meant to abstain.
    """
    markers = [m for m in MARKER_TO_CELL_TYPE if m in expr.columns]
    z = (expr[markers] - expr[markers].mean()) / expr[markers].std().replace(0, 1)
    winning_marker = z.idxmax(axis=1)
    return winning_marker.map(MARKER_TO_CELL_TYPE).values


def split_reference(ref: sc.AnnData) -> tuple[np.ndarray, np.ndarray, np.ndarray, LabelEncoder]:
    """Stratified train/test split on cell indices, computed before any
    feature selection touches the reference -- so HVG selection can be
    restricted to the train split only and never sees held-out cells.
    """
    y_raw = ref.obs[LABEL_COL].astype(str).values
    encoder = LabelEncoder()
    y = encoder.fit_transform(y_raw)
    idx_train, idx_test = train_test_split(
        np.arange(len(y)), test_size=TEST_SIZE, random_state=RANDOM_STATE, stratify=y
    )
    return idx_train, idx_test, y, encoder


def train_and_evaluate(ref: sc.AnnData, feature_genes: list[str], idx_train: np.ndarray,
                        idx_test: np.ndarray, y: np.ndarray,
                        encoder: LabelEncoder) -> lgb.Booster:
    X = ref[:, feature_genes].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)

    X_train, X_test = X[idx_train], X[idx_test]
    y_train, y_test = y[idx_train], y[idx_test]

    train_set = lgb.Dataset(X_train, label=y_train)
    test_set = lgb.Dataset(X_test, label=y_test, reference=train_set)
    params = {
        "objective": "multiclass",
        "num_class": len(encoder.classes_),
        "metric": "multi_logloss",
        "verbosity": -1,
        "seed": RANDOM_STATE,
    }
    booster = lgb.train(
        params, train_set, num_boost_round=200,
        valid_sets=[test_set], callbacks=[lgb.early_stopping(20), lgb.log_evaluation(0)],
    )

    y_pred_proba = booster.predict(X_test, num_iteration=booster.best_iteration)
    y_pred = np.argmax(y_pred_proba, axis=1)
    clf_acc = accuracy_score(y_test, y_pred)
    clf_f1 = f1_score(y_test, y_pred, average="macro")

    expr_test_df = pd.DataFrame(X_test, columns=feature_genes)
    baseline_pred_labels = naive_baseline_predict(expr_test_df)
    known_labels = set(encoder.classes_)
    baseline_pred_labels = np.array(
        [lbl if lbl in known_labels else "__unmapped__" for lbl in baseline_pred_labels]
    )
    y_test_labels = encoder.inverse_transform(y_test)
    baseline_acc = accuracy_score(y_test_labels, baseline_pred_labels)
    baseline_f1 = f1_score(y_test_labels, baseline_pred_labels, average="macro")

    print(f"\nHeld-out split: {len(y_train)} train / {len(y_test)} test cells "
          f"(stratified {TEST_SIZE:.0%}, random_state={RANDOM_STATE})")
    print(f"\n{'':>20s} {'accuracy':>10s} {'macro-F1':>10s}")
    print(f"{'LightGBM classifier':>20s} {clf_acc:10.4f} {clf_f1:10.4f}")
    print(f"{'naive marker baseline':>20s} {baseline_acc:10.4f} {baseline_f1:10.4f}")
    verdict = "BEATS" if (clf_acc > baseline_acc and clf_f1 > baseline_f1) else "DOES NOT BEAT"
    print(f"\nGo/no-go: classifier {verdict} the naive baseline on held-out reference data.")

    print(f"\nPer-class LightGBM classifier performance (held-out test, {len(encoder.classes_)} classes):")
    clf_precision, clf_recall, clf_per_class_f1, clf_support = precision_recall_fscore_support(
        y_test, y_pred, labels=range(len(encoder.classes_)), zero_division=0
    )
    for cls, p, r, f1, n in zip(encoder.classes_, clf_precision, clf_recall, clf_per_class_f1, clf_support):
        print(f"  {cls:>22s}  precision={p:.4f}  recall={r:.4f}  f1={f1:.4f}  n={n}")

    # Training-set class ratio for Ductal1 vs Ductal2, to distinguish a
    # genuine reference class imbalance from a classifier detection-
    # sensitivity issue when interpreting the local Ductal1:Ductal2 skew.
    train_labels = encoder.inverse_transform(y_train)
    d1_train = int((train_labels == "Ductal cell type 1").sum())
    d2_train = int((train_labels == "Ductal cell type 2").sum())
    d1_idx = list(encoder.classes_).index("Ductal cell type 1")
    d1_recall = clf_recall[d1_idx]
    print(f"\nDuctal1 detection-sensitivity check: held-out recall={d1_recall:.4f} "
          f"(other classes range {clf_recall.min():.4f}-{clf_recall.max():.4f}); "
          f"training-set Ductal1:Ductal2 = {d1_train}:{d2_train} ({d1_train / (d1_train + d2_train):.1%} Ductal1)")

    print(f"\nPer-class baseline performance ({len(known_labels)} classes in reference):")
    per_class_f1 = f1_score(y_test_labels, baseline_pred_labels, average=None, labels=sorted(known_labels))
    per_class_acc = pd.Series(y_test_labels == baseline_pred_labels, index=y_test_labels).groupby(level=0).mean()
    for cls, f1 in zip(sorted(known_labels), per_class_f1):
        print(f"  {cls:>22s}  acc={per_class_acc.get(cls, float('nan')):.4f}  f1={f1:.4f}"
              f"{'  (unmapped in baseline)' if cls not in MARKER_TO_CELL_TYPE.values() else ''}")

    return booster


def apply_to_local(booster: lgb.Booster, encoder: LabelEncoder, local: sc.AnnData,
                    feature_genes: list[str]) -> pd.Series:
    X_local = local.raw[:, feature_genes].X
    X_local = np.asarray(X_local.todense()) if hasattr(X_local, "todense") else np.asarray(X_local)
    proba = booster.predict(X_local, num_iteration=booster.best_iteration)
    pred_idx = np.argmax(proba, axis=1)
    pred_labels = encoder.inverse_transform(pred_idx)
    # `predicted_cell_type == "Ductal cell type 2"` is the intended malignant-cell
    # filter for Step 7 (subtype scoring): Peng et al. 2019 (Cell Res 29:725-738)
    # identify ductal subtype 2 as the malignant/CNV-bearing ductal population,
    # vs. subtype 1 as normal ductal -- this label, not a marker threshold, is
    # what should gate malignant-cell selection downstream.
    predictions = pd.Series(pred_labels, index=local.obs_names, name="predicted_cell_type")

    print(f"\nApplied classifier to {local.n_obs} local cells.")
    print(predictions.value_counts())

    print("\nPer-sample predicted cell-type composition (row-normalized %):")
    comp = pd.crosstab(local.obs["sample"], predictions, normalize="index") * 100
    print(comp.round(1))
    return predictions


CANONICAL_CHECK_MARKERS = {
    "PTPRC": "T cell",
    "COL1A1": "Fibroblast cell",
    "PECAM1": "Endothelial cell",
    "CD3D": "T cell",
    "CD68": "Macrophage cell",
}

# EPCAM deliberately excluded from the per-marker check above: it's
# pan-ductal/pan-epithelial by definition (Ductal1 and Ductal2 are both
# epithelial), so it cannot discriminate the two classes and passing/failing
# on it says nothing about whether the classifier's Ductal1-vs-Ductal2 split
# is correct. Peng et al. 2019 (Cell Res 29:725-738) distinguish Ductal2
# (malignant, CNV-defined) from Ductal1 (normal-like) primarily via CNV
# inference, but they also report, from their own DE/text (not this repo's
# assumption -- standard proliferation genes like MKI67/TOP2A turned out to
# be specific to a further malignant subgroup, not the Ductal1-vs-2 split
# itself) that type 2 cells show "much higher expression of reported poor
# prognosis PDAC markers, such as CEACAM1/5/6 and KRT19" relative to type 1.
MALIGNANT_DUCTAL_GENES = ["CEACAM1", "CEACAM5", "CEACAM6", "KRT19"]


def face_validity_check(local: sc.AnnData, predictions: pd.Series) -> None:
    """Sanity check on the local predictions themselves (not the reference):
    for each of the canonical markers, do cells the classifier assigned to
    that marker's cell type actually show elevated expression of it in the
    local data? Catches a wholesale domain-shift failure (e.g. FNA-specific
    technical differences causing systematic mislabeling) before Step 6
    aggregates per-patient composition off these labels.
    """
    print("\nFace-validity check: local marker expression by predicted cell type")
    markers = [m for m in CANONICAL_CHECK_MARKERS if m in local.raw.var_names]
    expr = local.raw[:, markers].X
    expr = np.asarray(expr.todense()) if hasattr(expr, "todense") else np.asarray(expr)
    expr_df = pd.DataFrame(expr, columns=markers, index=local.obs_names)
    expr_df["predicted_cell_type"] = predictions.values

    for marker, target_type in CANONICAL_CHECK_MARKERS.items():
        if marker not in markers:
            print(f"  {marker}: not found in local data, skipping")
            continue
        by_type = expr_df.groupby("predicted_cell_type")[marker].mean().sort_values(ascending=False)
        target_mean = by_type.get(target_type, float("nan"))
        top_type = by_type.index[0]
        flag = "" if top_type == target_type else "  *** target type is not top expressor ***"
        print(f"  {marker} -> {target_type}: mean expr in target = {target_mean:.3f}, "
              f"highest-expressing predicted type = {top_type} ({by_type.iloc[0]:.3f}){flag}")

    ductal1_vs_ductal2_check(local, predictions)


def ductal1_vs_ductal2_check(local: sc.AnnData, predictions: pd.Series) -> None:
    """Replaces the (invalid) EPCAM leg of the face-validity check. Scores
    Peng et al.'s own malignant-ductal gene panel (CEACAM1/5/6, KRT19 --
    pulled from the paper's text, see comment above) in locally-predicted
    Ductal cell type 2 vs Ductal cell type 1 cells. If the classifier
    recovered the malignant/normal-like split correctly, predicted Ductal2
    cells should show higher expression of this panel than predicted Ductal1
    cells in the local data itself.
    """
    genes = [g for g in MALIGNANT_DUCTAL_GENES if g in local.raw.var_names]
    print(f"\nDuctal cell type 1 vs 2 check (Peng et al. 2019 malignant-ductal panel: {', '.join(genes)}):")
    if not genes:
        print("  none of the panel genes found in local data -- cannot run this check")
        return

    expr = local.raw[:, genes].X
    expr = np.asarray(expr.todense()) if hasattr(expr, "todense") else np.asarray(expr)
    expr_df = pd.DataFrame(expr, columns=genes, index=local.obs_names)
    expr_df["predicted_cell_type"] = predictions.values

    d1 = expr_df[expr_df["predicted_cell_type"] == "Ductal cell type 1"]
    d2 = expr_df[expr_df["predicted_cell_type"] == "Ductal cell type 2"]
    print(f"  predicted Ductal1 n={len(d1)}, predicted Ductal2 n={len(d2)}")
    for gene in genes:
        m1, m2 = d1[gene].mean(), d2[gene].mean()
        flag = "" if m2 > m1 else "  *** Ductal2 does NOT exceed Ductal1 ***"
        print(f"  {gene}: mean(Ductal1)={m1:.3f}  mean(Ductal2)={m2:.3f}{flag}")


def main() -> None:
    ref = load_reference()
    local = ad.read_h5ad(LOCAL_PATH)

    idx_train, idx_test, y, encoder = split_reference(ref)
    feature_genes = select_feature_genes(ref[idx_train], local)
    booster = train_and_evaluate(ref, feature_genes, idx_train, idx_test, y, encoder)

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    booster.save_model(str(MODEL_PATH))
    print(f"\nSaved model to {MODEL_PATH}")

    predictions = apply_to_local(booster, encoder, local, feature_genes)
    face_validity_check(local, predictions)
    local.obs["predicted_cell_type"] = predictions
    local.write_h5ad(OUT_PATH)
    predictions.to_frame().join(local.obs["sample"]).to_csv(PRED_CSV)
    print(f"Saved classified AnnData to {OUT_PATH}")
    print(f"Saved per-cell predictions to {PRED_CSV}")


if __name__ == "__main__":
    main()
