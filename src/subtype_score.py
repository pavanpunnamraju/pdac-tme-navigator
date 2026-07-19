"""Step 7: Moffitt classical/basal-like subtype scoring.

Spec: docs/step_7_subtype_scoring_spec.md (logic-agent). Scores the 25-gene
classical and 25-gene basal Moffitt signatures (Moffitt RA et al., Nature
Genetics 2015, PMID 26343385, via the pdacR package) on Step 3's normalized
log-expression data, restricted to cells the Step 5 classifier called
"Ductal cell type 2" (the malignant/CNV-bearing ductal population per Peng
et al. 2019) -- the only population a classical/basal identity is meaningful
for.

Scoring uses `adata.raw` from normalized.h5ad (full-gene-set, shifted-log
expression) rather than normalized.h5ad's top-level `.X`, which is scaled and
restricted to the 2000 Harmony-integration HVGs -- most of the 50 signature
genes aren't in that HVG set at all, and scaled values are the wrong input
for score_genes' control-gene-set correction.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parent.parent
NORMALIZED_PATH = REPO_ROOT / "data" / "processed" / "normalized.h5ad"
PRED_CSV = REPO_ROOT / "data" / "processed" / "local_cell_type_predictions.csv"
PER_CELL_OUT = REPO_ROOT / "data" / "processed" / "subtype_scores_per_cell.csv"
PER_PATIENT_OUT = REPO_ROOT / "data" / "processed" / "subtype_scores.csv"
GENE_RESOLUTION_OUT = REPO_ROOT / "data" / "processed" / "gene_resolution.json"

DUCTAL2_LABEL = "Ductal cell type 2"

BASAL_GENES = [
    "ANXA8L2", "AREG", "CST6", "CTSL2", "DHRS9", "FAM83A", "FGFBP1", "GPR87",
    "KRT15", "KRT17", "KRT6A", "KRT6C", "KRT7", "LEMD1", "LY6D", "S100A2",
    "SCEL", "SERPINB3", "SERPINB4", "SLC2A1", "SPRR1B", "SPRR3", "TNS4",
    "UCA1", "VGLL1",
]

CLASSICAL_GENES = [
    "AGR2", "AGR3", "ANXA10", "ATAD4", "BTNL8", "CDH17", "CEACAM6", "CLRN3",
    "CTSE", "CYP3A7", "FAM3D", "KRT20", "LGALS4", "LOC400573", "LYZ", "MYO1A",
    "PLA2G10", "REG4", "SPINK4", "ST6GALNAC1", "TFF1", "TFF2", "TFF3",
    "TSPAN8", "VSIG2",
]

# Current-symbol remapping for 2015-vintage names, checked against
# adata.raw.var_names (see docs/step_7_subtype_scoring_spec.md section 1).
# CTSL2 -> CTSV and ATAD4 -> PRR15L are confirmed HGNC renamings (both
# resolve to genes present in this reference). LOC400573 (an Entrez
# placeholder, not a real symbol) resolves to TMEM238L per NCBI, but
# TMEM238L is not present in this 10x reference's gene set either, so it
# still drops. ANXA8L2 -> ANXA8L1 is a tentative paralog resolution, not a
# strict renaming -- GENCODE versions disagree on which of the two
# near-identical ANXA8 paralogs is annotated, so this mapping is flagged as
# lower-confidence in the report below. UCA1 has no resolvable current
# symbol in this reference's annotation (checked LINC00178, UCAT1, CUDR).
SYMBOL_REMAP = {
    "CTSL2": "CTSV",
    "ANXA8L2": "ANXA8L1",
    "ATAD4": "PRR15L",
    "LOC400573": "TMEM238L",
}
TENTATIVE_REMAPS = {"ANXA8L2"}

# Reasons for genes that remain unresolved after the remap table above --
# recorded so the provenance JSON explains *why*, not just *that*, a gene
# was dropped.
DROP_REASONS = {
    "LOC400573": "Entrez placeholder resolves to TMEM238L (per NCBI), but "
                  "TMEM238L is not present in this reference's gene set either",
    "UCA1": "no resolvable current symbol found in this reference's "
            "annotation (checked LINC00178, UCAT1, CUDR)",
}

# n >= 43 floor: 95% CI half-width on a binomial proportion at worst-case
# 50/50 is SE = sqrt(0.25/n); requiring half-width <= 15pp solves to
# n >= 1.96^2 * 0.25 / 0.15^2 ~= 43 (spec section 4).
MIN_CELLS_FOR_CALL = 43
# Secondary, same-derivation threshold for "clears the floor but barely":
# half-width <= 10pp solves to n >= 1.96^2 * 0.25 / 0.10^2 ~= 97. Below this,
# flag the per-patient call as low-confidence rather than reporting it at
# the same confidence level as the comfortably-powered samples.
LOW_CONFIDENCE_CELL_THRESHOLD = 97


def resolve_signature(genes: list[str], var_names: set[str], label: str) -> tuple[list[str], dict]:
    resolved, matched_direct, remapped, dropped = [], [], [], []
    for gene in genes:
        if gene in var_names:
            resolved.append(gene)
            matched_direct.append(gene)
        elif gene in SYMBOL_REMAP and SYMBOL_REMAP[gene] in var_names:
            mapped = SYMBOL_REMAP[gene]
            resolved.append(mapped)
            remapped.append({
                "original": gene,
                "remapped_to": mapped,
                "confidence": "tentative" if gene in TENTATIVE_REMAPS else "confident",
            })
        else:
            dropped.append({
                "gene": gene,
                "reason": DROP_REASONS.get(gene, "no current symbol found in adata.raw.var_names"),
            })

    print(f"\n{label} signature: {len(genes)} genes total, {len(resolved)} matched, "
          f"{len(dropped)} dropped")
    if remapped:
        remap_strs = []
        for r in remapped:
            suffix = " (tentative)" if r["confidence"] == "tentative" else ""
            remap_strs.append(f"{r['original']}->{r['remapped_to']}{suffix}")
        print(f"  remapped: {', '.join(remap_strs)}")
    if dropped:
        print(f"  dropped (no current symbol found in adata.raw.var_names): "
              f"{', '.join(d['gene'] for d in dropped)}")

    record = {
        "total_genes": len(genes),
        "n_matched": len(resolved),
        "n_dropped": len(dropped),
        "matched_direct": matched_direct,
        "remapped": remapped,
        "dropped": dropped,
    }
    return resolved, record


def load_ductal2_expression() -> sc.AnnData:
    norm = sc.read_h5ad(NORMALIZED_PATH)
    expr = norm.raw.to_adata()
    expr.obs["sample"] = norm.obs["sample"]

    preds = pd.read_csv(PRED_CSV, index_col=0)
    ductal2_barcodes = preds.index[preds["predicted_cell_type"] == DUCTAL2_LABEL]
    missing = ductal2_barcodes.difference(expr.obs_names)
    if len(missing) > 0:
        print(f"WARNING: {len(missing)} Ductal2 barcodes from {PRED_CSV.name} "
              f"not found in {NORMALIZED_PATH.name}, dropping them")
    ductal2_barcodes = ductal2_barcodes.intersection(expr.obs_names)

    subset = expr[ductal2_barcodes].copy()
    print(f"\nDuctal2 cells for scoring: {subset.n_obs} "
          f"({subset.n_obs} of {len(preds)} total classified local cells)")
    return subset


def score_and_call(adata: sc.AnnData) -> tuple[pd.DataFrame, dict]:
    var_names = set(adata.var_names)
    classical_genes, classical_record = resolve_signature(CLASSICAL_GENES, var_names, "Classical")
    basal_genes, basal_record = resolve_signature(BASAL_GENES, var_names, "Basal")

    sc.tl.score_genes(adata, gene_list=classical_genes, score_name="classical_score")
    sc.tl.score_genes(adata, gene_list=basal_genes, score_name="basal_score")

    result = adata.obs[["sample", "classical_score", "basal_score"]].copy()
    result["score_diff"] = result["classical_score"] - result["basal_score"]
    result["call"] = np.where(result["score_diff"] > 0, "classical", "basal")

    resolution = {"classical": classical_record, "basal": basal_record}
    return result, resolution


def report_score_diff_distribution(per_cell: pd.DataFrame) -> None:
    diff = per_cell["score_diff"]
    print("\nPer-cell classical-minus-basal score-difference distribution "
          "(no ambiguity threshold applied; reported for reviewer visibility):")
    print(diff.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).to_string())


def per_patient_rollup(per_cell: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, group in per_cell.groupby("sample", observed=True):
        n = len(group)
        n_classical = int((group["call"] == "classical").sum())
        pct_classical = n_classical / n * 100

        if n < MIN_CELLS_FOR_CALL:
            majority_call = "insufficient cells, no call"
            pct_display = np.nan
            confidence = "n/a"
        else:
            majority_call = "classical" if n_classical > n - n_classical else "basal"
            pct_display = pct_classical
            confidence = "low (n barely clears floor)" if n < LOW_CONFIDENCE_CELL_THRESHOLD else "standard"

        rows.append({
            "sample": sample,
            "n_ductal2": n,
            "n_classical": n_classical,
            "n_basal": n - n_classical,
            "pct_classical": pct_display,
            "majority_call": majority_call,
            "confidence": confidence,
        })

    return pd.DataFrame(rows).sort_values("n_ductal2", ascending=False).reset_index(drop=True)


def main() -> None:
    ductal2 = load_ductal2_expression()
    per_cell, gene_resolution = score_and_call(ductal2)
    report_score_diff_distribution(per_cell)

    PER_CELL_OUT.parent.mkdir(parents=True, exist_ok=True)
    per_cell.to_csv(PER_CELL_OUT)
    print(f"\nSaved per-cell scores to {PER_CELL_OUT}")

    with open(GENE_RESOLUTION_OUT, "w") as f:
        json.dump(gene_resolution, f, indent=2)
    print(f"Saved gene matched/remapped/dropped provenance to {GENE_RESOLUTION_OUT}")

    rollup = per_patient_rollup(per_cell)
    rollup.to_csv(PER_PATIENT_OUT, index=False)
    print(f"Saved per-patient rollup to {PER_PATIENT_OUT}")

    print(f"\nPer-patient subtype call (n>={MIN_CELLS_FOR_CALL} floor, "
          f"n>={LOW_CONFIDENCE_CELL_THRESHOLD} for standard confidence):")
    print(rollup.to_string(index=False))


if __name__ == "__main__":
    main()
