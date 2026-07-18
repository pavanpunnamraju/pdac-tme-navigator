# PDAC TME Navigator — Design Doc

Date: 2026-07-17
Status: Approved, entering build phase

## Problem

Given single-cell RNA-seq data from a patient's pancreatic tumor biopsy, answer three
things: what's in the tumor microenvironment, what molecular subtype is the tumor, and
what does that combination suggest therapeutically.

## Data

- Source: PaCMEN-downloaded PDAC FNA scRNA-seq data,
  `FNA_scRNA _JJL/` (local, not committed to git — patient-derived data).
- 9 samples total. 8 in standard 10x CellRanger format (`matrix.mtx`, `barcodes.tsv`,
  `features.tsv`): `MK362`, `AM67`, `MK364`, `MK336`, `BW21`, `MK371`, `MK359`, `MK447`.
- 1 sample (`RS01`) is ddSeq/Dropseq-style CSV output — different schema, **excluded from
  v1** to avoid a format-conversion detour; revisit as future work.
- No cell-type annotations or clinical metadata bundled — raw UMI counts only.

## Why cell-level, not patient-level

9 patients is not enough for any patient-level supervised model (outcome prediction,
subtype-from-patient-features, etc.) — that's a statistical non-starter regardless of
method. At the cell level, the 8 usable samples give tens of thousands of cells, enough
for a real trained classifier. The ML core of this project operates at the cell level;
patient-level outputs are aggregated from cell-level predictions.

## Reference dataset for training labels

[Zenodo 6024273](https://zenodo.org/records/6024273) — "Establishment of a reference
single-cell RNA sequencing dataset for human pancreatic adenocarcinoma" (iScience, 2022).
A pre-integrated, pre-annotated reference combining 5 public PDAC scRNA-seq cohorts
(Peng et al./GSE111672, GSE155698, GSE154778, GSM4293555), with cell-type annotations,
malignant-cell subtype calls, and CAF subtype calls already assigned. Used purely as the
labeled training set for the cell-type classifier — not part of the pipeline's own
identity or claims.

## Pipeline

| Step | What happens | Tooling |
|---|---|---|
| 1. Ingest | Load the 8 CellRanger-format samples | scanpy/anndata |
| 2. QC | Per-sample filtering: min genes/cell, max mito %, doublet detection | scanpy, scrublet |
| 3. Normalize | Total-count norm, log1p, HVG selection, scale | scanpy |
| 4. Integrate | Batch-correct across the 8 patients into one joint embedding | Harmony |
| 5. Cell-type classifier (ML core) | Train supervised classifier on Zenodo 6024273 reference labels, apply to local cells. Eval via held-out split of the reference (F1/accuracy); baseline = naive marker-threshold heuristic | scikit-learn/LightGBM |
| 6. Composition | Aggregate predicted cell types into per-patient TME composition breakdown | pandas |
| 7. Subtype scoring | Score malignant cells against Moffitt classical/basal-like gene signatures | `scanpy.tl.score_genes` |
| 8. Therapy rules | Lookup table: (subtype, dominant TME pattern) → literature-backed therapeutic considerations, cited | plain Python rules engine |
| 9. Report | Per-patient summary: composition chart + subtype call + therapy notes | notebook/script output, no UI in v1 |

## Success criteria

- Classifier beats the naive marker-threshold baseline on held-out reference data
  (F1/accuracy).
- Full pipeline runs end-to-end on the 8 local samples, producing a per-patient report.
- Every non-obvious claim (subtype signature, therapy association) is traceable to a
  citation — no invented associations.
- Portfolio-ready repo: README, reproducible pipeline, clear metrics.

## Explicitly out of scope for v1

- `RS01` sample (format mismatch).
- Any patient-level outcome/survival prediction (insufficient n).
- App UI / LLM layer (future phase).
- Bulk RNA-seq deconvolution or external clinical-outcome cohorts (TCGA-PAAD, ICGC, etc.)
  — not needed for this scope; revisit only if a future phase adds outcome prediction.

## Milestones (self-paced, ~1-2 weeks)

- **Wk1**: ingest → QC → integrate → classifier trained & evaluated (go/no-go gate: beats
  baseline?).
- **Wk1.5**: composition aggregation + subtype scoring.
- **Wk2**: therapy rules layer + per-patient report + README/portfolio packaging.

## Delegation

- **Codebase agent**: ingestion, QC, normalization, integration, classifier
  training/eval scripts.
- **Logic agent**: subtype scoring logic, therapy rules engine.
- **Git agent**: repo init, commits per milestone.
- **Lead/big-picture agent**: checks each milestone's metric gate before advancing.
- **UI agent**: idle this phase — no UI in the MVP.
