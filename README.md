# PDAC TME Navigator

Given single-cell RNA-seq from a pancreatic ductal adenocarcinoma (PDAC) biopsy, this
pipeline answers three things: what's in the tumor microenvironment (TME), what
molecular subtype the tumor is, and what that combination suggests therapeutically —
with every non-obvious claim traceable to a citation. It is **not** a diagnostic tool
and makes **no outcome/survival predictions**: 8 patients is not enough to support that,
and it isn't attempted (see [design doc](docs/plans/2026-07-17-pdac-tme-navigator-design.md)).

**Deliverable**: [`reports/pdac_tme_report.html`](reports/pdac_tme_report.html) — a
per-patient report (TME composition chart, subtype call, cited therapy considerations)
for all 8 samples.

## Pipeline

| Step | What it does | Code |
|---|---|---|
| 1. Ingest | Load the 8 CellRanger-format samples into AnnData | [`src/ingest.py`](src/ingest.py) |
| 2. QC | Per-sample adaptive filtering (MAD-based mito cutoff, min genes/cell) + doublet detection | [`src/qc.py`](src/qc.py) |
| 3. Normalize | Total-count normalization, log1p, HVG selection, scaling | [`src/normalize.py`](src/normalize.py) |
| 4. Integrate | Batch-correct across the 8 patients (PCA → Harmony) into one joint embedding | [`src/integrate.py`](src/integrate.py) |
| 5. Cell-type classifier | LightGBM trained on the Zenodo 6024273 reference, applied to local cells; evaluated against a naive marker-threshold baseline on held-out reference data | [`src/classify.py`](src/classify.py) |
| 6. Composition | Aggregate predicted cell types into per-patient TME composition (counts + proportions, 10 cell types) | [`src/compose.py`](src/compose.py) |
| 7. Subtype scoring | Score malignant (Ductal2) cells against Moffitt classical/basal-like signatures, roll up to a per-patient call | [`src/subtype_score.py`](src/subtype_score.py) |
| 8. Therapy rules | Lookup (subtype, dominant TME pattern) → literature-backed, citation-verified therapy considerations | [`src/therapy_rules.py`](src/therapy_rules.py) |
| 9. Report | Combine steps 6-8 into a static per-patient HTML report | [`src/report.py`](src/report.py) |

## Results / metrics

- **Cell-type classifier (the project's one hard go/no-go gate): beats baseline.**
  - LightGBM: **0.9812 accuracy / 0.9713 macro-F1** (held-out split of the reference data)
  - Naive marker-threshold baseline: **0.4758 accuracy / 0.4472 macro-F1**
- 8 local samples, **26,804 integrated cells**, 10-category cell-type taxonomy (Ductal
  type 1/2, Acinar, Endocrine, Fibroblast, Stellate, Endothelial, Macrophage, T cell,
  B cell).
- Full pipeline runs end-to-end on all 8 samples, producing the per-patient report.

## How to reproduce

Run in order; each script reads the previous step's output from `data/processed/` and
writes its own output there (large intermediates — `.h5ad`, model files — are gitignored;
small per-sample/per-patient CSVs and JSON are committed).

```
python src/ingest.py          # FNA_scRNA _JJL/<SAMPLE>_filtered_feature_bc_matrix/ -> data/processed/ingested.h5ad
python src/qc.py              # ingested.h5ad -> qc_filtered.h5ad, qc_annotated.h5ad, docs/qc_report.md
python src/normalize.py       # qc_filtered.h5ad -> normalized.h5ad
python src/integrate.py       # normalized.h5ad -> integrated.h5ad
python src/reference_prep/export_pk_all.R   # Zenodo 6024273 pk_all.rds -> mtx export (requires R + Seurat)
python src/reference_prep/build_reference.py  # mtx export -> data/reference/pk_all.h5ad
python src/classify.py        # integrated.h5ad + pk_all.h5ad -> classified.h5ad, local_cell_type_predictions.csv
python src/compose.py         # local_cell_type_predictions.csv -> data/processed/composition.csv
python src/subtype_score.py   # normalized.h5ad + local_cell_type_predictions.csv -> subtype_scores.csv, subtype_scores_per_cell.csv, gene_resolution.json
python src/therapy_rules.py   # composition.csv + subtype_scores.csv -> therapy_notes.json, therapy_notes.csv, therapy_notes_report.md
python src/report.py          # composition.csv + subtype_scores.csv + therapy_notes.json -> reports/pdac_tme_report.html
```

Raw input (`FNA_scRNA _JJL/`, patient-derived) and the Zenodo reference are not part of
this repo and must be obtained separately; see the design doc's Data section.

## Reference dataset attribution

[Zenodo 6024273](https://zenodo.org/records/6024273) — "Establishment of a reference
single-cell RNA sequencing dataset for human pancreatic adenocarcinoma" (iScience, 2022).
Used **only** as labeled training data for the Step 5 classifier — it is not part of
this pipeline's own findings or claims about the local 8 patients.

## Known limitations

These are documented in detail in [`docs/qc_report.md`](docs/qc_report.md), the agent
docs under [`docs/agents/`](docs/agents/), and inline in the
[report](reports/pdac_tme_report.html) itself — this section points to them rather than
re-litigating them:

- **Doublet detection found no separable doublet population in any sample** (unimodal
  scrublet score distributions, no bimodal shoulder) — a known scrublet blind spot for
  homotypic doublets and low-diversity samples, not evidence the data is doublet-free.
  See [`docs/qc_report.md`](docs/qc_report.md).
- **BW21's B-cell and endothelial proportions**: cell identity is solid, but the
  magnitude carries a Harmony batch-imbalance caveat (Harmony has a documented weak spot
  for cell types imbalanced across integration batches). See the cohort-caveats section
  of the [report](reports/pdac_tme_report.html).
- **Ductal cell type 1 is likely undercounted** — local samples show 0.5-4.4%
  Ductal1-of-ductal-cells vs. ~25.6% implied by the reference's training ratio; a
  prevalence-conditioned recall check confirms partial calibration bias (magnitude
  unresolved), compounded by a plausible independent biological cause (FNA needles
  target tumor mass, effacing normal ductal tissue). Not used in Step 7 subtype scoring.
- **Fibroblast/Stellate counts are too sparse (max 87 cells/sample) for CAF-subtype
  resolution** (myCAF/iCAF/apCAF) — only bulk fibroblast composition is reported;
  CAF-subtype-conditioned therapy claims are out of scope for this cohort.
- **MK336 has no subtype call** (13 Ductal2 cells, below the confidence floor); **MK447's
  call is flagged low-confidence** (67 cells, clears the base floor but not the stricter
  tier). Both are reported plainly, not hidden.
- **Cohort-wide Ductal2-vs-immune bimodality** across the 8 samples is attributed to FNA
  sampling-site heterogeneity (needle placement relative to tumor mass), not a
  per-sample anomaly.

## Status

Steps 1-9 complete. Portfolio-ready.
