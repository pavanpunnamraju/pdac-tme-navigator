# Step 7 — Moffitt Classical/Basal-like Subtype Scoring Spec

Logic-agent spec for codebase-agent to implement. Defines gene lists, cell
selection, scoring method, and per-patient rollup logic only — no pipeline
code here.

## 1. Moffitt gene signatures — verified against primary source

Pulled from the `pdacR` package (Moffitt lab's own curation repo,
https://github.com/rmoffitt/pdacR), which parses these directly from the
original supplementary table and cites **Moffitt RA et al., *Nature
Genetics* 2015, PMID 26343385**. Extracted directly from the shipped
`gene_lists.rds` in that repo (not reconstructed from memory) and
cross-checked several distinctive genes (GPR87, VGLL1, ANXA8L2 as basal;
LYZ, TFF1/2 as classical) against independent secondary literature.

**Basal-like (25 genes):**
```
ANXA8L2, AREG, CST6, CTSL2, DHRS9, FAM83A, FGFBP1, GPR87, KRT15, KRT17,
KRT6A, KRT6C, KRT7, LEMD1, LY6D, S100A2, SCEL, SERPINB3, SERPINB4, SLC2A1,
SPRR1B, SPRR3, TNS4, UCA1, VGLL1
```

**Classical (25 genes):**
```
AGR2, AGR3, ANXA10, ATAD4, BTNL8, CDH17, CEACAM6, CLRN3, CTSE, CYP3A7,
FAM3D, KRT20, LGALS4, LOC400573, LYZ, MYO1A, PLA2G10, REG4, SPINK4,
ST6GALNAC1, TFF1, TFF2, TFF3, TSPAN8, VSIG2
```

**Flag for codebase-agent — symbol currency check required before scoring:**
Three symbols are 2015-vintage and may not match current gene symbols in
`adata.var_names`:
- `CTSL2` — HGNC has since split this from CTSL; likely needs mapping to a
  current symbol.
- `ANXA8L2` — naming inconsistent across GENCODE versions (may appear as
  ANXA8L1/ANXA8B depending on annotation).
- `LOC400573` — Entrez placeholder ID, not a real symbol; needs lookup
  against whatever reference GTF was used for alignment.

Before running `score_genes`, check the actual overlap between these 50
symbols and `adata.var_names`. If any silently drop to zero matches, that's
a scoring integrity issue — a 25-gene signature losing even 2–3 genes to a
name mismatch measurably changes the score. Report the matched/dropped gene
count per signature in the Step 7 output.

## 2. Which cells to score

**Ductal cell type 2 only** — not Ductal1 + Ductal2.

This is already correctly specified in `src/classify.py` (lines 207–211):
`predicted_cell_type == "Ductal cell type 2"` is the malignant/CNV-bearing
population per **Peng et al. 2019, *Cell Res* 29:725-738**, vs. Ductal1
being CNV-defined normal-like ductal tissue. Scoring Ductal1 against a
malignant-tumor-derived subtype signature would be a category error —
normal ductal cells don't have a classical/basal identity; that axis is
specific to transformed epithelium.

Side note carried over from the Step 6 composition review: this also gives
Ductal1's near-absence in local samples (0.5–4.4% of ductal cells) a third,
genuinely biological explanation beyond the previously-flagged
prevalence-dependent classifier recall degradation — FNA needles are aimed
at tumor mass, and normal ductal tissue is frequently effaced by malignant
expansion in PDAC. This adds to, but doesn't change, the "partial
calibration bias, magnitude unresolved" conclusion already reached for
Ductal1 — moot for Step 7 either way since Ductal1 isn't used downstream.

## 3. Scoring method and call logic

- Run `sc.tl.score_genes` twice on the **Step 3 normalized/log data** (not
  raw counts), restricted to Ductal2 cells: `classical_score` (25-gene
  classical list) and `basal_score` (25-gene basal list), default
  control-gene-set correction.
- **Per-cell call:** `classical` if `classical_score > basal_score`, else
  `basal`. Use the relative comparison, not an absolute threshold on either
  score alone — both scores share the same background-correction noise, so
  the difference cancels shared bias in a way a fixed cutoff on one score
  wouldn't. This is the standard approach in the Moffitt/PurIST lineage of
  work, not an invented rule.
- Do **not** impose an arbitrary margin/ambiguity threshold on
  `|classical_score - basal_score|` to declare cells "indeterminate" —
  there's no principled null distribution here to set that cutoff, and per
  this project's baseline domain knowledge, subtype is a spectrum, so a
  small margin isn't itself evidence of poor data quality. Report the score
  difference distribution in the output so a reviewer can see it, but don't
  hard-filter cells on it.
- **Per-patient rollup:** majority label among that patient's Ductal2 cells
  (`% classical` = classical-call cells / total Ductal2 cells for that
  patient). Report both the majority call **and** the percentage split —
  not just a binary label. A patient at 55/45 and a patient at 95/5 both get
  "classical" under majority-only rollup but represent very different
  biological pictures; collapsing that away overclaims precision.

## 4. Insufficient-cell threshold for per-patient calls

Per-sample Ductal2 (malignant) cell counts, from
`data/processed/local_cell_type_predictions.csv`:

| sample | Ductal2 n |
|---|---|
| MK362 | 4118 |
| MK371 | 2373 |
| BW21 | 2003 |
| AM67 | 1642 |
| MK359 | 1180 |
| MK364 | 526 |
| **MK447** | **67** |
| **MK336** | **13** |

Don't pick a round-number cutoff arbitrarily — derive it from the actual
statistic being reported. A per-patient `% classical` is a binomial
proportion; its sampling noise at the true 50/50 (worst case) point is
`SE = sqrt(0.25/n)`. Requiring a 95% CI half-width no wider than ±15
percentage points solves to `n ≥ 1.96² × 0.25 / 0.15² ≈ 43`.

- **MK336 (n=13):** fails by a wide margin — CI half-width ≈ ±27pp at
  50/50. This sample's whole ductal compartment is thin to begin with (13
  of 1,361 total cells; the same sample that was 89.6% T-cell-dominant in
  the Step 6 composition review). Report as "insufficient cells, no
  subtype call reported" rather than a number that looks precise but isn't.
- **MK447 (n=67):** clears the n≥43 floor but barely — CI half-width ≈
  ±12pp. Report the call but flag it as lower-confidence in the output
  rather than presenting it at the same confidence level as the other six
  samples.
- The other five samples (526–4,118 cells) are comfortably powered.

## Handoff

This spec (gene lists, cell selection, scoring method, rollup logic,
confidence floor) is ready for codebase-agent to implement in Step 7.
Logic-agent does not write pipeline code — implementation, testing, and
output file generation are codebase-agent's responsibility once this spec
is confirmed.
