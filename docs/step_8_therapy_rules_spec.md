# Step 8 — Therapy Rules Engine Spec

Logic-agent spec for codebase-agent to implement. Defines the (subtype call, TME
pattern) lookup table, the operational definition of "dominant TME pattern," the
missing/low-confidence-data handling for MK336/MK447, and every citation backing an
association. No pipeline code here.

## 0. Inputs

- `data/processed/subtype_scores.csv` — per-sample `majority_call`
  (classical/basal-like/no call) and `confidence` (standard/low/n-a), from Step 7.
- `data/processed/composition.csv` — per-sample cell counts and proportions across
  10 coarse cell types, from Step 6.

## 1. Scope reminder (per design doc, Step 8 row)

Coarse TME patterns only — bulk fibroblast/macrophage/T-cell composition. No
myCAF/iCAF/apCAF or SPP1+ macrophage resolution: the local classifier's taxonomy
(`src/classify.py`, Zenodo 6024273 reference) only outputs coarse "Fibroblast cell" /
"Macrophage cell" labels, and Step 6 already found fibroblast/stellate counts too
sparse for finer resolution in this cohort (max 87 cells/sample, `MK362`
stellate=10). Any therapy note that implicitly requires CAF- or macrophage-subtype
resolution (e.g. the Hwang et al. SPP1+/myCAF high-risk ecotype) must be explicitly
flagged as **not assessable from this data**, not silently invoked.

## 2. Operationalizing "dominant TME pattern"

### 2.1 Which compartments compete

The design doc names three candidate axes: fibroblast, macrophage, T-cell
proportion. Stellate cells are excluded from the dominance contest — the design doc
doesn't name them as a fourth axis, and Step 6 found stellate counts even sparser
than fibroblast counts across this cohort (0–87 cells/sample, mostly single digits).
Ductal/acinar/endocrine/endothelial/B-cell compartments are not TME-pattern axes by
design-doc scope (they're tumor or bystander populations, not the
immune/stromal-composition question Step 8 is answering).

### 2.2 Reliability floor — reuse Step 6's own threshold, don't invent a new one

`src/compose.py` already defines `LOW_COUNT_FLAG_THRESHOLD = 20` — below this count, a
cell type's proportion for a sample is "too noisy to trust at face value (small-N
denominator effect)." Reusing that exact threshold here (rather than picking a new
number for Step 8) keeps the pipeline internally consistent and avoids the appearance
of threshold-shopping.

**Per-sample compartment counts** (from `composition.csv`):

| sample | Fibroblast n | Macrophage n | T-cell n |
|---|---|---|---|
| AM67 | 3 | 163 | 134 |
| BW21 | 50 | 940 | 5063 |
| MK336 | 1 | 63 | 1220 |
| MK359 | 7 | 224 | 248 |
| MK362 | 18 | 1067 | 827 |
| MK364 | 10 | 535 | 1051 |
| MK371 | 19 | 510 | 490 |
| MK447 | 0 | 75 | 545 |

**Fibroblast count clears n≥20 in exactly one sample (BW21) out of eight.** This is
the same finding the design doc already flagged for stellate cells, just less severe
for fibroblast — worth stating plainly rather than re-deriving each time: in this
cohort, bulk fibroblast content is essentially never assessable as a *dominant*-axis
call. It still gets reported per-sample (count + explicit "insufficient cells" flag
when n<20), because a fibroblast/stromal caveat is scientifically relevant regardless
of whether it can win the dominance contest (§4.3) — but the primary dominant-pattern
classification in this cohort will, empirically, always resolve to a macrophage/T-cell
question. Report this as a stated limitation of the input data, not a hidden
implementation detail.

**Rule:** an axis is *eligible* for the dominance contest only if its count ≥ 20 for
that sample. An ineligible axis is reported with its raw count and proportion but
carries a `"insufficient cells to assess"` flag and is excluded from the dominant-call
computation — never silently dropped from the output, never silently treated as zero.

### 2.3 Deciding "dominant" vs. "mixed" among eligible axes

Don't pick an arbitrary percentage-point margin (e.g. "top axis wins if it's 10pp
ahead") — Step 7 already established the project's convention for this exact
judgment call: don't impose an unprincipled cutoff when a real statistical test is
available. Reuse Step 7's own method, extended pairwise:

- Take the two highest-count eligible axes (`n_1 ≥ n_2`). Treat this as a two-category
  proportion test: `p = n_1 / (n_1 + n_2)`, `n = n_1 + n_2`, worst-case
  `SE = sqrt(0.25 / n)`, 95% CI half-width `= 1.96 × SE`.
- **Dominant** if `|p − 0.5|` exceeds that half-width (the top axis's lead is not
  explainable by sampling noise at a 50/50 null).
- **Mixed / co-dominant** otherwise — report both compartments, don't force a winner.
- If a third axis is also eligible (only possible for BW21 in this cohort, where
  fibroblast n=50), it is still reported in full but only enters the "dominant vs.
  mixed" decision if it out-counts the current top-2; the pairwise test always runs on
  whichever two axes have the highest counts. Flag any case where the smallest
  eligible axis is within the same CI band as the top axis as a "3-way ambiguous"
  case for manual review — none occur in the current 8-sample cohort, but codebase-agent
  shouldn't assume that generalizes to future samples.

**Computed calls for the current 8 samples** (fibroblast excluded from all but BW21
per §2.2; math shown so codebase-agent can unit-test against it):

| sample | eligible axes | dominant call |
|---|---|---|
| AM67 | Mac=163, T=134 | **mixed** (p=0.549, 95% CI half-width=5.7pp — gap not significant) |
| BW21 | Fib=50, Mac=940, T=5063 | **T-cell dominant** (T vs Mac: p=0.843, half-width=1.3pp) |
| MK336 | Mac=63, T=1220 | **T-cell dominant** (p=0.951, half-width=2.7pp) |
| MK359 | Mac=224, T=248 | **mixed** (p=0.525, half-width=4.5pp) |
| MK362 | Mac=1067, T=827 | **macrophage dominant** (p=0.563, half-width=2.2pp) |
| MK364 | Mac=535, T=1051 | **T-cell dominant** (p=0.663, half-width=2.5pp) |
| MK371 | Mac=510, T=490 | **mixed** (p=0.51, half-width=3.1pp) |
| MK447 | Mac=75, T=545 | **T-cell dominant** (p=0.879, half-width=3.9pp) |

Note this gives three possible dominant-pattern buckets in practice —
**T-cell-dominant**, **macrophage-dominant**, **mixed (macrophage/T-cell)** — plus a
fibroblast/stromal caveat layer that's reported alongside any of the three (§4.3), not
a fourth bucket of its own, because it structurally can't win the contest in this
cohort (§2.2).

## 3. Handling MK336 (no subtype call) and MK447 (low-confidence call)

Per Step 7, `subtype_scores.csv` already carries this distinction — Step 8 should
propagate it, not re-decide it:

- **MK336** (`majority_call = "insufficient cells, no call"`, n=13 Ductal2 cells,
  below Step 7's n≥43 floor): the *TME-pattern* axis is still fully computable (it's
  driven by the 1,220-cell T-cell / 63-cell macrophage counts, not the 13-cell ductal
  compartment) — compute and report it normally (T-cell dominant, per §2.3). But the
  output row must **not** synthesize a subtype-conditioned therapy note. Any
  association keyed on `(subtype, TME pattern)` is unavailable for MK336 by
  construction — only the TME-pattern-only associations (§4.2, §4.3) apply. Output
  text: `"No subtype-conditioned therapy note: malignant subtype unavailable
  (insufficient ductal cells, n=13). TME-pattern-only considerations below."` This is
  the "don't silently drop" requirement — the row exists, it's just visibly missing
  one axis, with the reason stated.
- **MK447** (`majority_call = "classical"`, `confidence = "low (n barely clears
  floor)"`, 85% of n=67, 95% CI ≈ ±12pp per Step 7 §4): report the full
  subtype-conditioned association, but every sentence of it must carry the confidence
  qualifier inline — not as a separate footnote a reader can miss. Output text
  prefix: `"Classical call, low confidence (n=67, 95% CI ±12pp) — treat as
  directional:"` before the FOLFIRINOX-preference note. This is the "don't silently
  treat as full-confidence" requirement.
- All other six samples: `confidence = "standard"` in Step 7's output — report
  associations without a qualifier prefix.

## 4. The lookup table

Two independent axes feed the output, not a single combined key: (a) subtype
associations, keyed on `majority_call` alone; (b) TME-pattern associations, keyed on
the §2.3 dominant-pattern call. A sample's full therapy note is the union of whichever
rows apply to it, each individually cited. This is deliberate, not a simplification —
collapsing them into one 6-cell `(subtype × pattern)` grid would imply the literature
has established joint associations for each specific combination, which it hasn't;
the two literatures (subtype-vs-chemo-response, TME-composition-vs-immune-therapy)
are largely independent bodies of work that happen to both apply to the same patient.

### 4.1 Subtype axis (keyed on `majority_call`)

**Correction from the first draft of this spec:** the original §4.1 asserted a clean
"classical → prefer FOLFIRINOX" recommendation as if the regimen-choice question were
settled. Full-text verification found that's only true of the *prognostic* claim
(basal-like does worse than classical, consistently), not the *regimen-choice* claim
(which drug is better for which subtype) — the two highest-quality 2025 sources on
regimen choice directly disagree with each other on the classical arm. Splitting these
two claims apart, rather than collapsing them, is the fix.

**(a) Prognosis (well-supported, consistent across independent cohorts):** basal-like
tumors show worse response and survival on FOLFIRINOX specifically than classical
tumors do, across three independent data sources:
- Linehan phase 1b trial data: basal ORR 0%, DCR 33% vs. classical ORR 40%, DCR 100%
  (on FOLFIRINOX arms).
- COMPASS trial data: basal ORR 10%, DCR 50% vs. classical ORR 36.7%, DCR 100%.
- Real-world retrospective validation (n=931 advanced PDAC, first-line FFX or GnP):
  among FFX-treated patients, basal median OS 7.0mo vs. classical 11.8mo (HR 1.86,
  95% CI 1.49–2.33, P<.001).

Citations: Moffitt RA et al., *Nat Genet* 2015, PMID 26343385 (defines the subtypes;
general resected-cohort prognosis HR 1.89 for basal-like, not FOLFIRINOX-specific —
cited here only for the subtype definition, not the chemo-response numbers above);
Rashid NU et al., "Purity Independent Subtyping of Tumors (PurIST)...", *Clin Cancer
Res* 2020, PMC6942634 (source of the Linehan/COMPASS ORR/DCR figures above — verified
by full-text fetch, not abstract-only); Real-world validation cohort, *JCO Precision
Oncology* 2025, DOI 10.1200/PO-25-00197 / PMC12419025 (source of the n=931 HR — also
full-text verified).

**(b) Regimen choice for classical tumors (contested — report both, don't pick a
side):** the n=931 retrospective cohort above found FFX "conferred a substantial
survival benefit over GnP" among classical-subtype patients with good performance
status. But **PASS-01** — the one randomized, head-to-head trial that specifically
tested regimen selection by molecular subtype (Knox et al., *J Clin Oncol* 2025,
PMID 40929627; molecular-subtype correlative results, *JCO* 2024 LBA4004) — found the
**opposite** for classical patients: GnP significantly outperformed mFOLFIRINOX (OS
13.9mo vs. 9.7mo, P=.047). Overall trial ITT population also favored GnP (HR 1.57,
95% CI 1.08–2.28, P=.017). Randomized evidence generally outranks retrospective
cohort evidence for a causal treatment-choice question, but PASS-01 is a smaller
phase II trial and this is a single trial versus a single large retrospective cohort —
not a settled question either way. **Output text should present both findings with
both citations and should not assert a regimen preference for classical patients** —
this is the one place in the spec where "no invented associations" specifically means
not resolving a live conflict in the source literature on the engine's behalf.

**(c) Regimen choice for basal-like tumors:** no regimen has been shown to
significantly outperform another for basal-like disease. PASS-01: mFFX 7.5mo vs. GnP
8.9mo, P=.75 (not significant; numerically favors GnP but underpowered to conclude
that). Report basal-like's worse prognosis (4.1a) without asserting a specific
alternative-regimen recommendation — the original draft's "some evidence favors
gemcitabine + nab-paclitaxel" line overstated what PASS-01 actually shows (a
non-significant numerical trend, not "evidence favors").

| `majority_call` | Association | Citation |
|---|---|---|
| `classical` | Consistently better response/survival on FOLFIRINOX than basal-like (4.1a). Whether classical patients themselves do better on FFX or on GnP is contested between the largest retrospective cohort (favors FFX) and the only randomized subtype-stratified trial (favors GnP, P=.047) — report both, no recommendation (4.1b). | See 4.1a/4.1b citations above. |
| `basal-like` (called `basal` in this pipeline's output) | Worse response/survival on FOLFIRINOX specifically, consistent across three cohorts (4.1a). No regimen shown to significantly outperform another for basal-like disease — PASS-01 found a non-significant numerical trend favoring GnP (P=.75); don't present this as an established preference (4.1c). | See 4.1a/4.1c citations above. |
| no call (MK336) | Not applicable — see §3. | — |

### 4.2 TME-pattern axis, T-cell-dominant

| Association | Citation |
|---|---|
| Whole-sample T-cell dominance by dissociated-cell composition is **not** evidence of exploitable checkpoint-inhibitor sensitivity. PDAC checkpoint-inhibitor monotherapy has failed even in patients with measurable immune infiltration — full-text confirms zero RECIST responders among 27 patients on single-agent anti-CTLA-4, with the authors' own conclusion stating the regimen is "ineffective for the treatment of advanced pancreas cancer." This pipeline's composition data — dissociated, non-spatial cell counts — cannot distinguish true intratumoral T-cell infiltration from T cells sequestered in the stroma and physically excluded from tumor nests, which full-text confirms is a PDAC-specific (not generic solid-tumor) mechanism in the second citation below: FAP+ stromal cells and their matrix "minimize T cell recruitment and mediate immune exclusion," and depleting them allows T-cell penetration into tumor nests. | Royal RE et al., "Phase 2 trial of single agent ipilimumab (anti-CTLA-4) for locally advanced or metastatic pancreatic adenocarcinoma," *J Immunother* 2010;33(8):828-833 (PMC7322622) — full-text verified; "Desmoplastic stroma restricts T cell extravasation and mediates immune exclusion and immunosuppression in solid tumors," *Nat Commun* 2023, PMC10120701 — full-text verified to be PDAC-focused (with cross-validation in lung and KPC mouse models), not merely "solid tumors" generically as the title implies. |
| Higher T-cell content, where it does reflect genuine intratumoral infiltration (not assessable from this pipeline's data alone), is associated with improved survival in PDAC. Full-text check found this association holds specifically for T cells spatially proximal to cancer cells, and — notably — **holds regardless of desmoplasia/collagen density**, i.e. the prognostic signal isn't just a proxy for "less stroma." That's a meaningfully stronger and more specific claim than "T-cell content correlates with survival" alone, but it also underlines why this pipeline's non-spatial composition data can only gesture at the association, not measure the thing that's actually been shown to matter (proximity, not raw count). | Carstens JL et al., "Spatial computation of intratumoral T cells correlates with survival of patients with pancreatic cancer," *Nat Commun* 2017;8:15095 (PMID 28447602) — full-text verified via search-indexed results (direct fetch blocked by a login wall; corroborated across three independent secondary summaries of the same paper, including the original PubMed listing). |

### 4.3 TME-pattern axis, macrophage-dominant

| Association | Citation |
|---|---|
| Macrophage/myeloid-enriched TME is associated with worse outcomes in PDAC — full-text confirms this with two independent, more specific findings than the original draft's generic framing. (i) Comparing myeloid-enriched vs. adaptive(T/B-cell)-enriched patients in the APACT trial cohort, myeloid enrichment had shorter overall survival. (ii) A separate cohort found "M0 macrophages" (unpolarized/undifferentiated, the closest match to this pipeline's undifferentiated "Macrophage cell" label) associated with worse OS across the full cohort (HR 1.23, P=1.6×10⁻⁹) — and, full-text confirms, this association is *stronger in basal-like tumors specifically* (HR 1.36, P=1.7×10⁻⁴) than in classical tumors (HR 1.14, P=3.6×10⁻³). That subtype-interaction detail is a real, citable finding this pipeline's independent-axes design (§4, intro) deliberately doesn't act on — noted here as a limitation: a macrophage-dominant + basal-like combination may carry a stronger risk signal than either axis alone suggests, but this pipeline reports the axes separately rather than modeling their interaction, and codebase-agent should not synthesize a combined-risk claim beyond what's stated here. | (i) "Distinct immune cell infiltration patterns in pancreatic ductal adenocarcinoma (PDAC) exhibit divergent immune cell selection and immunosuppressive mechanisms," *Nat Commun* 2024, DOI 10.1038/s41467-024-55424-2 — corroborated via three independent secondary summaries (direct fetch blocked by login wall). (ii) "Relevance of Immune Infiltration and Clinical Outcomes in Pancreatic Ductal Adenocarcinoma Subtypes," PMC7815939 — full-text verified. |
| **Not claimed:** this pipeline cannot resolve SPP1+ macrophage subsets, so the specific high-risk SPP1+-macrophage/myCAF co-occurrence ecotype cannot be invoked here even though it's a real, citable pattern in the literature. | Hwang WL et al., *Nat Genet* 2022;54:1178-1191 — cited explicitly as **not applicable** to this pipeline's output; flag, don't invoke, per §1. |

### 4.4 Mixed (macrophage/T-cell) pattern

No independent citation of its own — report §4.2 and §4.3's notes together with a
prefix stating neither compartment statistically dominates (from §2.3's CI test), so
neither a purely myeloid-suppression framing nor a purely T-cell-engagement framing
should be presented as the sample's defining TME feature.

### 4.5 Fibroblast/stromal caveat layer (attached whenever fibroblast n≥20, i.e. BW21 only in this cohort; all others get the "insufficient cells" flag from §2.2)

| Association | Citation |
|---|---|
| Dense desmoplastic stroma is a documented physical/biochemical barrier to T-cell infiltration — high stromal content is not simply neutral background. Full-text confirms this citation is specifically a PDAC finding (FAP+ stromal cells/matrix), not a generic cross-tumor-type inference. | Same *Nat Commun* 2023 desmoplasia citation as §4.2 (PMC10120701). |
| "High stroma → deplete the stroma" is **not** a safe inference: CAF/stromal depletion strategies (Hedgehog pathway inhibition) produced *more* aggressive, less differentiated tumors and *shortened* survival. Full-text adds two details the original draft's summary lost: (i) both findings are **genetically-engineered mouse models** (KPC/PKT strains) with a **human-tumor correlation** (lower myofibroblast content in resected human PDAC also correlated with worse survival) rather than a human interventional trial result — state it as mouse-model evidence with human correlational support, not as a human trial finding; (ii) Özdemir et al. also found the worsened-outcome phenotype could be *rescued* by anti-CTLA-4 checkpoint blockade (~60% survival extension) but *not* by gemcitabine alone — i.e., the depletion-is-bad finding is specifically an immunosuppression story, not a generic "stroma is protective" claim, and pairs with (rather than contradicts) the checkpoint-inhibitor-futility note in §4.2, since that note is about *un*-depleted, naturally T-cell-excluded tumors. Rhim et al.'s independent Shh-deletion model reached the same qualitative result (more undifferentiated, more metastatic, faster death) and additionally found the combination of stromal depletion + a Smoothened inhibitor (IPI-926) *accelerated* disease versus depletion alone — reinforcing, not just replicating, the same conclusion. | Özdemir BC et al., "Depletion of Carcinoma-Associated Fibroblasts and Fibrosis Induces Immunosuppression and Accelerates Pancreas Cancer with Reduced Survival," *Cancer Cell* 2014;25(6):719-734 (PMC4180632) — full-text verified; Rhim AD et al., "Stromal Elements Act to Restrain, Rather Than Support, Pancreatic Ductal Adenocarcinoma," *Cancer Cell* 2014;25(6):735-747 (PMC4821630) — full-text verified. |
| Bulk "Fibroblast cell" count alone cannot say whether the stroma is functionally protective (iCAF-dominant) or tumor-promoting (myCAF-dominant). **Correction from the first draft:** full-text check found Öhlund et al. 2017 does *not* itself make a prognosis/outcome claim for either subtype — it characterizes myCAF (αSMA-high, myofibroblastic, tumor-cell-contact-activated) vs. iCAF (IL-6-high, inflammatory, paracrine-activated) phenotypically and spatially, and only speculates, as a stated hypothesis for future work, that "certain CAF subtypes might have protumorigenic properties, whereas others might have antitumorigenic features" — not an empirical protective-vs-promoting finding. The spec's limitation claim should be "this pipeline cannot resolve which CAF subtype dominates, and even the primary literature defining these subtypes treats their differential prognostic impact as an open hypothesis, not an established result" — a *stronger* reason for the v1 scope limitation than originally stated, not a weaker one. Separately, full-text confirms apCAFs (MHC-II+, lacking classical co-stimulatory molecules CD40/CD80/CD86) directly induce naive CD4+ T cells into Tregs in an antigen-specific manner — this specific immunosuppressive mechanism *is* an empirical finding in Elyada et al. 2019, not speculative, and remains accurately citable on its own terms (independent of the myCAF/iCAF prognosis point above). | Öhlund D et al., "Distinct populations of inflammatory fibroblasts and myofibroblasts in pancreatic cancer," *J Exp Med* 2017, DOI 10.1084/jem.20162024 — full-text verified, claim scope corrected as above; Elyada E et al., "Cross-Species Single-Cell Analysis of Pancreatic Ductal Adenocarcinoma Reveals Antigen-Presenting Cancer-Associated Fibroblasts," *Cancer Discov* 2019;9(8):1102-1123 — Treg-induction claim confirmed via full-text-derived secondary sources; corroborated by Kerdidani/Huang et al., "Mesothelial cell-derived antigen-presenting cancer-associated fibroblasts induce expansion of regulatory T cells in pancreatic cancer," *Cancer Cell* 2022 (PMC9197998), which extends the same finding in vivo. |

## 5. Citation verification methodology note

**This section documents a full-text verification pass, superseding an earlier draft
that only did abstract/metadata-level verification (title, journal, year, and
PMID/PMC/DOI existence checks) — that pass was not sufficient and produced at least
one materially wrong claim (§4.1), caught only once the actual results sections were
read.**

Every citation in §4 was checked this round by fetching the full text (or, where a
publisher paywall or login wall blocked direct fetch — PurIST's real-world validation
full text initially, Elyada 2019 apCAF, Carstens 2017, s41467-024-55424-2 — by
cross-referencing multiple independent secondary sources that themselves quote or
closely paraphrase the paper's results/discussion sections, not just its abstract or
topic) against the specific finding attributed to it in the spec. Every citation in §4
now carries an inline note on how it was checked (full-text fetch vs. corroborated
secondary sources) so a future editor can see which citations still warrant a harder
verification pass if they become more load-bearing.

**What this pass changed:**

- **§4.1 (major correction):** the original draft's "classical → FOLFIRINOX preferred"
  recommendation conflated a well-supported prognostic claim (basal-like does worse on
  FOLFIRINOX, confirmed across three independent cohorts) with a regimen-choice claim
  that turned out to be actively contested in the literature. Pulling the actual
  PASS-01 trial results (a randomized, subtype-stratified trial the first draft hadn't
  found) showed GnP significantly outperforming FOLFIRINOX in classical patients
  (P=.047) — the opposite of what the retrospective cohort the first draft relied on
  showed. The spec no longer asserts a regimen preference for either subtype; it
  reports the conflict. This also means the "well-established, treat as settled"
  framing for this association in `docs/agents/logic-agent.md`'s domain-knowledge
  baseline is overstated for the regimen-choice half of the claim (the prognostic half
  still holds) — worth flagging to the human/orchestrator as a baseline-doc update,
  separate from this spec.
- **§4.3:** replaced a generic "TAM dominance associated with immunosuppression"
  framing with the two papers' actual, more specific quantitative findings, including
  a real basal-like-specific interaction effect (§4.3) that the spec explicitly
  declines to act on given its independent-axes design.
- **§4.5:** corrected an overclaim — Öhlund et al. 2017 does not itself establish
  myCAF/iCAF differential prognosis; that's stated in the primary paper as an open
  hypothesis, not a finding. The spec's limitation language is now stronger (less can
  be inferred from bulk fibroblast counts than the first draft implied) rather than
  weaker. The Özdemir/Rhim stromal-depletion citations held up under full-text check
  but needed the mouse-model/human-correlation distinction made explicit, since the
  original wording could be read as describing a human interventional trial.
- **§4.2:** held up under full-text check; tightened wording to reflect two details
  the first draft's abstract-level pass had missed (Carstens' "regardless of
  desmoplasia" finding; the desmoplasia-barrier paper being PDAC-specific rather than
  generically "solid tumors" as its title alone suggests).

Citations not yet fetched at full-text depth (blocked by login walls, corroborated
instead via ≥2 independent secondary sources quoting results-level text) should still
get a direct full-text read before this ships in a user-facing report, noted inline
per-citation above.

## 6. Output format

For each of the 8 samples, the rules engine should emit:

```
sample, subtype_call, subtype_confidence, tme_pattern, tme_pattern_basis,
fibroblast_status, therapy_notes[]
```

where `therapy_notes` is a list of `{text, citation}` pairs (never bare text without
an attached citation), assembled from whichever of §4.1–§4.5 apply per §3's
confidence-qualifier rules. A sample with a "no call" subtype (MK336) still gets a
row — just with `subtype_call = "not called"` and only TME-pattern-axis entries in
`therapy_notes`, per §3.

## Handoff

This spec (dominant-pattern definition with worked numbers, MK336/MK447 handling,
full lookup table with citations, output schema) is ready for codebase-agent to
implement as a plain-Python rules engine in `src/therapy_rules.py`. Logic-agent does
not write pipeline code — implementation, testing, and output file generation are
codebase-agent's responsibility once this spec is confirmed. Flag back to logic-agent
if the actual `composition.csv`/`subtype_scores.csv` at implementation time differ
from the values used in §2.3's worked table (e.g. if the pipeline is rerun and counts
shift) — the dominance *method* is fixed by this spec, but the specific per-sample
calls in the table above are a snapshot, not hardcoded truth.
