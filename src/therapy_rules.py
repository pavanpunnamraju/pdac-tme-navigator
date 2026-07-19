"""Step 8: therapy rules engine.

Spec: docs/step_8_therapy_rules_spec.md (logic-agent, full-text citation-verified
per spec section 5). Combines Step 6's per-sample TME composition
(data/processed/composition.csv) and Step 7's per-sample subtype calls
(data/processed/subtype_scores.csv) into a per-sample set of cited therapy
considerations.

Two independent axes (spec section 4, intro): subtype associations keyed on
`majority_call` alone, and TME-pattern associations keyed on the dominant-
pattern call from the section 2.3 two-proportion test. They are not combined
into a joint (subtype x pattern) lookup -- the underlying literatures are
largely independent bodies of work, and collapsing them would imply a joint
association the literature hasn't established.
"""

import csv
import json
import math
from pathlib import Path

from compose import LOW_COUNT_FLAG_THRESHOLD

REPO_ROOT = Path(__file__).resolve().parent.parent
COMPOSITION_CSV = REPO_ROOT / "data" / "processed" / "composition.csv"
SUBTYPE_CSV = REPO_ROOT / "data" / "processed" / "subtype_scores.csv"
OUT_CSV = REPO_ROOT / "data" / "processed" / "therapy_notes.csv"
OUT_JSON = REPO_ROOT / "data" / "processed" / "therapy_notes.json"
OUT_REPORT = REPO_ROOT / "data" / "processed" / "therapy_notes_report.md"

# Composition.csv's column names for the three dominance-contest axes
# (section 2.1 -- fibroblast/macrophage/T-cell only; stellate and the
# tumor/bystander compartments are not TME-pattern axes by design-doc scope).
AXIS_COLUMNS = {
    "fibroblast": "count_Fibroblast cell",
    "macrophage": "count_Macrophage cell",
    "tcell": "count_T cell",
}
# Title-case labels for prose (basis sentences); PATTERN_LABELS is the
# lowercase form spec section 2.3 uses in its "X-dominant" bucket names.
AXIS_LABELS = {"fibroblast": "Fibroblast", "macrophage": "Macrophage", "tcell": "T-cell"}
PATTERN_LABELS = {"fibroblast": "fibroblast", "macrophage": "macrophage", "tcell": "T-cell"}

Z_95 = 1.96


# ---------------------------------------------------------------------------
# Section 2.2-2.3: dominance test
# ---------------------------------------------------------------------------

def two_proportion_dominance(n1: int, n2: int) -> tuple[str, float, float]:
    """p, half-width, and 'dominant'/'mixed' for the top-2-count axes, per
    spec section 2.3: p = n1/(n1+n2), SE = sqrt(0.25/n), half-width = 1.96*SE,
    dominant iff |p-0.5| exceeds the half-width, else mixed.
    """
    n = n1 + n2
    p = n1 / n
    half_width = Z_95 * math.sqrt(0.25 / n)
    call = "dominant" if abs(p - 0.5) > half_width else "mixed"
    return call, p, half_width


def compute_dominance(counts: dict[str, int], threshold: int = LOW_COUNT_FLAG_THRESHOLD) -> dict:
    """Returns a dict describing the TME-pattern call for one sample:
    - eligible: {axis: n} for axes clearing the reliability floor (section 2.2)
    - ineligible: {axis: n} for axes below it (still reported, never dropped)
    - pattern: "T-cell dominant" / "macrophage dominant" / "mixed (macrophage/T-cell)"
      / a degenerate single-/zero-axis label if fewer than 2 axes are eligible
      (not expected in the current 8-sample cohort, but handled rather than
      assumed away)
    - basis: human-readable string with the p/half-width numbers, for the
      tme_pattern_basis output column
    - three_way_ambiguous: bool, per section 2.3's third-axis check
    """
    eligible = {axis: n for axis, n in counts.items() if n >= threshold}
    ineligible = {axis: n for axis, n in counts.items() if n < threshold}
    ranked = sorted(eligible.items(), key=lambda kv: kv[1], reverse=True)

    result = {
        "eligible": eligible,
        "ineligible": ineligible,
        "three_way_ambiguous": False,
    }

    if len(ranked) < 2:
        if len(ranked) == 1:
            axis, n = ranked[0]
            result["pattern"] = f"{AXIS_LABELS[axis]} only eligible axis (n={n}) -- no contest possible"
            result["basis"] = "only one axis clears n>=%d; dominance test requires two" % threshold
        else:
            result["pattern"] = "not assessable -- no axis clears the reliability floor"
            result["basis"] = "no axis clears n>=%d" % threshold
        return result

    (axis1, n1), (axis2, n2) = ranked[0], ranked[1]
    call, p, half_width = two_proportion_dominance(n1, n2)
    basis = (f"{AXIS_LABELS[axis1]} vs {AXIS_LABELS[axis2]}: n1={n1}, n2={n2}, "
             f"p={p:.3f}, 95% CI half-width={half_width * 100:.1f}pp")

    if call == "dominant":
        result["pattern"] = f"{PATTERN_LABELS[axis1]} dominant"
    else:
        result["pattern"] = "mixed (macrophage/T-cell)" if {axis1, axis2} == {"macrophage", "tcell"} \
            else f"mixed ({AXIS_LABELS[axis1]}/{AXIS_LABELS[axis2]})"
    result["basis"] = basis
    result["dominant_axis"] = axis1 if call == "dominant" else None
    result["axes_in_contest"] = frozenset({axis1, axis2})

    # Section 2.3: a third eligible axis is still fully reported, but only
    # enters the dominant-vs-mixed decision if it out-counts the current
    # top-2 (i.e. it would already be axis1/axis2 above, since `ranked` is
    # sorted by count). Separately, flag "3-way ambiguous" if the smallest
    # eligible axis is within the same CI band as the top axis.
    if len(ranked) >= 3:
        axis3, n3 = ranked[-1]
        _, p13, half_width13 = two_proportion_dominance(n1, n3) if n1 >= n3 else two_proportion_dominance(n3, n1)
        if abs(p13 - 0.5) <= half_width13:
            result["three_way_ambiguous"] = True
            result["basis"] += (f"; 3-way ambiguous: smallest eligible axis "
                                 f"({AXIS_LABELS[axis3]}, n={n3}) within CI band of top axis "
                                 f"({AXIS_LABELS[axis1]}, n={n1}) -- flag for manual review")

    return result


# ---------------------------------------------------------------------------
# Section 4: lookup table, as data
# ---------------------------------------------------------------------------

_MOFFITT_PRIMARY_CITE = "Moffitt RA et al., Nat Genet 2015, PMID 26343385 (subtype definition only)"
_PURIST_CITE = "Rashid NU et al., \"PurIST\", Clin Cancer Res 2020, PMC6942634 (Linehan/COMPASS ORR/DCR figures)"
_RETRO_931_CITE = "Real-world validation cohort, JCO Precision Oncology 2025, DOI 10.1200/PO-25-00197 / PMC12419025 (n=931 retrospective cohort)"
_PASS01_CITE = "Knox et al., \"PASS-01\", J Clin Oncol 2025, PMID 40929627 / JCO 2024 LBA4004 (randomized, subtype-stratified trial)"

_PROGNOSIS_NOTE = {
    "text": (
        "Basal-like tumors show worse response and survival on FOLFIRINOX specifically "
        "than classical tumors do, consistent across three independent data sources: "
        "Linehan phase 1b (basal ORR 0%, DCR 33% vs. classical ORR 40%, DCR 100%); "
        "COMPASS (basal ORR 10%, DCR 50% vs. classical ORR 36.7%, DCR 100%); n=931 "
        "real-world retrospective, FFX-treated (basal median OS 7.0mo vs. classical "
        "11.8mo, HR 1.86, 95% CI 1.49-2.33, P<.001)."
    ),
    "citation": f"{_MOFFITT_PRIMARY_CITE}; {_PURIST_CITE}; {_RETRO_931_CITE}.",
}

_CLASSICAL_REGIMEN_NOTE = {
    "text": (
        "Regimen choice for classical tumors is contested, not settled: the n=931 "
        "retrospective cohort found FOLFIRINOX conferred a substantial survival "
        "benefit over gemcitabine+nab-paclitaxel (GnP) among classical patients with "
        "good performance status, but PASS-01 -- the only randomized, subtype-"
        "stratified trial -- found the opposite for classical patients (GnP OS 13.9mo "
        "vs. mFOLFIRINOX 9.7mo, P=.047; ITT population also favored GnP, HR 1.57, 95% "
        "CI 1.08-2.28, P=.017). No regimen preference is asserted for classical "
        "patients."
    ),
    "citation": f"{_RETRO_931_CITE}; {_PASS01_CITE}.",
}

_BASAL_REGIMEN_NOTE = {
    "text": (
        "No regimen has been shown to significantly outperform another for "
        "basal-like disease. PASS-01: mFOLFIRINOX 7.5mo vs. GnP 8.9mo, P=.75 -- not "
        "significant; numerically favors GnP but underpowered to conclude a "
        "preference. No alternative-regimen recommendation is made for basal-like "
        "disease."
    ),
    "citation": _PASS01_CITE + ".",
}

# Section 4.1: keyed on majority_call. "classical" gets prognosis + the
# contested-regimen note (no recommendation); "basal" gets prognosis + the
# no-significant-alternative note. Neither call gets a regimen recommendation.
SUBTYPE_NOTES = {
    "classical": [_PROGNOSIS_NOTE, _CLASSICAL_REGIMEN_NOTE],
    "basal": [_PROGNOSIS_NOTE, _BASAL_REGIMEN_NOTE],
}

# MK447's exact confidence-qualifier prefix (spec section 3), applied inline
# to every subtype-conditioned sentence, not as a separate footnote.
MK447_CONFIDENCE_PREFIX = (
    "Classical call, low confidence (n=67, 95% CI ±12pp) — treat as directional: "
)

# MK336's exact placeholder text (spec section 3), verbatim.
MK336_PLACEHOLDER_TEXT = (
    "No subtype-conditioned therapy note: malignant subtype unavailable "
    "(insufficient ductal cells, n=13). TME-pattern-only considerations below."
)
MK336_PLACEHOLDER_NOTE = {
    "text": MK336_PLACEHOLDER_TEXT,
    "citation": "N/A -- data-availability note, not a literature association (Step 8 spec section 3).",
}

# Section 4.2: T-cell-dominant.
TCELL_NOTES = [
    {
        "text": (
            "Whole-sample T-cell dominance by dissociated-cell composition is not "
            "evidence of exploitable checkpoint-inhibitor sensitivity. PDAC "
            "checkpoint-inhibitor monotherapy has failed even in patients with "
            "measurable immune infiltration (zero RECIST responders among 27 "
            "patients on single-agent anti-CTLA-4; trial authors concluded the "
            "regimen is \"ineffective for the treatment of advanced pancreas "
            "cancer\"). This pipeline's dissociated, non-spatial composition data "
            "cannot distinguish true intratumoral T-cell infiltration from T cells "
            "sequestered in the stroma and physically excluded from tumor nests -- "
            "a PDAC-specific immune-exclusion mechanism (FAP+ stromal cells/matrix "
            "\"minimize T cell recruitment and mediate immune exclusion\")."
        ),
        "citation": (
            "Royal RE et al., J Immunother 2010;33(8):828-833, PMC7322622; "
            "\"Desmoplastic stroma restricts T cell extravasation and mediates "
            "immune exclusion and immunosuppression in solid tumors\" (PDAC-focused), "
            "Nat Commun 2023, PMC10120701."
        ),
    },
    {
        "text": (
            "Higher T-cell content, where it reflects genuine intratumoral "
            "infiltration (not assessable from this pipeline's non-spatial data "
            "alone), is associated with improved PDAC survival. The association "
            "holds specifically for T cells spatially proximal to cancer cells and "
            "holds regardless of desmoplasia/collagen density -- i.e. it is not "
            "simply a proxy for \"less stroma.\" This pipeline's composition data can "
            "only gesture at the association; it cannot measure the spatial "
            "proximity that has actually been shown to matter."
        ),
        "citation": "Carstens JL et al., Nat Commun 2017;8:15095, PMID 28447602.",
    },
]

# Section 4.3: macrophage-dominant.
MACROPHAGE_NOTES = [
    {
        "text": (
            "Macrophage/myeloid-enriched TME is associated with worse PDAC "
            "outcomes. In the APACT trial cohort, myeloid-enriched patients had "
            "shorter overall survival than adaptive (T/B-cell)-enriched patients. "
            "Separately, \"M0\" (unpolarized/undifferentiated) macrophages -- the "
            "closest match to this pipeline's undifferentiated \"Macrophage cell\" "
            "label -- were associated with worse OS across a full cohort (HR 1.23, "
            "P=1.6e-9), an association stronger in basal-like tumors specifically "
            "(HR 1.36, P=1.7e-4) than in classical tumors (HR 1.14, P=3.6e-3). This "
            "subtype-interaction detail is not acted on here: this pipeline reports "
            "the subtype and TME-pattern axes independently rather than modeling "
            "their interaction, so a macrophage-dominant + basal-like combination "
            "is not asserted to carry a combined risk signal beyond what each axis "
            "states on its own."
        ),
        "citation": (
            "\"Distinct immune cell infiltration patterns in pancreatic ductal "
            "adenocarcinoma (PDAC) exhibit divergent immune cell selection and "
            "immunosuppressive mechanisms\", Nat Commun 2024, "
            "DOI 10.1038/s41467-024-55424-2; \"Relevance of Immune Infiltration and "
            "Clinical Outcomes in Pancreatic Ductal Adenocarcinoma Subtypes\", "
            "PMC7815939."
        ),
    },
    {
        "text": (
            "Not claimed: this pipeline cannot resolve SPP1+ macrophage subsets, "
            "so the specific high-risk SPP1+-macrophage/myCAF co-occurrence "
            "ecotype reported elsewhere in the literature is not invoked here."
        ),
        "citation": (
            "Hwang WL et al., Nat Genet 2022;54:1178-1191 -- cited explicitly as "
            "not applicable to this pipeline's output; flagged per spec section 1, "
            "not invoked."
        ),
    },
]

# Section 4.4: mixed pattern -- no citation of its own, just a framing prefix
# that both the macrophage and T-cell notes get reported together under.
MIXED_PREFIX_NOTE = {
    "text": (
        "Neither macrophage nor T-cell compartment statistically dominates in "
        "this sample (per the section 2.3 two-proportion CI test). Report the "
        "macrophage-dominant and T-cell-dominant considerations below together: "
        "neither a purely myeloid-suppression framing nor a purely "
        "T-cell-engagement framing should be presented as this sample's defining "
        "TME feature."
    ),
    "citation": "No independent citation -- framing note only (Step 8 spec section 4.4).",
}

# Section 4.5: fibroblast/stromal caveat, attached only when fibroblast n
# clears LOW_COUNT_FLAG_THRESHOLD (BW21 only in this cohort).
FIBROBLAST_CAVEAT_NOTES = [
    {
        "text": (
            "Dense desmoplastic stroma is a documented physical/biochemical "
            "barrier to T-cell infiltration in PDAC specifically (FAP+ stromal "
            "cells/matrix), not simply neutral background."
        ),
        "citation": (
            "\"Desmoplastic stroma restricts T cell extravasation and mediates "
            "immune exclusion and immunosuppression in solid tumors\" (PDAC-focused), "
            "Nat Commun 2023, PMC10120701."
        ),
    },
    {
        "text": (
            "\"High stroma -> deplete the stroma\" is not a safe inference: in "
            "genetically-engineered mouse models (KPC/PKT strains), CAF/stromal "
            "depletion (Hedgehog pathway inhibition) produced more aggressive, "
            "less-differentiated tumors and shortened survival, with a "
            "human-tumor correlation (lower myofibroblast content in resected "
            "human PDAC also correlated with worse survival) -- mouse-model "
            "evidence with human correlational support, not a human interventional "
            "trial result. The worsened-outcome phenotype could be rescued by "
            "anti-CTLA-4 checkpoint blockade (~60% survival extension) but not by "
            "gemcitabine alone -- an immunosuppression-specific story, not a "
            "generic \"stroma is protective\" claim; it pairs with, rather than "
            "contradicts, the checkpoint-inhibitor-futility note above (which "
            "concerns un-depleted, naturally T-cell-excluded tumors). An "
            "independent Shh-deletion model reached the same qualitative result "
            "and found stromal depletion plus a Smoothened inhibitor accelerated "
            "disease versus depletion alone."
        ),
        "citation": (
            "Ozdemir BC et al., Cancer Cell 2014;25(6):719-734, PMC4180632; "
            "Rhim AD et al., Cancer Cell 2014;25(6):735-747, PMC4821630."
        ),
    },
    {
        "text": (
            "Bulk \"Fibroblast cell\" count alone cannot say whether the stroma "
            "is functionally protective (iCAF-dominant) or tumor-promoting "
            "(myCAF-dominant). The primary paper characterizing myCAF vs. iCAF "
            "does not itself make a prognosis/outcome claim for either subtype -- "
            "it treats their differential prognostic impact as an open hypothesis "
            "for future work, not an established result. Separately (and this "
            "part is an empirical finding, not a hypothesis): antigen-presenting "
            "CAFs (apCAFs, MHC-II+, lacking classical co-stimulatory molecules "
            "CD40/CD80/CD86) directly induce naive CD4+ T cells into Tregs in an "
            "antigen-specific manner."
        ),
        "citation": (
            "Ohlund D et al., J Exp Med 2017, DOI 10.1084/jem.20162024; Elyada E "
            "et al., Cancer Discov 2019;9(8):1102-1123 (apCAF/Treg finding); "
            "corroborated by Kerdidani/Huang et al., Cancer Cell 2022, PMC9197998."
        ),
    },
]


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def subtype_notes_for_sample(majority_call: str, confidence: str) -> list[dict]:
    if majority_call not in SUBTYPE_NOTES:
        return [MK336_PLACEHOLDER_NOTE]

    notes = [dict(n) for n in SUBTYPE_NOTES[majority_call]]
    if confidence.startswith("low"):
        for note in notes:
            note["text"] = MK447_CONFIDENCE_PREFIX + note["text"]
    return notes


def tme_pattern_notes_for_sample(pattern: str, axes_in_contest: frozenset | None) -> list[dict]:
    """Section 4.2-4.4 only define notes for macrophage-dominant, T-cell-
    dominant, and macrophage/T-cell mixed. `axes_in_contest` (the two axes
    the section 2.3 pairwise test actually ran on) gates the mixed case
    explicitly, rather than assuming mac/T-cell by default -- a mixed call
    between any other axis pair (possible in a future cohort where
    fibroblast reaches the top 2) has no section 4 entry and returns no
    notes rather than silently reusing the mac/T-cell ones.
    """
    if pattern == "macrophage dominant":
        return [dict(n) for n in MACROPHAGE_NOTES]
    if pattern == "T-cell dominant":
        return [dict(n) for n in TCELL_NOTES]
    if pattern.startswith("mixed"):
        if axes_in_contest == frozenset({"macrophage", "tcell"}):
            return [dict(MIXED_PREFIX_NOTE)] + [dict(n) for n in MACROPHAGE_NOTES] + [dict(n) for n in TCELL_NOTES]
        # Mixed call between a pair involving fibroblast: no section 4 entry
        # defined for this combination (not expected in the current cohort,
        # but not assumed away for future ones either).
        return []
    # Degenerate cases (fewer than 2 eligible axes) don't have a lookup entry;
    # none occur in the current 8-sample cohort (section 2.3).
    return []


def fibroblast_status_and_notes(fib_n: int, fib_prop: float) -> tuple[str, list[dict]]:
    if fib_n >= LOW_COUNT_FLAG_THRESHOLD:
        status = f"n={fib_n} ({fib_prop:.1%}), eligible -- see stromal caveat notes"
        return status, [dict(n) for n in FIBROBLAST_CAVEAT_NOTES]
    status = f"n={fib_n} ({fib_prop:.1%}), insufficient cells to assess (n<{LOW_COUNT_FLAG_THRESHOLD})"
    return status, []


def build_therapy_notes(composition_row: dict, subtype_row: dict | None) -> dict:
    sample = composition_row["sample"]
    counts = {axis: int(float(composition_row[col])) for axis, col in AXIS_COLUMNS.items()}
    dominance = compute_dominance(counts)

    fib_n = counts["fibroblast"]
    fib_prop = float(composition_row["prop_Fibroblast cell"])
    fibroblast_status, fibroblast_notes = fibroblast_status_and_notes(fib_n, fib_prop)

    raw_majority_call = subtype_row["majority_call"] if subtype_row is not None else "insufficient cells, no call"
    raw_confidence = subtype_row["confidence"] if subtype_row is not None else "n/a"
    is_called = raw_majority_call in SUBTYPE_NOTES
    subtype_call = raw_majority_call if is_called else "not called"
    subtype_confidence = raw_confidence

    notes = []
    notes.extend(subtype_notes_for_sample(raw_majority_call, raw_confidence))
    notes.extend(tme_pattern_notes_for_sample(dominance["pattern"], dominance.get("axes_in_contest")))
    notes.extend(fibroblast_notes)

    return {
        "sample": sample,
        "subtype_call": subtype_call,
        "subtype_confidence": subtype_confidence,
        "tme_pattern": dominance["pattern"],
        "tme_pattern_basis": dominance["basis"],
        "fibroblast_status": fibroblast_status,
        "three_way_ambiguous": dominance["three_way_ambiguous"],
        "therapy_notes": notes,
    }


def load_composition() -> list[dict]:
    with open(COMPOSITION_CSV, newline="") as f:
        return list(csv.DictReader(f))


def load_subtype_scores() -> dict[str, dict]:
    with open(SUBTYPE_CSV, newline="") as f:
        return {row["sample"]: row for row in csv.DictReader(f)}


def write_csv(rows: list[dict]) -> None:
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sample", "subtype_call", "subtype_confidence", "tme_pattern",
                  "tme_pattern_basis", "fibroblast_status", "three_way_ambiguous",
                  "therapy_notes"]
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out_row = dict(row)
            out_row["therapy_notes"] = json.dumps(row["therapy_notes"])
            writer.writerow(out_row)


def write_report(rows: list[dict]) -> None:
    lines = ["# Step 8 Therapy Notes\n"]
    for row in rows:
        lines.append(f"## {row['sample']}\n")
        lines.append(f"- **Subtype call:** {row['subtype_call']} "
                      f"(confidence: {row['subtype_confidence']})")
        lines.append(f"- **TME pattern:** {row['tme_pattern']} -- {row['tme_pattern_basis']}")
        lines.append(f"- **Fibroblast status:** {row['fibroblast_status']}")
        if row["three_way_ambiguous"]:
            lines.append("- **3-way ambiguous:** flagged for manual review (section 2.3)")
        lines.append("\n**Therapy notes:**\n")
        for note in row["therapy_notes"]:
            lines.append(f"- {note['text']}\n  - *Citation:* {note['citation']}")
        lines.append("")
    OUT_REPORT.write_text("\n".join(lines))


def main() -> None:
    composition_rows = load_composition()
    subtype_by_sample = load_subtype_scores()

    results = []
    for row in composition_rows:
        subtype_row = subtype_by_sample.get(row["sample"])
        results.append(build_therapy_notes(row, subtype_row))

    write_csv(results)
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    write_report(results)

    print(f"Wrote {OUT_CSV}, {OUT_JSON}, {OUT_REPORT} ({len(results)} samples)")
    print("\nPer-sample summary:")
    for row in results:
        print(f"  {row['sample']:>6s}  subtype={row['subtype_call']:<10s} "
              f"({row['subtype_confidence']:<25s})  tme={row['tme_pattern']:<28s}  "
              f"fibroblast={row['fibroblast_status']}")


if __name__ == "__main__":
    main()
