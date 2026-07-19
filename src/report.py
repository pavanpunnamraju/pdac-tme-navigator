"""Step 9: Per-patient report — aggregates Step 6 (composition), Step 7
(subtype scoring), and Step 8 (therapy notes) into one static HTML file.

Pure presentation layer: no new modeling, scoring, or interpretive judgment.
Therapy note text/citations are reproduced verbatim from therapy_notes.json
(Step 8 citation-verified them word-for-word) — only reformatted into HTML,
never reworded.
"""

import base64
import html
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
COMPOSITION_CSV = REPO_ROOT / "data" / "processed" / "composition.csv"
SUBTYPE_CSV = REPO_ROOT / "data" / "processed" / "subtype_scores.csv"
THERAPY_JSON = REPO_ROOT / "data" / "processed" / "therapy_notes.json"
OUT_HTML = REPO_ROOT / "reports" / "pdac_tme_report.html"

# Same order as src/compose.py's CELL_TYPES — kept identical across every
# patient's chart so a cell type sits at the same position on every bar,
# making cross-patient comparison possible at a glance.
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

CELL_TYPE_COLORS = {
    "Ductal cell type 1": "#9e9e9e",
    "Ductal cell type 2": "#8b1e3f",
    "Acinar cell": "#c9a227",
    "Endocrine cell": "#7fb069",
    "Fibroblast cell": "#4a7c59",
    "Stellate cell": "#2e8b8b",
    "Endothelial cell": "#3f6fa8",
    "Macrophage cell": "#e07a5f",
    "T cell": "#5c4d8a",
    "B cell": "#c76ba3",
}

# Cohort-level caveats surfaced during upstream steps that don't have an
# obvious per-field home in composition.csv / subtype_scores.csv /
# therapy_notes.json, but are load-bearing per project convention (see
# docs/agents/big-picture-agent.md Step 6/Step 4-5 entries). Reproduced here
# rather than dropped because assembling this report is the first place all
# three steps' outputs sit side by side.
COHORT_CAVEATS = [
    (
        "Ductal cell type 1 is likely undercounted (magnitude unresolved).",
        "Local samples show 0.5-4.4% Ductal1-of-ductal-cells vs. ~25.6% implied "
        "by the reference training ratio. A prevalence-conditioned recall check "
        "(35 reference patients, Pearson r=0.502) confirms the classifier's "
        "pooled 96.08% Ductal1 recall overstates performance at the local, "
        "much-rarer regime. Resolution reached: “partial calibration bias, "
        "magnitude unresolved” — not a hard systematic-undercount claim, "
        "the evidence doesn't support that strength at this n. A second, "
        "independent biological explanation also applies (FNA needles target "
        "tumor mass; malignant expansion effaces normal ductal tissue in PDAC), "
        "additive to, not a replacement for, the calibration-bias finding. "
        "Ductal1 is not used in Step 7 subtype scoring (Ductal2-only, per Peng "
        "et al. 2019), so this is a composition-chart caveat only.",
    ),
    (
        "BW21's B-cell and endothelial proportions: identity is solid, "
        "magnitude is not reliably estimated.",
        "BW21's elevated B-cell (9.3%, n=872) and endothelial (2.7%, n=253) "
        "counts were confirmed as real, sample-specific biology (clean marker "
        "programs, no QC/artifact red flags) compounded by Harmony's documented "
        "weak spot for cell types imbalanced across batches. Two BW21-dominated "
        "Leiden clusters were not force-merged via Harmony re-tuning. Any "
        "downstream quantitative claim about BW21's B-cell/endothelial "
        "abundance (or their ratio to other patients) should be treated as "
        "not reliably estimated from this pipeline, even though the cell "
        "identity calls themselves are trustworthy.",
    ),
    (
        "Cohort-wide Ductal2-vs-immune bimodality is FNA sampling-site "
        "heterogeneity, not a BW21-specific finding.",
        "4 samples are high-Ductal2/low-immune (AM67 84%, MK359 70%, MK371 "
        "69%, MK362 67%); 4 are low-Ductal2/immune-dominant (BW21, MK336, "
        "MK447, MK364). This is a cohort-level pattern attributed to FNA "
        "needle placement relative to tumor mass, not a property of any one "
        "sample.",
    ),
    (
        "Fibroblast/Stellate counts are low across all 8 samples (max 87 "
        "cells/sample) — CAF-subtype-resolved claims (myCAF/iCAF/apCAF) are "
        "not supported by this cohort at current resolution.",
        "Plausibly the same FNA limitation (aspirates favor cellular/epithelial "
        "content over stroma). Only BW21 (n=50) clears the n≥20 threshold "
        "used elsewhere in this pipeline for any fibroblast-conditioned "
        "statement; the other 7 samples' fibroblast/stellate figures are "
        "reported but should not be read as precise.",
    ),
]


def load_inputs():
    comp = pd.read_csv(COMPOSITION_CSV)
    # keep_default_na=False so the literal string "n/a" in the confidence
    # column (MK336's "not applicable") isn't parsed as NaN; na_values=[""]
    # keeps genuinely blank cells (e.g. MK336's pct_classical) as NaN.
    subtype = pd.read_csv(SUBTYPE_CSV, keep_default_na=False, na_values=[""])
    therapy = {d["sample"]: d for d in json.loads(THERAPY_JSON.read_text())}
    return comp, subtype, therapy


def make_composition_chart(row: pd.Series) -> str:
    """Horizontal bar chart of the 10 cell-type proportions for one sample.
    Returns a base64-encoded PNG data URI.
    """
    props = [row[f"prop_{ct}"] * 100 for ct in CELL_TYPES]
    counts = [int(row[f"count_{ct}"]) for ct in CELL_TYPES]
    colors = [CELL_TYPE_COLORS[ct] for ct in CELL_TYPES]

    fig, ax = plt.subplots(figsize=(6.5, 3.2), dpi=130)
    y_pos = range(len(CELL_TYPES))
    ax.barh(y_pos, props, color=colors)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(CELL_TYPES, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("% of total cells", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    for i, (p, c) in enumerate(zip(props, counts)):
        label = f"{p:.1f}%  (n={c})"
        ax.text(p + max(props) * 0.02, i, label, va="center", fontsize=7)
    ax.set_xlim(0, max(props) * 1.35 if max(props) > 0 else 1)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def esc(s) -> str:
    return html.escape(str(s)) if s is not None else ""


def subtype_block(srow: pd.Series | None, sample: str) -> str:
    if srow is None:
        return "<p><em>No subtype scoring row found for this sample.</em></p>"

    call = srow["majority_call"]
    conf = srow["confidence"]
    pct = srow["pct_classical"]
    n_ductal2 = int(srow["n_ductal2"])

    parts = [
        f'<p><span class="badge badge-{esc(call).replace(" ", "-")}">{esc(call)}</span> '
        f"&nbsp; confidence: <strong>{esc(conf)}</strong></p>"
    ]
    if pd.notna(pct):
        parts.append(
            f"<p>{pct:.1f}% classical / {100 - pct:.1f}% basal among "
            f"{n_ductal2} Ductal cell type 2 (malignant) cells.</p>"
        )
    else:
        parts.append(
            f"<p>No call: only {n_ductal2} Ductal cell type 2 cells "
            f"(below the confidence-floor threshold) &mdash; malignant subtype "
            f"cannot be reliably assessed for this sample.</p>"
        )

    if sample == "MK447":
        parts.append(
            '<p class="caveat"><strong>Low-confidence flag:</strong> n=67 '
            "clears the base insufficiency floor (n≥43 for ±15pp @ 95% CI) "
            "but not the stricter second tier (n≈97 for ±10pp @ 95% CI) "
            "used to flag lower-confidence calls. Reported as a call with an "
            "explicit lower-confidence flag, not withheld.</p>"
        )
    if sample == "MK336":
        parts.append(
            '<p class="caveat"><strong>No call:</strong> n=13 Ductal2 cells is '
            "below even the base insufficiency floor (n≥43). Same sample "
            "flagged T-cell-dominant in the Step 6 composition review. "
            "Therapy notes below are TME-pattern-only, with no "
            "subtype-conditioned regimen note.</p>"
        )
    return "\n".join(parts)


def therapy_block(tnotes: dict | None) -> str:
    if tnotes is None:
        return "<p><em>No therapy notes found for this sample.</em></p>"

    parts = [
        f'<p><strong>TME pattern:</strong> {esc(tnotes["tme_pattern"])} '
        f'&mdash; {esc(tnotes["tme_pattern_basis"])}</p>',
        f'<p><strong>Fibroblast status:</strong> {esc(tnotes["fibroblast_status"])}</p>',
    ]
    if tnotes.get("three_way_ambiguous"):
        parts.append(
            '<p class="caveat"><strong>Three-way ambiguous TME call.</strong></p>'
        )

    parts.append('<div class="therapy-notes">')
    for note in tnotes["therapy_notes"]:
        parts.append(
            "<div class='note'>"
            f"<p class='note-text'>{esc(note['text'])}</p>"
            f"<p class='note-citation'><em>Citation:</em> {esc(note['citation'])}</p>"
            "</div>"
        )
    parts.append("</div>")
    return "\n".join(parts)


PAGE_CSS = """
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
       max-width: 900px; margin: 0 auto; padding: 24px 20px 80px;
       color: #1a1a1a; line-height: 1.5; }
h1 { font-size: 1.6rem; margin-bottom: 4px; }
.subtitle { color: #555; margin-top: 0; margin-bottom: 24px; }
h2.patient-header { font-size: 1.3rem; border-top: 3px solid #333;
       padding-top: 18px; margin-top: 40px; }
h3 { font-size: 1.02rem; margin-bottom: 6px; color: #333; }
.toc { columns: 4; margin-bottom: 28px; }
.toc a { text-decoration: none; color: #1a4d7a; }
.chart-img { max-width: 100%; }
.badge { display: inline-block; padding: 2px 10px; border-radius: 12px;
       font-weight: 600; font-size: 0.85rem; color: #fff; }
.badge-classical { background: #2e6f40; }
.badge-basal { background: #a4373a; }
.badge-not-called { background: #777; }
.caveat { background: #fff6e0; border-left: 4px solid #d8a300;
       padding: 8px 12px; font-size: 0.92rem; }
.cohort-caveats { background: #f4f6f9; border: 1px solid #d7dde5;
       border-radius: 6px; padding: 14px 18px; margin-bottom: 30px; }
.cohort-caveats h4 { margin-bottom: 4px; }
.cohort-caveats p { margin-top: 2px; font-size: 0.92rem; color: #333; }
.therapy-notes .note { border-bottom: 1px solid #e5e5e5; padding: 10px 0; }
.note-text { margin: 0 0 4px; }
.note-citation { margin: 0; font-size: 0.85rem; color: #555; }
.grid { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }
@media (max-width: 700px) { .grid { grid-template-columns: 1fr; } .toc { columns: 2; } }
footer { margin-top: 50px; color: #888; font-size: 0.8rem; }
"""


def build_report(comp: pd.DataFrame, subtype: pd.DataFrame, therapy: dict) -> str:
    samples = sorted(comp["sample"].tolist())

    toc = "".join(f'<li><a href="#{s}">{s}</a></li>' for s in samples)

    cohort_html = ['<div class="cohort-caveats"><h3>Cohort-level caveats</h3>']
    for title, body in COHORT_CAVEATS:
        cohort_html.append(f"<h4>{esc(title)}</h4><p>{esc(body)}</p>")
    cohort_html.append("</div>")

    sections = []
    for s in samples:
        comp_row = comp[comp["sample"] == s].iloc[0]
        subtype_rows = subtype[subtype["sample"] == s]
        srow = subtype_rows.iloc[0] if len(subtype_rows) else None
        tnotes = therapy.get(s)

        chart_uri = make_composition_chart(comp_row)

        bw21_note = ""
        if s == "BW21":
            bw21_note = (
                '<p class="caveat">See cohort-level caveat above: B-cell and '
                "endothelial proportions for this sample carry a known "
                "Harmony batch-integration magnitude caveat.</p>"
            )

        sections.append(
            f"""
<h2 class="patient-header" id="{s}">{s}</h2>
<div class="grid">
  <div>
    <h3>TME composition ({int(comp_row['total_cells'])} cells)</h3>
    <img class="chart-img" src="{chart_uri}" alt="Composition chart for {s}">
    {bw21_note}
  </div>
  <div>
    <h3>Moffitt subtype call</h3>
    {subtype_block(srow, s)}
  </div>
</div>
<h3>Therapy considerations</h3>
{therapy_block(tnotes)}
"""
        )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>PDAC TME Navigator — Per-Patient Report</title>
<style>{PAGE_CSS}</style>
</head>
<body>
<h1>PDAC TME Navigator</h1>
<p class="subtitle">Per-patient composition, Moffitt subtype call, and
citation-verified therapy considerations &mdash; {len(samples)} samples.
Aggregation only: no new modeling or interpretive judgment made at this
step; all calls and citations carried forward verbatim from Steps 6-8.</p>

<ul class="toc">{toc}</ul>

{"".join(cohort_html)}

{"".join(sections)}

<footer>Generated by src/report.py (Step 9). Sources: composition.csv,
subtype_scores.csv, therapy_notes.json.</footer>
</body>
</html>
"""


def main():
    comp, subtype, therapy = load_inputs()
    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(build_report(comp, subtype, therapy))
    print(f"Wrote {OUT_HTML}")


if __name__ == "__main__":
    main()
