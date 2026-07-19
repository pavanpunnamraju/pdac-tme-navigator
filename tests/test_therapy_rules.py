"""Regression test for src/therapy_rules.py's dominance test against the
worked numbers in docs/step_8_therapy_rules_spec.md section 2.3 -- catches
drift if composition.csv is regenerated and the per-sample calls change
silently. The dominance *method* is fixed by the spec; this test pins the
specific 8-sample snapshot values the spec worked through by hand.

Run directly: .venv/bin/python tests/test_therapy_rules.py
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from therapy_rules import compute_dominance  # noqa: E402

# sample -> (counts, expected_pattern, expected_p, expected_half_width_pp)
# Counts and expected values transcribed from spec section 2.2 (per-sample
# compartment counts) and section 2.3 (worked dominance table).
CASES = {
    "AM67": ({"fibroblast": 3, "macrophage": 163, "tcell": 134},
             "mixed (macrophage/T-cell)", 0.549, 5.7),
    "BW21": ({"fibroblast": 50, "macrophage": 940, "tcell": 5063},
              "T-cell dominant", 0.843, 1.3),
    "MK336": ({"fibroblast": 1, "macrophage": 63, "tcell": 1220},
              "T-cell dominant", 0.951, 2.7),
    "MK359": ({"fibroblast": 7, "macrophage": 224, "tcell": 248},
              "mixed (macrophage/T-cell)", 0.525, 4.5),
    "MK362": ({"fibroblast": 18, "macrophage": 1067, "tcell": 827},
               "macrophage dominant", 0.563, 2.2),
    "MK364": ({"fibroblast": 10, "macrophage": 535, "tcell": 1051},
               "T-cell dominant", 0.663, 2.5),
    "MK371": ({"fibroblast": 19, "macrophage": 510, "tcell": 490},
              "mixed (macrophage/T-cell)", 0.51, 3.1),
    "MK447": ({"fibroblast": 0, "macrophage": 75, "tcell": 545},
              "T-cell dominant", 0.879, 3.9),
}


class TestDominanceAgainstSpecTable(unittest.TestCase):
    def test_worked_table(self):
        for sample, (counts, expected_pattern, expected_p, expected_hw_pp) in CASES.items():
            with self.subTest(sample=sample):
                result = compute_dominance(counts)
                self.assertEqual(result["pattern"], expected_pattern,
                                  f"{sample}: pattern mismatch")

                # Recover p / half-width from the basis string for a numeric check.
                basis = result["basis"]
                p_str = basis.split("p=")[1].split(",")[0]
                hw_str = basis.split("half-width=")[1].split("pp")[0]
                self.assertAlmostEqual(float(p_str), expected_p, places=3,
                                        msg=f"{sample}: p mismatch")
                # delta, not decimal places: MK362's true half-width (2.2518pp)
                # sits almost exactly on a .5 rounding boundary between the
                # spec's reported 2.2 and this code's 2.3 -- a tolerance keeps
                # the test about catching real drift, not a rounding artifact.
                self.assertAlmostEqual(float(hw_str), expected_hw_pp, delta=0.15,
                                        msg=f"{sample}: half-width mismatch")

    def test_bw21_fibroblast_eligible_but_not_in_top2(self):
        # BW21 is the one sample where fibroblast clears n>=20; the spec
        # requires the pairwise test still run on the top-2-by-count axes
        # (T-cell, macrophage), not fibroblast.
        result = compute_dominance({"fibroblast": 50, "macrophage": 940, "tcell": 5063})
        self.assertIn("fibroblast", result["eligible"])
        self.assertEqual(result["dominant_axis"], "tcell")
        self.assertFalse(result["three_way_ambiguous"])

    def test_ineligible_axis_reported_not_dropped(self):
        result = compute_dominance({"fibroblast": 3, "macrophage": 163, "tcell": 134})
        self.assertEqual(result["ineligible"], {"fibroblast": 3})
        self.assertNotIn("fibroblast", result["eligible"])


if __name__ == "__main__":
    unittest.main()
