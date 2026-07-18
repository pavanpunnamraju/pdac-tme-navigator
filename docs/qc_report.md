# QC Report (Step 2)

Filter applied: per-sample adaptive mito cutoff (median + 3 MAD, raw pct scale) + `n_genes_by_counts >= 250`. No doublet-score filtering.

## Doublet detection limitation

> Doublet detection (scrublet, run per-sample) found no separable doublet population in any of the 8 samples — score distributions were unimodal with no bimodal shoulder, even in smaller samples. This is a known scrublet limitation for homotypic doublets and low cell-type-diversity samples, not evidence the data is doublet-free. Reported per-sample doublet counts (0-8 cells) reflect an auto-threshold artifact and should not be read as a doublet-burden estimate. Undetected doublets may contribute noise to per-patient composition estimates; they do not affect classifier training, which uses an external reference.

## Before / after counts

| Sample | Before | After | Dropped |
|---|---|---|---|
| MK362 | 7479 | 6140 | 1339 |
| AM67 | 2773 | 1946 | 827 |
| MK364 | 3087 | 2178 | 909 |
| MK336 | 1618 | 1361 | 257 |
| BW21 | 10927 | 9363 | 1564 |
| MK371 | 4188 | 3437 | 751 |
| MK359 | 2093 | 1680 | 413 |
| MK447 (low-confidence, depth-limited) | 1465 | 699 | 766 |
| **Total** | **33630** | **26804** | **6826** |
