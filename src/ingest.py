"""Step 1: ingest the 8 CellRanger-format PDAC FNA samples into a single AnnData.

RS01 is excluded (ddSeq/Dropseq CSV format, different schema, out of scope per
docs/plans/2026-07-17-pdac-tme-navigator-design.md).
"""

from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse

DATA_DIR = Path(__file__).resolve().parent.parent / "FNA_scRNA _JJL"
OUT_PATH = Path(__file__).resolve().parent.parent / "data" / "processed" / "ingested.h5ad"

SAMPLES = ["MK362", "AM67", "MK364", "MK336", "BW21", "MK371", "MK359", "MK447"]


def _resolve(path: Path, name: str) -> Path:
    """Prefer the .gz variant if present, else the plain file.

    Compression is inconsistent across samples in this dataset (some ship
    both, MK447 ships only .gz) — pandas/scipy both transparently handle
    gzip, so picking whichever exists is enough.
    """
    gz = path / f"{name}.gz"
    return gz if gz.exists() else path / name


def load_sample(sample: str) -> ad.AnnData:
    """Read a CellRanger v3-style mtx triplet directly.

    scanpy's read_10x_mtx assumes v3-format features.tsv implies gzip
    compression, which doesn't hold consistently here — so we parse the
    triplet ourselves instead of fighting its format detection.
    """
    path = DATA_DIR / f"{sample}_filtered_feature_bc_matrix"
    matrix = scipy.io.mmread(_resolve(path, "matrix.mtx")).T.tocsr()
    features = pd.read_csv(_resolve(path, "features.tsv"), header=None, sep="\t")
    barcodes = pd.read_csv(_resolve(path, "barcodes.tsv"), header=None)

    adata = ad.AnnData(X=scipy.sparse.csr_matrix(matrix))
    adata.var_names = pd.Index(features[1].values)
    adata.var["gene_ids"] = features[0].values
    adata.var["feature_types"] = features[2].values
    adata.obs_names = barcodes[0].values

    adata.var_names_make_unique()
    adata.obs["sample"] = sample
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    return adata


def ingest() -> ad.AnnData:
    per_sample = {}
    for sample in SAMPLES:
        adata = load_sample(sample)
        per_sample[sample] = adata
        print(f"{sample}: {adata.n_obs} cells x {adata.n_vars} genes")

    combined = ad.concat(per_sample, label="sample_batch_key", index_unique=None, join="outer")
    combined.obs.drop(columns=["sample_batch_key"], inplace=True)
    print(f"\nCombined: {combined.n_obs} cells x {combined.n_vars} genes across {len(SAMPLES)} samples")
    return combined


if __name__ == "__main__":
    combined = ingest()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(OUT_PATH)
    print(f"\nSaved to {OUT_PATH}")
