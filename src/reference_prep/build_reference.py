"""Assemble the Zenodo 6024273 reference (pk_all.rds, exported by
export_pk_all.R into an mtx triplet + meta.csv) into a single AnnData.

Run after export_pk_all.R:
    Rscript src/reference_prep/export_pk_all.R data/reference/pk_all.rds data/reference/pk_all_export
    python src/reference_prep/build_reference.py
"""

from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
EXPORT_DIR = REPO_ROOT / "data" / "reference" / "pk_all_export"
OUT_PATH = REPO_ROOT / "data" / "reference" / "pk_all.h5ad"


def build_reference() -> ad.AnnData:
    matrix = scipy.io.mmread(EXPORT_DIR / "matrix.mtx").T.tocsr()  # genes x cells -> cells x genes
    genes = pd.read_csv(EXPORT_DIR / "genes.tsv", header=None)[0].values
    barcodes = pd.read_csv(EXPORT_DIR / "barcodes.tsv", header=None)[0].values
    meta = pd.read_csv(EXPORT_DIR / "meta.csv", index_col=0)

    adata = ad.AnnData(X=scipy.sparse.csr_matrix(matrix))
    adata.var_names = pd.Index(genes)
    adata.obs_names = pd.Index(barcodes)
    adata.var_names_make_unique()

    meta = meta.loc[adata.obs_names]
    adata.obs = meta

    print(f"Reference: {adata.n_obs} cells x {adata.n_vars} genes")
    label_cols = [c for c in adata.obs.columns if "type" in c.lower() or "cell" in c.lower()]
    print(f"Candidate label columns: {label_cols}")
    for col in label_cols:
        print(f"\n{col} value counts:")
        print(adata.obs[col].value_counts(dropna=False))
    return adata


if __name__ == "__main__":
    reference = build_reference()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    reference.write_h5ad(OUT_PATH)
    print(f"\nSaved reference AnnData to {OUT_PATH}")
