"""
H5AD file export for standalone script execution.

This module handles export of single-cell data to H5AD format (AnnData),
including count matrix generation, metadata preparation, and file I/O
for downstream analysis with Scanpy and other tools.
"""

import gzip
import json
import os

try:
    from mfsflow.scripts.path_layout import expression_dir
except ImportError:
    from path_layout import expression_dir


MATRIX_TYPES = (
    "exon.umi",
    "intron.umi",
    "inex.umi",
    "exon.read",
    "intron.read",
    "inex.read",
)


def export_h5ad(out_dir, project, config=None, main_matrix="inex.umi"):
    try:
        import anndata as ad
        import pandas as pd
        import scipy.io
    except ImportError as e:
        raise RuntimeError(
            "H5AD export requires anndata, pandas and scipy. Install requirements.txt first."
        ) from e

    expression_output_dir = expression_dir(out_dir)
    matrices = {}
    genes = None
    barcodes = None

    for matrix_type in MATRIX_TYPES:
        matrix_dir = os.path.join(expression_output_dir, f"{project}.{matrix_type}")
        if not _matrix_dir_exists(matrix_dir):
            continue

        matrix, cur_genes, cur_barcodes = _read_mex(matrix_dir)
        if genes is None:
            genes = cur_genes
            barcodes = cur_barcodes
        else:
            if cur_genes != genes:
                raise ValueError(f"Gene order differs in {matrix_dir}")
            if cur_barcodes != barcodes:
                raise ValueError(f"Barcode order differs in {matrix_dir}")

        matrices[matrix_type.replace(".", "_")] = matrix.T.tocsr()

    if not matrices:
        raise FileNotFoundError(f"No MEX matrices found under {expression_output_dir}")

    main_key = main_matrix.replace(".", "_")
    if main_key not in matrices:
        main_key = "exon_umi" if "exon_umi" in matrices else next(iter(matrices))

    obs = pd.DataFrame(index=barcodes)
    obs.index.name = "barcode"

    var = pd.DataFrame(
        {
            "gene_id": [g[0] for g in genes],
            "gene_name": [g[1] for g in genes],
            "feature_type": [g[2] for g in genes],
        },
        index=[g[0] for g in genes],
    )
    var.index.name = "gene_id"

    adata = ad.AnnData(X=matrices[main_key], obs=obs, var=var)
    for key, matrix in matrices.items():
        adata.layers[key] = matrix

    adata.uns["project"] = project
    adata.uns["main_matrix"] = main_key
    if config is not None:
        adata.uns["pipeline_config"] = _json_safe(config)

    out_path = os.path.join(expression_output_dir, f"{project}.h5ad")
    adata.write_h5ad(out_path, compression="gzip")
    return out_path


def _matrix_dir_exists(matrix_dir):
    return (
        os.path.exists(os.path.join(matrix_dir, "matrix.mtx.gz"))
        and os.path.exists(os.path.join(matrix_dir, "features.tsv.gz"))
        and os.path.exists(os.path.join(matrix_dir, "barcodes.tsv.gz"))
    )


def _read_mex(matrix_dir):
    import scipy.io

    matrix_path = os.path.join(matrix_dir, "matrix.mtx.gz")
    features_path = os.path.join(matrix_dir, "features.tsv.gz")
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")

    matrix = scipy.io.mmread(matrix_path).tocsr()
    genes = []
    with gzip.open(features_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 1:
                parts = [parts[0], parts[0], "Gene Expression"]
            elif len(parts) == 2:
                parts.append("Gene Expression")
            genes.append(tuple(parts[:3]))

    with gzip.open(barcodes_path, "rt") as f:
        barcodes = [line.strip() for line in f if line.strip()]

    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match features/barcodes in {matrix_dir}"
        )

    return matrix, genes, barcodes


def _json_safe(value):
    try:
        json.dumps(value)
        return value
    except TypeError:
        if isinstance(value, dict):
            return {str(k): _json_safe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(v) for v in value]
        return str(value)
