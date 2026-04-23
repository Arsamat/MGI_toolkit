"""
Standalone script for running a single-cell cNMF job.
Invoked via subprocess.Popen from single_cell_async_util so it runs
completely outside the uvicorn process — no lock inheritance, no memory sharing.
"""
import argparse
import io
import os
import shutil
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd
from pathlib import Path

# Make backend/ importable regardless of cwd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nmf.cNMF import cNMF_consensus

_RUNNER_DIR = Path(__file__).resolve().parent


def _apply_gene_mapping(df_expr, mouse=False):
    """Map Ensembl IDs → gene symbols on the DataFrame index."""
    if mouse:
        import mygene
        mg = mygene.MyGeneInfo()
        res = mg.querymany(
            list(df_expr.index.astype(str)),
            scopes="ensembl.gene",
            fields="symbol",
            species="mouse",
            as_dataframe=True,
            returnall=False,
            verbose=False,
        )
        if isinstance(res, pd.DataFrame) and "symbol" in res.columns:
            sym = res["symbol"].dropna()
            sym = sym[~sym.index.duplicated(keep="first")]
            raw_mapping = sym.to_dict()
        else:
            raw_mapping = {}
    else:
        genes = pd.read_csv(_RUNNER_DIR / "gene_maps.csv")
        raw_mapping = dict(zip(genes["Gene stable ID"], genes["Gene name"]))

    seen = set()
    result = []
    for val in df_expr.index:
        mapped = raw_mapping.get(val, val)
        if not mapped or str(mapped).lower() == "nan":
            mapped = val
        if mapped in seen:
            mapped = val
        seen.add(mapped)
        result.append(mapped)

    df_expr.index = result
    return df_expr
from infra.s3_utils import s3, BUCKET
from infra.job_store import update_job
from infra.email_utils import send_job_complete, send_job_failed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sc_job_id",      required=True)
    parser.add_argument("--job_id",         required=True)
    parser.add_argument("--k",              type=int, required=True)
    parser.add_argument("--hvg",            type=int, required=True)
    parser.add_argument("--design_factor",  required=True)
    parser.add_argument("--metadata_index", required=True)
    parser.add_argument("--batch_vars",     default="")   # comma-separated or empty
    parser.add_argument("--gene_column",    required=True)
    parser.add_argument("--gene_symbols",   required=True)  # "True" / "False"
    parser.add_argument("--email",          required=True)
    parser.add_argument("--metadata_path",  required=True)  # temp file written by caller
    args = parser.parse_args()

    sc_job_id     = args.sc_job_id
    job_id        = args.job_id
    k             = args.k
    hvg           = args.hvg
    design_factor = args.design_factor
    metadata_index = args.metadata_index
    batch_vars    = [v for v in args.batch_vars.split(",") if v]
    gene_column   = args.gene_column
    gene_symbols  = args.gene_symbols.lower() == "true"
    email         = args.email
    metadata_path = args.metadata_path

    tmp_dir = tempfile.mkdtemp(prefix="preprocess_")
    out_dir = tempfile.mkdtemp(prefix="cnmf_")

    try:
        print(f"[sc_job_runner] Starting job {sc_job_id}", flush=True)

        # Download counts from S3
        s3_key = f"jobs/{job_id}/counts.csv"
        counts_path = os.path.join(tmp_dir, "counts.csv")
        print(f"[sc_job_runner] Downloading {s3_key}", flush=True)
        s3.download_file(BUCKET, s3_key, counts_path)

        # Filter cells and genes
        counts = pd.read_csv(counts_path)
        counts = counts.set_index(gene_column)

        min_umis_per_cell   = 1000
        gene_detect_frac_den = 500

        umi_per_cell = counts.sum(axis=0)
        keep_cells   = umi_per_cell >= min_umis_per_cell
        counts_f1    = counts.loc[:, keep_cells]
        print(f"[sc_job_runner] Cells after UMI filter: {counts_f1.shape[1]}", flush=True)

        n_cells             = counts_f1.shape[1]
        min_cells_detected  = int(np.ceil(n_cells / gene_detect_frac_den))
        detected_per_gene   = (counts_f1 > 0).sum(axis=1)
        keep_genes          = detected_per_gene >= min_cells_detected
        df_expr             = counts_f1.loc[keep_genes, :]
        print(f"[sc_job_runner] Genes after detection filter: {df_expr.shape[0]}", flush=True)

        if not gene_symbols:
            print("[sc_job_runner] Translating gene IDs", flush=True)
            df_expr = _apply_gene_mapping(df_expr)

        # Set index to gene names
        for col in [gene_column, "Geneid", "gene_name", "Unnamed: 0"]:
            if col in df_expr.columns:
                df_expr = df_expr.set_index(col)
                break
        df_expr.index.name = None

        metadata = pd.read_csv(metadata_path, sep="\t")

        print(f"[sc_job_runner] Running cNMF_consensus k={k}", flush=True)
        cNMF_consensus(
            k=k,
            hvg=hvg,
            df=df_expr,
            metadata=metadata,
            metadata_index=metadata_index,
            design_factor=design_factor,
            out_dir=out_dir,
            batch_vars=batch_vars,
            single_cell=True,
        )

        sc_run_name       = "BatchCorrected" if batch_vars else "cNMF"
        module_usage_path = os.path.join(out_dir, sc_run_name, f"{sc_run_name}.usages.k_{k}.dt_0_2.consensus.txt")
        gene_zscore_path  = os.path.join(out_dir, sc_run_name, f"{sc_run_name}.gene_spectra_score.k_{k}.dt_0_2.txt")

        out_zip = os.path.join(out_dir, "cnmf_bundle.zip")
        with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as z:
            for p in [module_usage_path, gene_zscore_path]:
                if os.path.isfile(p):
                    z.write(p, arcname=os.path.basename(p))
                else:
                    print(f"[sc_job_runner] WARNING: expected file not found: {p}", flush=True)

        s3_result_key = f"sc-cnmf-jobs/{sc_job_id}/results.zip"
        print(f"[sc_job_runner] Uploading results to s3://{BUCKET}/{s3_result_key}", flush=True)
        s3.upload_file(out_zip, BUCKET, s3_result_key)

        update_job(sc_job_id, "completed", s3_result_key=s3_result_key)
        send_job_complete(email, sc_job_id)
        print(f"[sc_job_runner] Job {sc_job_id} completed successfully", flush=True)

    except Exception as e:
        import traceback
        print(f"[sc_job_runner] ERROR: {e}", flush=True)
        traceback.print_exc()
        update_job(sc_job_id, "failed", error=str(e))
        try:
            send_job_failed(email, sc_job_id)
        except Exception as email_err:
            print(f"[sc_job_runner] Failed to send failure email: {email_err}", flush=True)

    finally:
        for d in [tmp_dir, out_dir]:
            if os.path.exists(d):
                shutil.rmtree(d)
        # Clean up the temp metadata file written by the caller
        if os.path.exists(metadata_path):
            os.remove(metadata_path)
        print("[sc_job_runner] Cleanup done", flush=True)


if __name__ == "__main__":
    main()
