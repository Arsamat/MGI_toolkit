import streamlit as st
import requests, io
import pandas as pd
from ui_theme import apply_custom_theme
from streamlit_autorefresh import st_autorefresh
from .upload_counts_helper_sc import upload_counts
from .results_display_sc import display_results_sc


def run_cnmf_sc():
    # ================================================================
    # APPLY CUSTOM THEME
    # ================================================================
    apply_custom_theme()

    # ================================================================
    # DEFAULT SESSION STATE
    # ================================================================
    DEFAULTS = {
        "job_id_sc": None,
        "cnmf_module_usages_sc": None,
        "cnmf_module_usages_transposed_sc": None,
        "cnmf_gene_loadings_sc": None,
        "cnmf_pdf_sc": None,
        "cnmf_sample_order_sc": None,
        "cnmf_sample_leaf_order_sc": None,
        "cnmf_sample_cluster_labels_sc": None,
        "cnmf_sample_dendogram_sc": None,
        "cnmf_sample_order_heatmap_sc": None,
        "cnmf_initial_sample_dendogram_sc": None,
        "cnmf_module_leaf_order_sc": None,
        "cnmf_module_cluster_labels_sc": None,
        "cnmf_initial_module_dendogram_sc": None,
        "cnmf_module_dendogram_sc": None,
        "cnmf_top_order_heatmap_sc": None,
        "cnmf_expression_heatmap_sc": None,
        "cnmf_preview_png_sc": None,
        "cnmf_annotations_default_sc": None,
        "cnmf_running_sc": False,
        "cnmf_zip_bytes_sc": None,
        "cnmf_sc_job_id": None,
        "display_scores_sc": None,
        "display_loadings_sc": None,
    }
    if "API_URL_sc" not in st.session_state:
        st.session_state["API_URL_sc"] = "http://18.218.84.81:8000/"

    for k, v in DEFAULTS.items():
        st.session_state.setdefault(k, v)

    # ================================================================
    # PAGE HEADER
    # ================================================================
    st.title("Run Single-Cell cNMF")
    st.markdown(
        "Run Consensus NMF (Kotliar et al.) on single-cell RNA-seq data. "
        "Jobs run on the server in the background — you will receive an email when your results are ready."
    )

    # ================================================================
    # REQUIRED INPUT VALIDATION
    # ================================================================
    if "meta_sc" not in st.session_state or st.session_state["meta_sc"] is None:
        st.error("Upload metadata first.")
        st.stop()

    meta = st.session_state["meta_sc"]
    metadata_index = st.session_state.get("metadata_index_sc", "")

    # ================================================================
    # DATA SETUP
    # ================================================================
    st.divider()
    st.header("Data Setup")

    st.subheader("Metadata")
    st.dataframe(meta.head(), use_container_width=True)

    st.subheader("Count Matrix")
    upload_counts()

    # ================================================================
    # PARAMETERS
    # ================================================================
    st.divider()
    st.header("Parameters")

    col1, col2 = st.columns(2)
    with col1:
        k = st.number_input("Number of programs (k)", min_value=2, max_value=50, value=7)
    with col2:
        hvg = st.number_input("Highly variable genes", min_value=100, max_value=20000, value=2000)

    st.session_state["gene_column_sc"] = st.text_input(
        "Gene name column",
        help="Name of the column in your counts file that contains gene identifiers",
    )

    st.session_state["gene_symbols_sc"] = st.checkbox(
        "Gene names are already symbols (not Ensembl IDs)",
        value=False,
    )

    st.subheader("Batch Correction")
    st.session_state["batch_sc"] = st.checkbox(
        "Enable batch correction (Harmony)",
        value=False,
    )
    if st.session_state["batch_sc"]:
        if st.session_state.get("meta_sc") is not None:
            meta_cols = list(st.session_state["meta_sc"].columns)
            st.session_state["batch_covars_sc"] = st.multiselect(
                "Metadata columns to use as batch covariates",
                options=meta_cols,
                default=[],
                help="Select variables that capture batch effects. Each must have values across at least two batches.",
            )
        else:
            st.warning("Upload metadata first to select batch covariates.")

    # ================================================================
    # SUBMIT JOB
    # ================================================================
    st.divider()
    st.header("Submit Job")

    email_input = st.text_input(
        "Email address",
        placeholder="you@example.com",
        help="You will be notified when your job completes or fails.",
    )

    if st.button("Submit cNMF Job", type="primary") and not st.session_state["cnmf_running_sc"]:
        if not email_input.strip():
            st.warning("Please enter your email address before submitting.")
        else:
            api_url = st.session_state["API_URL_sc"]

            meta_buf = io.StringIO()
            meta.to_csv(meta_buf, sep="\t", index=False)
            meta_bytes = meta_buf.getvalue().encode("utf-8")

            files = {
                "metadata": ("metadata.tsv", meta_bytes, "text/plain"),
            }
            data = {
                "k": int(k),
                "hvg": int(hvg),
                "design_factor": "Group",
                "metadata_index": metadata_index,
                "job_id": st.session_state["job_id_sc"],
                "batch_correct": ",".join(st.session_state.get("batch_covars_sc", [])),
                "gene_column": st.session_state["gene_column_sc"],
                "gene_symbols": st.session_state["gene_symbols_sc"],
                "email": email_input.strip(),
            }

            try:
                resp = requests.post(f"{api_url}cnmf_sc_submit", files=files, data=data, timeout=30)
                resp.raise_for_status()
                sc_job_id = resp.json()["sc_job_id"]
                st.session_state["cnmf_sc_job_id"] = sc_job_id
                st.session_state["cnmf_running_sc"] = True
                st.rerun()
            except Exception as e:
                st.error(f"Failed to submit job: {e}")

    if st.session_state.get("cnmf_sc_job_id"):
        st.info(f"**Job ID:** `{st.session_state['cnmf_sc_job_id']}` — save this to retrieve your results later.")

    # ================================================================
    # POLL JOB STATUS
    # ================================================================
    if st.session_state["cnmf_running_sc"]:
        sc_job_id = st.session_state["cnmf_sc_job_id"]
        st_autorefresh(interval=10000, key="cnmf_autorefresh_sc")

        try:
            resp = requests.get(f"{st.session_state['API_URL_sc']}job_status/{sc_job_id}", timeout=10)
            status_data = resp.json()
            status = status_data.get("status")

            if status == "completed":
                result_resp = requests.get(
                    f"{st.session_state['API_URL_sc']}job_results/{sc_job_id}", timeout=10
                )
                download_url = result_resp.json()["download_url"]
                zip_resp = requests.get(download_url, timeout=120)
                st.session_state["cnmf_zip_bytes_sc"] = zip_resp.content
                st.session_state["cnmf_running_sc"] = False
                st.success("cNMF completed! Results loaded.")
                st.rerun()
            elif status == "failed":
                st.session_state["cnmf_running_sc"] = False
                st.error(f"Job failed: {status_data.get('error', 'Unknown error')}")
            else:
                st.info("Job is running on the server. You will receive an email when it finishes — you can safely close this tab.")
        except Exception as e:
            st.warning(f"Could not check job status: {e}")

    # ================================================================
    # RESULTS (shared visualization)
    # ================================================================
    display_results_sc(meta, metadata_index)
