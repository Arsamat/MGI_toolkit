import streamlit as st
import requests

from ui_theme import apply_custom_theme
from .results_display_sc import display_results_sc


# Session state keys needed for results visualization
_DEFAULTS = {
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
    "cnmf_zip_bytes_sc": None,
    "cnmf_sc_job_id": None,
    "display_scores_sc": None,
    "display_loadings_sc": None,
}


def retrieve_job_sc():
    apply_custom_theme()

    if "API_URL_sc" not in st.session_state:
        st.session_state["API_URL_sc"] = "http://18.218.84.81:8000/"

    for k, v in _DEFAULTS.items():
        st.session_state.setdefault(k, v)

    # ----------------------------------------------------------------
    # Page header
    # ----------------------------------------------------------------
    st.title("Retrieve Job Results")
    st.markdown(
        "Enter a Job ID from a previous cNMF submission to load and visualize your results. "
        "Results are available for 24 hours after job completion."
    )

    # ----------------------------------------------------------------
    # Metadata check
    # ----------------------------------------------------------------
    if "meta_sc" not in st.session_state or st.session_state["meta_sc"] is None:
        st.error("Upload metadata first (use the **Upload Metadata** page).")
        st.stop()

    meta = st.session_state["meta_sc"]
    metadata_index = st.session_state.get("metadata_index_sc", "")

    # ----------------------------------------------------------------
    # Job retrieval
    # ----------------------------------------------------------------
    st.divider()
    st.header("Load Job")

    retrieve_id = st.text_input("Job ID", key="retrieve_job_id_input")

    if st.button("Load Results", type="primary", key="retrieve_job_btn"):
        if not retrieve_id.strip():
            st.warning("Please enter a Job ID.")
        else:
            api_url = st.session_state["API_URL_sc"]
            try:
                resp = requests.get(f"{api_url}job_status/{retrieve_id.strip()}", timeout=10)
                if resp.status_code == 404:
                    st.error("Job not found. Check your Job ID.")
                else:
                    status_data = resp.json()
                    status = status_data.get("status")
                    if status == "running":
                        st.info("Job is still running. Check back later.")
                    elif status == "failed":
                        st.error(f"Job failed: {status_data.get('error', 'Unknown error')}")
                    elif status == "completed":
                        result_resp = requests.get(f"{api_url}job_results/{retrieve_id.strip()}", timeout=10)
                        download_url = result_resp.json()["download_url"]
                        zip_resp = requests.get(download_url, timeout=120)
                        st.session_state["cnmf_zip_bytes_sc"] = zip_resp.content
                        st.session_state["cnmf_sc_job_id"] = retrieve_id.strip()
                        st.success("Results loaded successfully!")
                        st.rerun()
            except Exception as e:
                st.error(f"Error retrieving job: {e}")

    if st.session_state.get("cnmf_sc_job_id"):
        st.info(f"**Loaded Job:** `{st.session_state['cnmf_sc_job_id']}`")

    # ----------------------------------------------------------------
    # Results visualization (shared with Run cNMF page)
    # ----------------------------------------------------------------
    display_results_sc(meta, metadata_index)
