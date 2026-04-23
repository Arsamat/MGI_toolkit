"""
Shared results visualization for single-cell cNMF.
Used by both the Run cNMF page (after job completes) and the Retrieve Job page.
"""
import io
import json
import zipfile

import pandas as pd
import requests
import streamlit as st

from module_clustering import m_clustering
from module_heatmap import module_heatmap_ui
from make_expression_heatmap import get_expression_heatmap
from hypergeometric import hypergeom_ui


# ------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------
@st.cache_data
def _cached_feather_bytes(df):
    buf = io.BytesIO()
    df.reset_index(drop=False).to_feather(buf)
    buf.seek(0)
    return buf


@st.cache_data
def _convert_for_download(df):
    return df.to_csv().encode("utf-8")


def _toggle(key):
    st.session_state[key] = not st.session_state[key]


# ------------------------------------------------------------------
# Main display function
# ------------------------------------------------------------------
def display_results_sc(meta, metadata_index):
    """Parse the ZIP in session state and render all result visualisations.

    Expects ``st.session_state["cnmf_zip_bytes_sc"]`` to contain the raw ZIP
    bytes.  Does nothing (with an info message) when no ZIP is present.
    """

    zip_bytes = st.session_state.get("cnmf_zip_bytes_sc")
    if not zip_bytes:
        st.info("No results loaded yet. Submit a job or retrieve a previous one.")
        st.stop()

    # ------------------------------------------------------------------
    # Parse ZIP
    # ------------------------------------------------------------------
    zf = zipfile.ZipFile(io.BytesIO(zip_bytes))

    for name in zf.namelist():
        if name.endswith(".pdf"):
            st.session_state["cnmf_pdf_sc"] = zf.read(name)
        elif ".usages" in name.lower():
            df = pd.read_csv(zf.open(name), sep="\t", index_col=0)
            st.session_state["cnmf_module_usages_sc"] = df
            st.session_state["cnmf_module_usages_transposed_sc"] = df.T
        elif ".gene" in name.lower():
            df = pd.read_csv(zf.open(name), sep="\t", index_col=0)
            st.session_state["cnmf_gene_loadings_sc"] = df

    if st.session_state["cnmf_module_usages_sc"] is None or st.session_state["cnmf_gene_loadings_sc"] is None:
        st.error("ZIP loaded but usages / gene loadings not found.")
        st.stop()

    # ------------------------------------------------------------------
    # Results header
    # ------------------------------------------------------------------
    st.divider()
    st.header("Results")

    if st.session_state["cnmf_pdf_sc"] is not None:
        st.download_button(
            "Download Heatmap PDF",
            data=st.session_state["cnmf_pdf_sc"],
            file_name="cnmf_heatmap.pdf",
            mime="application/pdf",
        )

    # ------------------------------------------------------------------
    # Output matrices
    # ------------------------------------------------------------------
    st.subheader("Output Matrices")
    col_a, col_b = st.columns(2)

    with col_a:
        if st.button("Show / Hide Module Usages"):
            _toggle("display_scores_sc")
        if st.session_state.get("display_scores_sc", False):
            st.caption("Module Usage Matrix (cells x programs)")
            st.dataframe(st.session_state["cnmf_module_usages_sc"], use_container_width=True)
            st.download_button(
                "Download CSV",
                data=_convert_for_download(st.session_state["cnmf_module_usages_sc"]),
                file_name="module_usages.csv",
                mime="text/csv",
                icon=":material/download:",
            )

    with col_b:
        if st.button("Show / Hide Gene Loadings"):
            _toggle("display_loadings_sc")
        if st.session_state.get("display_loadings_sc", False):
            st.caption("Gene Spectra Score Matrix — first 20 rows x 20 columns")
            st.dataframe(st.session_state["cnmf_gene_loadings_sc"].iloc[:20, :20], use_container_width=True)
            st.download_button(
                "Download CSV",
                data=_convert_for_download(st.session_state["cnmf_gene_loadings_sc"]),
                file_name="gene_loadings.csv",
                mime="text/csv",
                icon=":material/download:",
            )

    # ------------------------------------------------------------------
    # Preview heatmap
    # ------------------------------------------------------------------
    st.divider()
    st.header("Preview Heatmap")

    module_usages_T = st.session_state["cnmf_module_usages_transposed_sc"]

    with st.form("cnmf_preview_form_sc"):
        annotation_cols = st.multiselect(
            "Metadata columns for annotation",
            options=meta.columns.tolist(),
            key="cnmf_annotation_widget_sc",
        )
        average_groups = st.checkbox("Average groups", value=False)
        submitted = st.form_submit_button("Generate Preview")

    if submitted:
        st.session_state["cnmf_annotations_default_sc"] = annotation_cols
        common_samples = [s for s in meta[metadata_index] if s in module_usages_T.columns]

        if not common_samples:
            st.warning("No overlapping samples between module usages and metadata.")
        else:
            module_bytes = _cached_feather_bytes(st.session_state["cnmf_module_usages_transposed_sc"])
            meta_bytes = _cached_feather_bytes(meta)

            files = {
                "df": ("modules.feather", module_bytes, "application/octet-stream"),
                "metadata": ("meta.feather", meta_bytes, "application/octet-stream"),
            }
            data = {
                "metadata_index": metadata_index,
                "average_groups": average_groups,
                "annotation_cols": json.dumps(annotation_cols),
            }
            resp = requests.post(
                st.session_state["API_URL_sc"] + "/initial_heatmap_preview",
                files=files,
                data=data,
            )
            st.session_state["cnmf_preview_png_sc"] = resp.content

    if st.session_state["cnmf_preview_png_sc"] is not None:
        st.image(st.session_state["cnmf_preview_png_sc"])
        st.download_button(
            "Download PNG",
            data=st.session_state["cnmf_preview_png_sc"],
            file_name="cnmf_preview.png",
            mime="image/png",
        )

    # ------------------------------------------------------------------
    # Sample clustering
    # ------------------------------------------------------------------
    st.divider()
    st.header("Sample Clustering")

    if st.checkbox("Cluster Samples"):
        module_bytes = _cached_feather_bytes(st.session_state["cnmf_module_usages_sc"])
        meta_bytes = _cached_feather_bytes(meta)

        st.subheader("Step 1 — Dendrogram")
        if st.button("Run Sample Dendrogram"):
            files = {
                "module_usages": ("modules.feather", module_bytes, "application/octet-stream"),
                "metadata": ("meta.feather", meta_bytes, "application/octet-stream"),
            }
            resp = requests.post(
                st.session_state["API_URL_sc"] + "/cluster_samples/",
                files=files,
                data={"metadata_index": metadata_index, "k": 0},
            )
            if resp.status_code == 200:
                st.session_state["cnmf_initial_sample_dendogram_sc"] = bytes.fromhex(resp.json()["dendrogram_png"])
            else:
                st.error(resp.text)

        if st.session_state["cnmf_initial_sample_dendogram_sc"] is not None:
            st.image(st.session_state["cnmf_initial_sample_dendogram_sc"])
            st.download_button(
                "Download PNG",
                data=st.session_state["cnmf_initial_sample_dendogram_sc"],
                file_name="cnmf_initial_sample_dendogram.png",
                mime="image/png",
            )

        st.subheader("Step 2 — Assign Clusters")
        k_samples = st.number_input("Number of sample clusters", 2, 10, 3)

        if st.button("Run Sample Clustering"):
            files = {
                "module_usages": ("modules.feather", module_bytes, "application/octet-stream"),
                "metadata": ("meta.feather", meta_bytes, "application/octet-stream"),
            }
            resp = requests.post(
                st.session_state["API_URL_sc"] + "/cluster_samples/",
                files=files,
                data={"metadata_index": metadata_index, "k": str(k_samples)},
            )
            if resp.status_code != 200:
                st.error(resp.text)
            else:
                payload = resp.json()
                st.session_state["cnmf_sample_leaf_order_sc"] = payload["leaf_order"]
                st.session_state["cnmf_sample_cluster_labels_sc"] = payload["cluster_labels"]
                st.session_state["cnmf_sample_order_sc"] = payload["sample_order"]
                st.session_state["cnmf_sample_dendogram_sc"] = bytes.fromhex(payload["dendrogram_png"])
                st.session_state["cnmf_sample_order_heatmap_sc"] = bytes.fromhex(payload["heatmap_png"])

        if st.session_state["cnmf_sample_dendogram_sc"] is not None:
            st.image(st.session_state["cnmf_sample_dendogram_sc"])
            st.download_button(
                "Download PNG",
                data=st.session_state["cnmf_sample_dendogram_sc"],
                file_name="cnmf_sample_dendogram.png",
                mime="image/png",
            )

            df_tmp = st.session_state["cnmf_module_usages_sc"].copy()
            df_tmp = df_tmp.reset_index()
            df_tmp.columns = ["Sample"] + list(df_tmp.columns[1:])
            df_tmp = df_tmp.set_index("Sample")
            sample_order = df_tmp.index.tolist()
            meta_aligned = meta.set_index(metadata_index).loc[sample_order].copy()
            meta_aligned["H_Clustering_Labels"] = st.session_state["cnmf_sample_cluster_labels_sc"]
            st.dataframe(meta_aligned, use_container_width=True)

        st.subheader("Step 3 — Annotated Heatmap")
        if st.session_state.get("cnmf_sample_leaf_order_sc") is not None:
            with st.form("cnmf_annotated_heatmap_form_sc"):
                annotation_cols_annot = st.multiselect(
                    "Annotation columns",
                    ["Cluster"] + meta.columns.tolist(),
                    default=st.session_state.get("cnmf_annotations_default_sc"),
                )
                submit_annot = st.form_submit_button("Generate Annotated Heatmap")

            if submit_annot:
                files = {
                    "module_usages": ("modules.feather", module_bytes, "application/octet-stream"),
                    "metadata": ("meta.feather", meta_bytes, "application/octet-stream"),
                }
                data = {
                    "metadata_index": metadata_index,
                    "leaf_order": json.dumps(st.session_state["cnmf_sample_leaf_order_sc"]),
                    "annotation_cols": json.dumps(annotation_cols_annot),
                    "cluster_labels": json.dumps(st.session_state["cnmf_sample_cluster_labels_sc"]),
                }
                resp = requests.post(
                    st.session_state["API_URL_sc"] + "/annotated_heatmap/",
                    files=files,
                    data=data,
                )
                if resp.status_code == 200:
                    st.session_state["cnmf_sample_order_heatmap_sc"] = resp.content
                else:
                    st.error(resp.text)

            if st.session_state["cnmf_sample_order_heatmap_sc"] is not None:
                st.image(st.session_state["cnmf_sample_order_heatmap_sc"])
                st.download_button(
                    "Download PNG",
                    data=st.session_state["cnmf_sample_order_heatmap_sc"],
                    file_name="cnmf_sample_order_heatmap.png",
                    mime="image/png",
                )

        if st.checkbox("Calculate Hypergeometric Values"):
            hypergeom_ui(
                meta_bytes,
                st.session_state["cnmf_module_usages_sc"],
                st.session_state["cnmf_sample_cluster_labels_sc"],
            )

        if st.checkbox("Cluster Modules"):
            st.subheader("Module Clustering — Step 1: Dendrogram")
            if st.button("Run Module Dendrogram"):
                dendro_png = m_clustering(
                    st.session_state["cnmf_module_usages_sc"],
                    st.session_state["cnmf_sample_order_sc"],
                    0,
                    cnmf=True,
                )
                st.session_state["cnmf_initial_module_dendogram_sc"] = dendro_png

            if st.session_state["cnmf_initial_module_dendogram_sc"]:
                st.image(st.session_state["cnmf_initial_module_dendogram_sc"])
                st.download_button(
                    "Download PNG",
                    data=st.session_state["cnmf_initial_module_dendogram_sc"],
                    file_name="cnmf_initial_module_dendogram.png",
                    mime="image/png",
                )

            st.subheader("Module Clustering — Step 2: Assign Clusters")
            n_mod = st.slider("Number of module clusters", 2, 12, 4)

            if st.button("Run Final Module Clustering"):
                dendro_png = m_clustering(
                    st.session_state["cnmf_module_usages_sc"],
                    st.session_state["cnmf_sample_order_sc"],
                    n_mod,
                    cnmf=True,
                )
                if dendro_png:
                    st.session_state["cnmf_module_dendogram_sc"] = dendro_png

            if st.session_state["cnmf_module_dendogram_sc"]:
                st.image(st.session_state["cnmf_module_dendogram_sc"])
                st.download_button(
                    "Download PNG",
                    data=st.session_state["cnmf_module_dendogram_sc"],
                    file_name="cnmf_module_dendogram.png",
                    mime="image/png",
                )

            module_heatmap_ui(
                meta_bytes,
                st.session_state["cnmf_module_usages_sc"],
                st.session_state["cnmf_sample_order_sc"],
                st.session_state["cnmf_module_leaf_order_sc"],
                st.session_state["cnmf_module_cluster_labels_sc"],
                cnmf=True,
                default_annotations=st.session_state["cnmf_annotations_default_sc"],
            )

    # ------------------------------------------------------------------
    # Top sample ordering
    # ------------------------------------------------------------------
    st.divider()
    st.header("Top Sample Ordering")

    if st.checkbox("Order by Top Samples"):
        with st.form("cnmf_top_samples_form_sc"):
            annotation_cols = st.multiselect(
                "Metadata columns for annotation",
                meta.columns.tolist(),
                default=st.session_state.get("cnmf_annotations_default_sc"),
            )
            submit_top = st.form_submit_button("Generate")

        if submit_top:
            files = {
                "module_usages": (
                    "modules.feather",
                    _cached_feather_bytes(st.session_state["cnmf_module_usages_sc"]),
                    "application/octet-stream",
                ),
                "metadata": (
                    "meta.feather",
                    _cached_feather_bytes(meta),
                    "application/octet-stream",
                ),
            }
            resp = requests.post(
                st.session_state["API_URL_sc"] + "/heatmap_top_samples/",
                files=files,
                data={"metadata_index": metadata_index, "annotation_cols": ",".join(annotation_cols)},
            )
            if resp.status_code == 200:
                st.session_state["cnmf_top_order_heatmap_sc"] = bytes.fromhex(resp.json()["heatmap_png"])
            else:
                st.error(resp.text)

        if st.session_state["cnmf_top_order_heatmap_sc"]:
            st.image(st.session_state["cnmf_top_order_heatmap_sc"])
            st.download_button(
                "Download PNG",
                data=st.session_state["cnmf_top_order_heatmap_sc"],
                file_name="cnmf_top_order_heatmap.png",
                mime="image/png",
            )

    # ------------------------------------------------------------------
    # Gene expression heatmap
    # ------------------------------------------------------------------
    st.divider()
    st.header("Gene Expression Heatmap")

    if st.checkbox("Show Gene Expression Matrix"):
        expr = get_expression_heatmap(
            st.session_state["cnmf_gene_loadings_sc"],
            st.session_state.get("cnmf_annotations_default_sc"),
        )
        if expr is not None:
            st.session_state["cnmf_expression_heatmap_sc"] = expr

        if st.session_state["cnmf_expression_heatmap_sc"] is not None:
            st.image(st.session_state["cnmf_expression_heatmap_sc"])
            st.download_button(
                "Download PNG",
                data=st.session_state["cnmf_expression_heatmap_sc"],
                file_name="cnmf_expression_heatmap.png",
                mime="image/png",
            )

    # ------------------------------------------------------------------
    # Navigation
    # ------------------------------------------------------------------
    st.divider()
    if st.button("Continue to Gene Descriptions"):
        st.session_state["_go_to_sc"] = "Get Gene Descriptions"
        st.rerun()
