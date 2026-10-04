import sys
import os
import platform
import asyncio
import warnings
from pathlib import Path

# 1. OS-agnostic path configuration
BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.append(str(BASE_DIR))

# 2. Windows-specific asyncio fix (silences Python 3.14+ deprecation warnings)
if platform.system() == "Windows":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        try:
            asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
        except (AttributeError, RuntimeError):
            pass

import streamlit as st
import pandas as pd
import plotly.express as px

# Imports from src package
from src.sequence import parse_fasta, validate_sequence, calculate_sequence_properties
from src.prediction import predict_secondary_structure
from src.structure_api import fetch_pdb_file, fetch_string_interactions

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="ProteinScope",
    page_icon="🧬",
    layout="wide"
)

# --- GLOBAL SESSION STATE ---
DEFAULT_FASTA = ">Sample_Protein\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"

if "sequence" not in st.session_state:
    st.session_state["sequence"] = ""
if "header" not in st.session_state:
    st.session_state["header"] = ""
if "raw_input" not in st.session_state:
    st.session_state["raw_input"] = DEFAULT_FASTA

# --- SIDEBAR NAVIGATION ---
st.sidebar.title("🧬 ProteinScope")
navigation = st.sidebar.radio(
    "Navigation",
    [
        "1. 🏠 Home & Input",
        "2. 📊 Sequence Analysis",
        "3. 🔮 Secondary Structure Prediction",
        "4. 🌐 3D Structure & APIs"
    ]
)

# --- SCREEN 1: HOME & INPUT ---
if navigation == "1. 🏠 Home & Input":
    st.title("🧬 ProteinScope: Protein Analysis Suite")
    st.markdown("Welcome! Paste a FASTA sequence or raw amino acid sequence below to get started.")

    raw_input = st.text_area(
        "Enter FASTA or Raw Sequence:",
        value=st.session_state["raw_input"],
        height=200
    )

    if st.button("Load & Validate Sequence"):
        header, seq = parse_fasta(raw_input)
        is_valid, msg = validate_sequence(seq)

        if is_valid:
            st.session_state["sequence"] = seq
            st.session_state["header"] = header
            st.session_state["raw_input"] = raw_input
            st.success(f"Sequence successfully loaded! Header: **{header}** (Length: {len(seq)} aa)")
        else:
            st.error(f"Validation Error: {msg}")

# --- SCREEN 2: SEQUENCE ANALYSIS ---
elif navigation == "2. 📊 Sequence Analysis":
    st.title("📊 Sequence Composition & Properties")

    if not st.session_state["sequence"]:
        st.info("No sequence loaded. Please return to 🏠 Home & Input first.")
    else:
        seq = st.session_state["sequence"]
        props = calculate_sequence_properties(seq)

        st.subheader(f"Protein: {st.session_state['header']}")

        col1, col2, col3 = st.columns(3)
        col1.metric("Molecular Weight", f"{props['molecular_weight']:.2f} Da")
        col2.metric("Isoelectric Point (pI)", f"{props['isoelectric_point']:.2f}")
        col3.metric("Length", f"{len(seq)} aa")

        st.subheader("Amino Acid Frequency")
        df_freq = pd.DataFrame(list(props["amino_acid_counts"].items()), columns=["Amino Acid", "Count"])
        fig = px.bar(df_freq, x="Amino Acid", y="Count", title="Amino Acid Distribution", color="Amino Acid")
        st.plotly_chart(fig, use_container_width=True)

# --- SCREEN 3: SECONDARY STRUCTURE PREDICTION ---
elif navigation == "3. 🔮 Secondary Structure Prediction":
    st.title("🔮 Secondary Structure Prediction")

    if not st.session_state["sequence"]:
        st.info("No sequence loaded. Please return to 🏠 Home & Input first.")
    else:
        seq = st.session_state["sequence"]
        preds = predict_secondary_structure(seq)

        st.subheader("Predicted Secondary States")
        c1, c2, c3 = st.columns(3)
        c1.metric("Helix (H)", f"{preds.get('helix_pct', 0):.1f}%")
        c2.metric("Sheet (E)", f"{preds.get('sheet_pct', 0):.1f}%")
        c3.metric("Coil (C)", f"{preds.get('coil_pct', 0):.1f}%")

        st.subheader("Predicted Sequence Map")
        st.code(preds.get("predicted_sequence", "N/A"), language="text")

# --- SCREEN 4: 3D STRUCTURE & APIS ---
elif navigation == "4. 🌐 3D Structure & APIs":
    st.title("🌐 Structural & Interaction Analysis")

    # Section 1: RCSB PDB
    st.header("1. Fetch PDB Structure")
    pdb_id = st.text_input("Enter 4-character PDB ID:", value="1TUP").strip()

    if st.button("Fetch PDB Data"):
        if pdb_id:
            with st.spinner(f"Fetching PDB data for '{pdb_id}'..."):
                pdb_data = fetch_pdb_file(pdb_id)

            if pdb_data:
                st.success(f"Successfully loaded PDB data for {pdb_id.upper()}!")
                with st.expander("Show Raw PDB Header & Data"):
                    st.code(pdb_data[:1500] + "\n... [truncated]", language="text")

                st.download_button(
                    label=f"Download {pdb_id.upper()}.pdb",
                    data=pdb_data,
                    file_name=f"{pdb_id.lower()}.pdb",
                    mime="chemical/x-pdb"
                )
            else:
                st.error(f"Could not find PDB entry for ID '{pdb_id}'. Please verify the ID on RCSB.org.")
        else:
            st.warning("Please enter a valid PDB ID.")

    st.divider()

    # Section 2: STRING DB
    st.header("2. Functional Interaction Network (STRING DB)")
    protein_name = st.text_input("Enter Protein Name (e.g., TP53, EGFR):", value="TP53").strip()

    if st.button("Fetch Interactions"):
        if protein_name:
            with st.spinner(f"Querying STRING DB for '{protein_name}'..."):
                interactions = fetch_string_interactions(protein_name)

            if interactions:
                st.success(f"Found {len(interactions)} interaction partner(s)!")
                df_interactions = pd.DataFrame(interactions)
                display_cols = [c for c in ['preferredName_A', 'preferredName_B', 'score', 'ncbiTaxonId'] if c in df_interactions.columns]
                st.dataframe(df_interactions[display_cols] if display_cols else df_interactions, use_container_width=True)
            else:
                st.error(f"No functional interaction network data found for '{protein_name}'.")
        else:
            st.warning("Please enter a protein name.")
