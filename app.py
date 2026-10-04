import sys
import os
import platform
import asyncio
from pathlib import Path

# 1. OS-agnostic path configuration
BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.append(str(BASE_DIR))

# 2. Windows-specific asyncio fix
if platform.system() == "Windows":
    try:
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    except AttributeError:
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
    page_title="ProteinScope Studio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- CUSTOM MODERN STYLING (CSS) ---
st.markdown("""
<style>
    /* Dark Theme Custom Accent Styling */
    .stApp {
        background-color: #0E1117;
    }
    
    /* Header Gradient Banner */
    .header-box {
        background: linear-gradient(135deg, #1E2640 0%, #0F172A 100%);
        border: 1px solid #334155;
        border-radius: 12px;
        padding: 24px;
        margin-bottom: 25px;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
    }
    
    /* Metric Card Styling */
    .metric-card {
        background: #1E293B;
        border-left: 4px solid #38BDF8;
        padding: 16px;
        border-radius: 8px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.2);
    }
    
    /* Sidebar Styling */
    [data-testid="stSidebar"] {
        background-color: #0F172A;
        border-right: 1px solid #1E293B;
    }

    /* Primary Action Buttons */
    .stButton>button {
        border-radius: 8px;
        font-weight: 600;
        transition: all 0.2s ease;
    }
</style>
""", unsafe_html=True)

# --- GLOBAL SESSION STATE ---
SAMPLE_PROTEINS = {
    "Hemoglobin Subunit Alpha (Human)": ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    "Tumor Protein p53 (Human)": ">sp|P04637|P53_HUMAN Cellular tumor antigen p53\nMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD",
    "Insulin (Human)": ">sp|P01308|INS_HUMAN Insulin\nMALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
}

if "sequence" not in st.session_state:
    st.session_state["sequence"] = ""
if "header" not in st.session_state:
    st.session_state["header"] = ""
if "raw_input" not in st.session_state:
    st.session_state["raw_input"] = SAMPLE_PROTEINS["Hemoglobin Subunit Alpha (Human)"]

# --- SIDEBAR NAVIGATION ---
st.sidebar.markdown("## 🧬 **ProteinScope**")
st.sidebar.caption("Bioinformatics Analysis Suite v2.0")
st.sidebar.divider()

navigation = st.sidebar.radio(
    "Select Screen:",
    [
        "🏠 Home & Input",
        "📊 Sequence Analysis",
        "🔮 Structure Prediction",
        "🌐 3D Structure & APIs"
    ]
)

st.sidebar.divider()
# Sidebar Status Widget
if st.session_state["sequence"]:
    st.sidebar.success(f"**Loaded Protein:**\n{st.session_state['header'][:25]}...")
    st.sidebar.info(f"**Length:** {len(st.session_state['sequence'])} amino acids")
else:
    st.sidebar.warning("No sequence loaded")

# --- SCREEN 1: HOME & INPUT ---
if navigation == "🏠 Home & Input":
    st.markdown("""
    <div class="header-box">
        <h1 style="color: #F8FAFC; margin:0;">🧬 ProteinScope Studio</h1>
        <p style="color: #94A3B8; margin-top:8px;">Load, validate, and analyze protein FASTA sequences seamlessly.</p>
    </div>
    """, unsafe_allow_true_style=True, unsafe_allow_html=True)

    col_presets, col_editor = st.columns([1, 2])

    with col_presets:
        st.subheader("📋 Preset Examples")
        st.caption("Select a standard benchmark protein to autofill:")
        selected_preset = st.selectbox("Choose Benchmark:", list(SAMPLE_PROTEINS.keys()))
        
        if st.button("Apply Selected Preset", use_container_width=True):
            st.session_state["raw_input"] = SAMPLE_PROTEINS[selected_preset]
            st.rerun()

    with col_editor:
        st.subheader("📝 Sequence Input")
        raw_input = st.text_area(
            "Enter FASTA or Raw Amino Acid Sequence:",
            value=st.session_state["raw_input"],
            height=220,
            help="Paste a standard FASTA sequence starting with '>' or plain sequence letters."
        )

        col_btn1, col_btn2 = st.columns(2)
        with col_btn1:
            if st.button("🚀 Load & Validate", type="primary", use_container_width=True):
                header, seq = parse_fasta(raw_input)
                is_valid, msg = validate_sequence(seq)

                if is_valid:
                    st.session_state["sequence"] = seq
                    st.session_state["header"] = header
                    st.session_state["raw_input"] = raw_input
                    st.success(f"Successfully loaded **{header}** ({len(seq)} aa)")
                else:
                    st.error(f"Validation Error: {msg}")
        
        with col_btn2:
            if st.button("🧹 Clear Input", use_container_width=True):
                st.session_state["raw_input"] = ""
                st.session_state["sequence"] = ""
                st.session_state["header"] = ""
                st.rerun()

# --- SCREEN 2: SEQUENCE ANALYSIS ---
elif navigation == "📊 Sequence Analysis":
    st.markdown("""
    <div class="header-box">
        <h1 style="color: #F8FAFC; margin:0;">📊 Sequence Composition Analysis</h1>
        <p style="color: #94A3B8; margin-top:8px;">Physicochemical properties & amino acid distribution.</p>
    </div>
    """, unsafe_allow_html=True)

    if not st.session_state["sequence"]:
        st.info("💡 No sequence currently loaded. Please visit **🏠 Home & Input** first.")
    else:
        seq = st.session_state["sequence"]
        props = calculate_sequence_properties(seq)

        st.subheader(f"Target: `{st.session_state['header']}`")
        
        # High-impact Metrics Display
        m1, m2, m3 = st.columns(3)
        m1.metric("Molecular Weight", f"{props['molecular_weight']:.2f} Da")
        m2.metric("Isoelectric Point (pI)", f"{props['isoelectric_point']:.2f}")
        m3.metric("Sequence Length", f"{len(seq)} residues")

        st.divider()

        # Visual Chart Breakdown
        st.subheader("Amino Acid Distribution")
        df_freq = pd.DataFrame(list(props["amino_acid_counts"].items()), columns=["Amino Acid", "Count"])
        
        fig = px.bar(
            df_freq, 
            x="Amino Acid", 
            y="Count", 
            color="Count",
            color_continuous_scale="Viridis",
            title="Residue Count Frequency",
            template="plotly_dark"
        )
        fig.update_layout(height=400, margin=dict(l=20, r=20, t=40, b=20))
        st.plotly_chart(fig, use_container_width=True)

# --- SCREEN 3: STRUCTURE PREDICTION ---
elif navigation == "🔮 Structure Prediction":
    st.markdown("""
    <div class="header-box">
        <h1 style="color: #F8FAFC; margin:0;">🔮 Secondary Structure Prediction</h1>
        <p style="color: #94A3B8; margin-top:8px;">Estimated Alpha-Helix, Beta-Sheet, and Coil compositions.</p>
    </div>
    """, unsafe_allow_html=True)

    if not st.session_state["sequence"]:
        st.info("💡 No sequence currently loaded. Please visit **🏠 Home & Input** first.")
    else:
        seq = st.session_state["sequence"]
        preds = predict_secondary_structure(seq)

        p1, p2, p3 = st.columns(3)
        p1.metric("Helix (H)", f"{preds.get('helix_pct', 0):.1f}%")
        p2.metric("Sheet (E)", f"{preds.get('sheet_pct', 0):.1f}%")
        p3.metric("Coil (C)", f"{preds.get('coil_pct', 0):.1f}%")

        st.divider()
        st.subheader("Secondary Structure Sequence Map")
        st.caption("Visual representation per amino acid position:")
        st.code(preds.get("predicted_sequence", "N/A"), language="text")

# --- SCREEN 4: 3D STRUCTURE & APIS ---
elif navigation == "🌐 3D Structure & APIs":
    st.markdown("""
    <div class="header-box">
        <h1 style="color: #F8FAFC; margin:0;">🌐 Structural & Interaction Databases</h1>
        <p style="color: #94A3B8; margin-top:8px;">Fetch structural entries from RCSB PDB and interaction networks from STRING DB.</p>
    </div>
    """, unsafe_allow_html=True)

    tab1, tab2 = st.tabs(["🏛️ RCSB PDB Loader", "🕸️ STRING Interaction Network"])

    with tab1:
        st.subheader("Fetch RCSB PDB Data")
        col_input, col_action = st.columns([3, 1])
        with col_input:
            pdb_id = st.text_input("Enter 4-Character PDB ID:", value="1TUP").strip()
        with col_action:
            st.write(" ")
            st.write(" ")
            fetch_pdb_btn = st.button("Fetch Structure", use_container_width=True, type="primary")

        if fetch_pdb_btn:
            if pdb_id:
                with st.spinner(f"Retrieving structure `{pdb_id}` from RCSB..."):
                    pdb_data = fetch_pdb_file(pdb_id)

                if pdb_data:
                    st.success(f"PDB structure for **{pdb_id.upper()}** retrieved!")
                    
                    st.download_button(
                        label=f"💾 Download {pdb_id.upper()}.pdb File",
                        data=pdb_data,
                        file_name=f"{pdb_id.lower()}.pdb",
                        mime="chemical/x-pdb",
                        use_container_width=True
                    )
                    
                    with st.expander("Inspect Raw PDB File Headers"):
                        st.code(pdb_data[:1800] + "\n... [truncated]", language="text")
                else:
                    st.error(f"Could not retrieve PDB entry for '{pdb_id}'.")

    with tab2:
        st.subheader("Query STRING Functional Network")
        col_str_input, col_str_btn = st.columns([3, 1])
        with col_str_input:
            protein_name = st.text_input("Enter Protein Name (e.g., TP53, EGFR):", value="TP53").strip()
        with col_str_btn:
            st.write(" ")
            st.write(" ")
            fetch_string_btn = st.button("Query STRING", use_container_width=True, type="primary")

        if fetch_string_btn:
            if protein_name:
                with st.spinner(f"Searching interaction partners for `{protein_name}`..."):
                    interactions = fetch_string_interactions(protein_name)

                if interactions:
                    st.success(f"Found {len(interactions)} functional interaction partner(s)!")
                    df_interactions = pd.DataFrame(interactions)
                    display_cols = [c for c in ['preferredName_A', 'preferredName_B', 'score', 'ncbiTaxonId'] if c in df_interactions.columns]
                    
                    st.dataframe(
                        df_interactions[display_cols] if display_cols else df_interactions,
                        use_container_width=True
                    )
                else:
                    st.error(f"No functional interaction network data found for '{protein_name}'.")