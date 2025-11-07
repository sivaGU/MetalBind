"""
GPCR Amino Acid Residue Analysis Explorer
A Streamlit application for exploring 369 GPCR receptor analysis results.
"""
import streamlit as st
import pandas as pd
import json
import py3Dmol
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from pathlib import Path
import utils
import re
from collections import defaultdict

# Page configuration
st.set_page_config(
    page_title="GPCR Amino Acid Residue Analysis Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .sidebar-header {
        font-size: 1.2rem;
        font-weight: bold;
        color: #1f77b4;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

@st.cache_data
def load_receptor_data():
    """Load all receptor data once and cache it."""
    return utils.collect_receptors()

def extract_ligand_info(pdb_file_path):
    """
    Extract ligand information from PDB file.
    Returns a dictionary with ligand properties.
    """
    try:
        with open(pdb_file_path, 'r') as f:
            lines = f.readlines()
        
        ligand_info = {
            'name': 'Unknown',
            'formula': 'Unknown',
            'weight': 'Unknown',
            'atom_count': 0,
            'bond_count': 0,
            'smiles': None
        }
        
        atoms = []
        atom_types = defaultdict(int)
        
        for line in lines:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                # Extract atom information
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                element = line[76:78].strip()
                
                if element:
                    atom_types[element] += 1
                    atoms.append({
                        'name': atom_name,
                        'res': res_name,
                        'element': element
                    })
        
        # Count atoms and estimate bonds
        ligand_info['atom_count'] = len(atoms)
        ligand_info['bond_count'] = len(atoms) - 1  # Rough estimate
        
        # Try to get ligand name from HET records
        for line in lines:
            if line.startswith('HET '):
                ligand_info['name'] = line[7:10].strip()
                break
        
        # Calculate molecular weight (rough estimate)
        total_weight = 0
        for element, count in atom_types.items():
            # Rough atomic weights
            weights = {
                'C': 12.01, 'H': 1.01, 'N': 14.01, 'O': 16.00,
                'S': 32.07, 'P': 30.97, 'F': 19.00, 'CL': 35.45,
                'BR': 79.90, 'I': 126.90
            }
            total_weight += weights.get(element.upper(), 12.01) * count
        
        if total_weight > 0:
            ligand_info['weight'] = f"{total_weight:.1f} g/mol"
        
        # Generate chemical formula
        formula_parts = []
        for element, count in sorted(atom_types.items()):
            if count == 1:
                formula_parts.append(element)
            else:
                formula_parts.append(f"{element}{count}")
        ligand_info['formula'] = ''.join(formula_parts)
        
        return ligand_info
        
    except Exception as e:
        st.error(f"Error extracting ligand info: {e}")
        return None

def main():
    # Load data
    receptors = load_receptor_data()
    
    if not receptors:
        st.error("No receptor data found. Please check the data directory path in utils.py")
        return
    
    # Header
    st.markdown('<h1 class="main-header">GPCR Amino Acid Residue Analysis Explorer</h1>', unsafe_allow_html=True)
    st.markdown(f"**Exploring {len(receptors)} GPCR receptor structures with comprehensive pocket analysis**")
    
    # Sidebar
    st.sidebar.markdown('<div class="sidebar-header">Receptor Selection</div>', unsafe_allow_html=True)
    
    # Receptor selection
    receptor_ids = sorted(receptors.keys())
    selected_receptor = st.sidebar.selectbox(
        "Select a receptor:",
        options=receptor_ids,
        index=0
    )
    
    # Get selected receptor data
    receptor_files = receptors[selected_receptor]
    
    # Sidebar info
    st.sidebar.markdown("---")
    st.sidebar.markdown('<div class="sidebar-header">Receptor Info</div>', unsafe_allow_html=True)
    
    # Load summary stats
    summary = utils.get_receptor_summary_stats(selected_receptor, receptor_files)
    
    st.sidebar.markdown(f"**Receptor ID:** {selected_receptor}")
    st.sidebar.markdown(f"**Pocket Residues:** {summary['pocket_residue_count']}")
    st.sidebar.markdown(f"**Conservation Score:** {summary['conservation_score']:.1f}")
    st.sidebar.markdown(f"**Hydrophobic:** {summary['hydrophobic_pct']:.1f}%")
    st.sidebar.markdown(f"**Polar:** {summary['polar_pct']:.1f}%")
    
    # Main content tabs
    tabs = st.tabs(["3D Receptor Visualization", "Pocket Analysis", "Conservation", "Data Download"])
    
    # Tab 1: 3D Structure
    with tabs[0]:
        st.subheader("3D Receptor Visualization")
        if receptor_files["pdb_receptor"] and receptor_files["pdb_receptor"].exists():
            col1, col2 = st.columns([2, 1])
            with col1:
                # Create 3D viewer
                view = py3Dmol.view(width=700, height=500)
                with open(receptor_files["pdb_receptor"], 'r') as f:
                    receptor_pdb = f.read()
                view.addModel(receptor_pdb, "pdb")
                view.setStyle({"cartoon": {"color": "spectrum"}})
                view.zoomTo()
                st.components.v1.html(view._make_html(), height=500, scrolling=False)
            
            with col2:
                # Show receptor info above ligand info
                st.markdown("**Receptor Info**")
                st.markdown(f"Receptor ID: {selected_receptor}")
                # Show ligand information if available
                ligand_info = None
                if receptor_files["pdb_ligand"] and receptor_files["pdb_ligand"].exists():
                    ligand_info = extract_ligand_info(receptor_files["pdb_ligand"])
                st.markdown("**Ligand Information**")
                if ligand_info:
                    st.markdown(f"**Chemical Formula:** {ligand_info.get('formula', 'Unknown')}")
                    st.markdown(f"**Molecular Weight:** {ligand_info.get('weight', 'Unknown')}")
                    st.markdown(f"**Number of Atoms:** {ligand_info.get('atom_count', 0)}")
                    st.markdown(f"**Number of Bonds:** {ligand_info.get('bond_count', 0)}")
                else:
                    st.info("No ligand information available.")
        else:
            st.error("Receptor PDB file not found.")
    
    # Tab 2: Pocket Analysis
    with tabs[1]:
        st.subheader("Pocket Analysis")
        
        # Load pocket data
        pocket_data = utils.load_pocket_data(receptor_files)
        
        if pocket_data is not None and not pocket_data.empty:
            col1, col2 = st.columns([2, 1])
            
            with col1:
                st.markdown("**Pocket Residues Table**")
                # Remove 'is_conserved_position' column if present
                pocket_data_display = pocket_data.drop(columns=['is_conserved_position'], errors='ignore')
                # Remove parenthesis and content from 'label' column for display only
                if 'label' in pocket_data_display.columns:
                    pocket_data_display = pocket_data_display.copy()
                    pocket_data_display['label'] = pocket_data_display['label'].str.replace(r'\s*\(.*\)', '', regex=True)
                st.dataframe(pocket_data_display, use_container_width=True, height=400)
                
                # Download button for pocket residues
                csv_data = pocket_data.to_csv(index=False)
                st.download_button(
                    label="Download Pocket Residues (CSV)",
                    data=csv_data,
                    file_name=f"{selected_receptor}_pocket_residues.csv",
                    mime="text/csv"
                )
            
            with col2:
                # Pocket analysis metrics
                analysis = utils.load_pocket_analysis(receptor_files)
                if analysis:
                    st.markdown("**Pocket Composition**")
                    comp = analysis.get("composition_by_type", {})
                    comp_pct = analysis.get("composition_by_type_pct", {})
                    
                    st.markdown(f"**Hydrophobic:** {comp.get('hydrophobic', 0)} residues ({comp_pct.get('hydrophobic', 0):.1f}%)")
                    st.markdown(f"**Polar:** {comp.get('polar', 0)} residues ({comp_pct.get('polar', 0):.1f}%)")
                    st.markdown(f"**Acidic:** {comp.get('acidic', 0)} residues ({comp_pct.get('acidic', 0):.1f}%)")
                    st.markdown(f"**Basic:** {comp.get('basic', 0)} residues ({comp_pct.get('basic', 0):.1f}%)")
                    
                    # Removed Distance Statistics display (min, max, mean)
                    # st.markdown("**Distance Statistics**")
                    # dist_stats = analysis.get("distance_statistics", {})
                    # st.markdown(f"**Min:** {dist_stats.get('min', 0):.2f} Å")
                    # st.markdown(f"**Max:** {dist_stats.get('max', 0):.2f} Å")
                    # st.markdown(f"**Mean:** {dist_stats.get('mean', 0):.2f} Å")
                    
                    # Pocket analysis download
                    pocket_analysis = utils.load_pocket_analysis(receptor_files)
                    if pocket_analysis:
                        json_data = json.dumps(pocket_analysis, indent=2)
                        st.download_button(
                            label="Download Pocket Analysis (TXT)",
                            data=json_data,
                            file_name=f"{selected_receptor}_pocket_analysis.txt",
                            mime="text/plain"
                        )
                else:
                    st.info("No pocket analysis data available.")
        else:
            st.error("No pocket data found for this receptor.")
    
    # Tab 3: Conservation Analysis
    with tabs[2]:
        st.subheader("Conservation Analysis")
        
        # Load conservation data
        conservation_data = utils.load_conservation_summary(receptor_files)
        
        if conservation_data:
            col1, col2 = st.columns([2, 1])
            
            with col1:
                # Conservation score distribution
                
                # Load per-residue conservation scores from pocket data
                cons_data = utils.load_pocket_data(receptor_files)
                if cons_data is not None and 'conservation_score' in cons_data.columns:
                    score_counts = cons_data['conservation_score'].value_counts().sort_index()
                    pie_df = pd.DataFrame({'Conservation Score': score_counts.index, 'Count': score_counts.values})
                    fig = px.pie(
                        pie_df,
                        names='Conservation Score',
                        values='Count',
                        title="Conservation Score Distribution (Pie Chart)",
                        color_discrete_sequence=px.colors.qualitative.Set3
                    )
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info("No conservation data available for plotting.")
            
            with col2:
                st.markdown("**Top Conserved Residues**")
                if "top_conserved_residues" in conservation_data:
                    top_residues = conservation_data["top_conserved_residues"]
                    if top_residues:
                        top_df = pd.DataFrame(top_residues)
                        st.dataframe(top_df, use_container_width=True, height=300)
                    else:
                        st.info("No conserved residues data available.")
                
                # Conservation summary metrics
                st.markdown("**Conservation Summary**")
                st.markdown(f"**Average Score:** {conservation_data.get('average_conservation_score', 0):.1f}")
                st.markdown(f"**Total Residues:** {conservation_data.get('total_residues', 0)}")
                
                # Conservation by level
                cons_by_level = conservation_data.get('conservation_by_level', {})
                if cons_by_level:
                    st.markdown("**Conservation by Level:**")
                    for level, count in cons_by_level.items():
                        st.markdown(f"**{level.title()}:** {count} residues")
                
                # Conservation data download
                if conservation_data:
                    cons_json_data = json.dumps(conservation_data, indent=2)
                    st.download_button(
                        label="Download Conservation Summary (TXT)",
                        data=cons_json_data,
                        file_name=f"{selected_receptor}_conservation_summary.txt",
                        mime="text/plain",
                        key=f"download_conservation_summary_{selected_receptor}"
                    )
        else:
            st.error("No conservation data found for this receptor.")
    
    # Tab 4: Data Download
    with tabs[3]:
        st.subheader("Data Download")
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.markdown("**Download Options**")
            
            # Pocket residues CSV
            # if receptor_files["pocket_csv"] and receptor_files["pocket_csv"].exists():
            #     with open(receptor_files["pocket_csv"], 'r') as f:
            #         csv_data = f.read()
            #     st.download_button(
            #         label="Download Pocket Residues (CSV)",
            #         data=csv_data,
            #         file_name=f"{selected_receptor}_pocket_residues.csv",
            #         mime="text/csv",
            #         key=f"download_pocket_csv_{selected_receptor}"
            #     )
            
            # Conservation CSV (renamed to Pocket Residue Data (CSV))
            if receptor_files["pocket_cons_csv"] and receptor_files["pocket_cons_csv"].exists():
                with open(receptor_files["pocket_cons_csv"], 'r') as f:
                    cons_csv_data = f.read()
                st.download_button(
                    label="Pocket Residue Data (CSV)",
                    data=cons_csv_data,
                    file_name=f"{selected_receptor}_conservation_data.csv",
                    mime="text/csv",
                    key=f"download_cons_csv_{selected_receptor}"
                )
            
            # Pocket analysis TXT (was JSON)
            if receptor_files["pocket_json"] and receptor_files["pocket_json"].exists():
                with open(receptor_files["pocket_json"], 'r') as f:
                    analysis_data = f.read()
                # Pretty-print JSON as TXT
                try:
                    analysis_txt = json.dumps(json.loads(analysis_data), indent=2)
                except Exception:
                    analysis_txt = analysis_data
                st.download_button(
                    label="Download Pocket Analysis (TXT)",
                    data=analysis_txt,
                    file_name=f"{selected_receptor}_pocket_analysis.txt",
                    mime="text/plain",
                    key=f"download_pocket_txt_{selected_receptor}"
                )
            
            # Conservation summary TXT (was JSON)
            if receptor_files["cons_json"] and receptor_files["cons_json"].exists():
                with open(receptor_files["cons_json"], 'r') as f:
                    cons_data = f.read()
                # Pretty-print JSON as TXT
                try:
                    cons_txt = json.dumps(json.loads(cons_data), indent=2)
                except Exception:
                    cons_txt = cons_data
                st.download_button(
                    label="Download Conservation Summary (TXT)",
                    data=cons_txt,
                    file_name=f"{selected_receptor}_conservation_summary.txt",
                    mime="text/plain",
                    key=f"download_cons_txt_{selected_receptor}"
                )
            
            # Receptor PDB
            if receptor_files["pdb_receptor"] and receptor_files["pdb_receptor"].exists():
                with open(receptor_files["pdb_receptor"], 'r') as f:
                    receptor_data = f.read()
                st.download_button(
                    label="Download Receptor PDB",
                    data=receptor_data,
                    file_name=f"{selected_receptor}_receptor.pdb",
                    mime="chemical/x-pdb",
                    key=f"download_receptor_pdb_{selected_receptor}"
                )
            
            # Ligand PDB
            if receptor_files["pdb_ligand"] and receptor_files["pdb_ligand"].exists():
                with open(receptor_files["pdb_ligand"], 'r') as f:
                    ligand_data = f.read()
                st.download_button(
                    label="Download Ligand PDB",
                    data=ligand_data,
                    file_name=f"{selected_receptor}_ligand.pdb",
                    mime="chemical/x-pdb",
                    key=f"download_ligand_pdb_{selected_receptor}"
                )
        
        with col2:
            pass

if __name__ == "__main__":
    main() 
