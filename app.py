import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from streamlit_ketcher import st_ketcher

# Function to calculate mean toxicity, uncertainty, and overall danger
def calculate_toxicity(smiles_list):
    data = []
    for idx, smiles in enumerate(smiles_list):
        mean_toxicity = 10  # Placeholder value
        uncertainty = 1     # Placeholder value
        overall_danger = mean_toxicity + uncertainty
        data.append({
            'Ranking': idx + 1,
            'SMILES': smiles,
            'Mean Toxicity': mean_toxicity,
            'Uncertainty': uncertainty,
            'Overall Danger': overall_danger
        })
    return pd.DataFrame(data)

# Initialize session state for smiles_list
if 'smiles_list' not in st.session_state:
    st.session_state.smiles_list = []


# Title and logo
col1, col2 = st.columns([1, 10])
with col1:
    st.image("figures/GCE_logo.png", use_column_width=True)
with col2:
    st.title("Chemical Toxicity Alert System")


# Dropdown to select input method
input_method = st.selectbox("Select input method", ["SMILES", "Draw Compound"])

# Input method: Draw Compound
if input_method == "Draw Compound":
    DEFAULT_SMILES = "CCO"  # Example default SMILES
    st.write("Draw your compound below:")
    smiles = st_ketcher(DEFAULT_SMILES)
    if st.button("Add compound(s)"):
        st.session_state.smiles_list.append(smiles)
        st.rerun()

# Input method: Comma separated SMILES
elif input_method == "SMILES":
    smiles = st.text_input("SMILES (comma separated)")
    if st.button("Add compound(s)"):
        smiles_list = [smile.strip() for smile in smiles.split(",")]
        st.session_state.smiles_list.extend(smiles_list)
        st.rerun()

# Display the list of compounds
st.subheader("List of Compounds")
if st.session_state.smiles_list:
    for idx, smiles in enumerate(st.session_state.smiles_list):
        st.write(f"{idx + 1}. {smiles}")


# print a line break
st.write("___")

# Display the list of compounds and their toxicity information
if st.button("Run"):
    if st.session_state.smiles_list:
        st.subheader("Ranked List of Compounds")
        df = calculate_toxicity(st.session_state.smiles_list)
        st.dataframe(df)

# Case Study 1
st.subheader("Case Study 1")
st.write("In the Environmental Protection sector of the government we have been able to use the Chemical Toxicity Alert System in order to know the toxicity of both residential and commercial illicit discharge. This allowed us to input a list of chemicals released in these spills and see the toxicity of the chemicals. We then can know if the area needs to be pumped out. Understanding the toxicity of chemicals is one of the most important parts of the job so having software that requires little chemical experience allows us to complete our jobs more effectively.")
st.image("figures/Guy Sampling.jpg")  # Placeholder image URL

# Case Study 2
st.subheader("Case Study 2")
st.write("Lorem ipsum, lorem ipsum")
st.image("https://via.placeholder.com/300")  # Placeholder image URL

# Optionally, visualize the molecules using rdkit
st.subheader("Molecule Visualizations")
for smiles in st.session_state.smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        st.image(Draw.MolToImage(mol), caption=smiles)
