import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from streamlit_ketcher import st_ketcher
import numpy as np
import pickle

# Function to calculate mean toxicity, uncertainty, and overall danger
def calculate_toxicity(smiles_list):
    nBits = [512, 128, 64]  # Different bit lengths for the fingerprints
    filenames = ['XGBRegressor_n512.pkl', 'XGBRegressor_n128.pkl', 'XGBRegressor_n64.pkl']
    models = [pickle.load(open(filename, 'rb')) for filename in filenames]

    all_predictions = []

    for nBit, model in zip(nBits, models):
        fps = []
        # Generate fingerprints for the input SMILES
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=2, nBits=nBit)
                vector = np.array(fp)
                fps.append(vector)
            else:
                fps.append(np.zeros((nBit,), dtype=int))  # Handle invalid SMILES strings

        X_input = np.array(fps)
        y_pred = model.predict(X_input)
        all_predictions.append(y_pred)

    all_predictions = np.array(all_predictions)
    mean_predictions = np.mean(all_predictions, axis=0)
    uncertainties = np.std(all_predictions, axis=0)

    data = []
    for idx, (smiles, mean_pred, uncertainty) in enumerate(zip(smiles_list, mean_predictions, uncertainties)):
        overall_danger = mean_pred + uncertainty
        data.append({
            'Ranking': idx + 1,
            'SMILES': smiles,
            'Mean Toxicity': mean_pred,
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

st.write("The Chemical Toxicity Alert System is a simple method for analyzing a compound or list of compounds. Simply input the chemicals you wish to analyze, and CTAS will not only generate a ranked list of the most toxic compounds, but also a numerical indicator of toxicity.")


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


# print a line break
st.write("___")

# Case Study 1
st.subheader("Case Study: Watershed Testing")
st.write("In the Environmental Protection sector of the government we have been able to use the Chemical Toxicity Alert System in order to know the toxicity of both residential and commercial illicit discharge. This allowed us to input a list of chemicals released in these spills and see the toxicity of the chemicals. We then can know if the area needs to be pumped out. Understanding the toxicity of chemicals is one of the most important parts of the job so having software that requires little chemical experience allows us to complete our jobs more effectively.")
st.image("figures/Guy Sampling.jpg")  # Placeholder image URL

# Case Study 2
st.subheader("Case Study: Reaction By-Products")
st.write("There are many potential cases for the generation of unknown chemicals during an experiment. Running a new reaction, scaling up a reaction, or optimizing the reaction can all result in potentially toxic byproducts that can be difficult to test or product. By separating an impure mixture via LC-MS/GC-MS into fragments, the Chemical Toxicity Alert System (CTAS) can detect the compounds that have the highest potential to be toxic and classify their toxicity from the generated list of impurities, allowing for rapid detection of potentially toxic byproducts.This AI-powered technology simplifies laboratory safety protocols and is a vital tool for any chemist looking to develop greener methods. ")
st.image("figures/mass spec.jpg")  

# Optionally, visualize the molecules using rdkit
st.subheader("Molecule Visualizations")
for smiles in st.session_state.smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        st.image(Draw.MolToImage(mol), caption=smiles)
