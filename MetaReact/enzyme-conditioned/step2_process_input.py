

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

def kekulize_smiles(smi):
    smi = smi.replace(' ','')
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags = True)
    smi = Chem.MolToSmiles(mol, canonical = False, kekuleSmiles = True)
    return smi

data_all= pd.read_csv("predict_argu10.csv")

data=data_all['substrate'].to_list()


data1=list(map(kekulize_smiles,data))

def space_out_letters(text):
    
    return ' '.join(list(text)) + " | "
input_smi = list(map(space_out_letters, data1)) 
enzyme_list = data_all['Enzyme'].tolist()
inputs = [f"{smi}{enzyme}" for smi, enzyme in zip(input_smi, enzyme_list)]


with open("argue_input.txt", "w", encoding="utf-8") as file:
    for line in inputs:
        file.write(line + "\n") 




