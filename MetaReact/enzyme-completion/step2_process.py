
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

data=pd.read_csv("test_argue10_group_new.tsv")
def kekulize_smiles(smi):
    smi = smi.replace(' ','')
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags = True)
    smi = Chem.MolToSmiles(mol, canonical = False, kekuleSmiles = True)
    return smi


data['substrate_smiles']=data['substrate'].apply(kekulize_smiles)

data=data['substrate_smiles'].to_list()

def space_out_letters(text):
    return ' '.join(list(text))
data1=list(map(space_out_letters,data))
result=pd.DataFrame(data1)
result.to_csv("test_enzyme_argue10_input_group.txt",index=False,header=None)



