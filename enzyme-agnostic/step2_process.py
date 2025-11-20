
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D




data=pd.read_csv("test_drugs_argu10.csv")

def kekulize_smiles(smi):
    smi = smi.replace(' ','')
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags = True)
    smi = Chem.MolToSmiles(mol, canonical = False, kekuleSmiles = True)
    return smi

data['substrate']=data['substrate'].apply(kekulize_smiles)

data=data['substrate'].to_list()

def space_out_letters(text):
    # 对字符串中的每个字符进行处理，用空格分隔
    # sections = text.split('|')
    # sections[0]
    return ' '.join(list(text))+" | <unk>"

data1=list(map(space_out_letters,data))
result=pd.DataFrame(data1)
result.to_csv("test_drug_input.txt",index=False,header=None)




