
from rdkit import Chem
from neutral import NeutraliseCharges
from multiprocessing import Pool
from rdkit import RDLogger
import pandas as pd
from rdkit.Chem import MolStandardize
RDLogger.DisableLog('rdApp.*')
import os
import pandas as pd
from e_smiles import get_e_smiles, merge_smiles, get_edit_from_e_smiles
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import AllChem
import csv
def canonicalize_smiles(smiles):
    if len(smiles)==0:
        return ''
    mol = Chem.MolFromSmiles(smiles)
    # lfc = MolStandardize.fragment.LargestFragmentChooser()
    
    if mol is not None:
        # mol2 = lfc.choose(mol)
        smi2=Chem.MolToSmiles(mol, isomericSmiles=True)
        smi,_=NeutraliseCharges(smi2)
        return smi 
    else:
        return ''


def canonicalize_predict(smiles):
    if len(smiles)==0:
        return ''
    mol = Chem.MolFromSmiles(smiles)
    lfc = MolStandardize.fragment.LargestFragmentChooser()
    
    if mol is not None:
        mol2 = lfc.choose(mol)
        smi2=Chem.MolToSmiles(mol2, isomericSmiles=True)
        smi,_=NeutraliseCharges(smi2)
        return smi
    else:
        return ''


def process_smiles(smiles):
    return merge_smiles(smiles)

def kekulize_smiles(smi):
    smi = smi.replace(' ','')
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags = True)
    smi = Chem.MolToSmiles(mol, canonical = False, kekuleSmiles = True)
    return smi

n_best = 10
substrate_lis = []
true_tgt_lis = []
with open('test_enzyme_argue10_input_group.txt') as f:
    for line in f.readlines():
        line = line.replace('\n','').replace(" ","")
        substrate_lis += [line] * n_best
        true_tgt_lis.append(line)


import csv

# 输入文件路径
input_file_path = 'test_enzyme_nbest10_argue10.txt'
# 输出文件路径，使用.tsv后缀
output_file_path = 'test_enzyme_argue10.tsv'

# 打开输出的TSV文件
with open(output_file_path, 'w', newline='') as tsvfile:
    # 定义列名
    fieldnames = ['enzyme', 'metabolite']
    # 创建CSV写入器，使用\t作为分隔符
    writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
    
    # 写入列名
    writer.writeheader()

    # 打开并读取输入文件
    with open(input_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if '|' in line:
                # 分割酶名称和化学结构
                enzyme, metabolite = line.split('|')
                enzyme = enzyme.strip()
                metabolite = metabolite.replace(' ', '').strip()  # 去掉metabolite中的所有空格
            else:
                # 如果没有 '|'，酶的名称设为 None，metabolite 保持不变
                enzyme = None
                metabolite = line.replace(' ', '').strip()  # 也去掉空格

            # 写入一行到TSV文件
            writer.writerow({'enzyme': enzyme, 'metabolite': metabolite})

print(f"Data has been written to {output_file_path}")



pred_tgt_lis = []
with open('test_enzyme_argue10.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        metabolite = row['metabolite'].replace(" ", "")  # 去除空格
        pred_tgt_lis.append(metabolite)


pred = [f"{x}>>>{y}" for x, y in zip(substrate_lis, pred_tgt_lis)]
pred_smi= Parallel(n_jobs=-1)(delayed(process_smiles)(x) for x in pred)


def remove_leavegroup(smi):
    frag = smi.split(".")
    result = [i for i in frag if 'Au' not in i]
    return ".".join(result)

predict_smiles_all=list(map(remove_leavegroup,pred_smi))
predict_cocan_smiles=list(map(canonicalize_smiles,predict_smiles_all))
data=pd.read_csv("test_enzyme_argue10.tsv",sep='\t')
data['smiles']=predict_cocan_smiles
data.to_csv("test_enzyme_argue10_add_smiles.tsv",index=False)
