
from rdkit import Chem
from neutral import NeutraliseCharges
from multiprocessing import Pool
from rdkit import RDLogger
import pandas as pd
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import MolStandardize
import os
import pandas as pd
from e_smiles import get_e_smiles, merge_smiles, get_edit_from_e_smiles
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import AllChem

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

if __name__ == "__main__":
    n_best = 10
    substrate_lis = []
    true_tgt_lis = []
    with open('argue_input.txt') as f:
        for line in f.readlines():
            line = line.replace('\n','').replace(" ","").split("|")[0]
            substrate_lis += [line] * n_best
            true_tgt_lis.append(line)

    pred_tgt_lis = []
    with open('output_nbest20.txt') as f:
        for line in f.readlines():
            line = line.replace('\n','').replace(" ","")
            pred_tgt_lis.append(line)

    pred = [f"{x}>>>{y}" for x, y in zip(substrate_lis, pred_tgt_lis)]
    pred_smi= Parallel(n_jobs=-1)(delayed(process_smiles)(x) for x in pred)

    def remove_leavegroup(smi):
        frag = smi.split(".")
        result = [i for i in frag if 'Au' not in i]
        return ".".join(result)

    predict_smiles_all=list(map(remove_leavegroup,pred_smi))

    predict_smiles=list(map(canonicalize_predict,predict_smiles_all))


    predict_cocan_smiles=list(map(canonicalize_smiles,predict_smiles))

    test=pd.read_csv("test_drugs.csv")

    # truth=test['truth']


    fold = 10
    num = int(len(true_tgt_lis) / 10)
    n = n_best #nbest

    score_1 = []
    for i in range(fold):
        for j in range(num):
            for k in range(1,n+1):
                score_1.append(1/k**2)       # square  

    vote_score_lis = [-1]*len(predict_cocan_smiles) 

    id=[]

    for j in range(num):
        smiles_score_dic = {}
        for i in range(fold):
            for k in range(n):
                idx = j * fold * n + i*n + k
                id.append(idx)
                smiles = predict_cocan_smiles[idx]
                smiles = smiles
                s = score_1[idx]
                if smiles == '':
                    s = 0
                
                if smiles not in smiles_score_dic:
                    smiles_score_dic[smiles] = s
                else:
                    smiles_score_dic[smiles] += s

        for i in range(fold):
            for k in range(n):
                idx = j * fold * n + i*n + k
                smiles = predict_cocan_smiles[idx]
                vote_score = smiles_score_dic[smiles]
                vote_score_lis[j * fold * n + i*n + k] = vote_score

    pre_smiles_50_lis = []
    pre_reactseq_50_lis = []

    for j in range(num):
        smiles_g = []
        vote_score_g = []
        for i in range(fold):
            for k in range(n):
                idx = j * fold * n + i*n + k
                smiles_g.append(predict_cocan_smiles[idx])
                vote_score_g.append(vote_score_lis[idx])

        zip_a_b_c = zip(smiles_g,  vote_score_g)
        sorted_zip = sorted(zip_a_b_c, key=lambda x:x[-1],reverse= True)
        smiles_g,  vote_score_g = zip(*sorted_zip)
        
        pre_smiles_50 = []
        # pre_reactseq_50 = []
        for i in smiles_g:
            if i not in pre_smiles_50 and len(pre_smiles_50) < 50:
                pre_smiles_50.append(i)
                # pre_reactseq_50.append(j)

        if len(pre_smiles_50) < 50:
            pre_smiles_50 = pre_smiles_50 + [''] * (50-len(pre_smiles_50))
            # pre_reactseq_50 = pre_reactseq_50 + [pre_reactseq_50[-1].split('>>>')[0]+'>>>'] * (50-len(pre_reactseq_50))
            
        pre_smiles_50_lis.append(pre_smiles_50)
        # pre_reactseq_50_lis.append(pre_reactseq_50)

    def get_similarity(smiles1, smiles2):
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        if mol1 and mol2:
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
            return TanimotoSimilarity(fp1, fp2)

    # num=129
    all=[]
    for i in range(num):
        data=[]
        for k in range(n_best):  #这里   
            data.append(pre_smiles_50_lis[i][k])
        all.append("|".join(data))

    data=pd.read_csv("test_drugs.csv")

    data['our_predict']=all

    df=data.copy()
    df.to_csv("pred_results.csv", index=False)

