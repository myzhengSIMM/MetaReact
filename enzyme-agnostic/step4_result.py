from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from collections import defaultdict, Counter
from multiprocessing import Pool
from joblib import Parallel, delayed
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import AllChem
from neutral import NeutraliseCharges
from e_smiles import get_e_smiles, merge_smiles, get_edit_from_e_smiles
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from PIL import Image, ImageDraw, ImageFont
import numpy as np
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import rdDepictor
rdDepictor.SetPreferCoordGen(True)
IPythonConsole.drawOptions.minFontSize=20
import pandas as pd

# 禁用 RDKit 日志
RDLogger.DisableLog('rdApp.*')

# ==================== Functions ====================

def canonicalize_smiles(smiles):
    """标准化 SMILES 表示"""
    if len(smiles) == 0:
        return ''
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        smi2 = Chem.MolToSmiles(mol, isomericSmiles=True)
        smi, _ = NeutraliseCharges(smi2)
        return smi
    else:
        return ''


def canonicalize_predict(smiles):
    """标准化 SMILES 表示并选择最大片段"""
    if len(smiles) == 0:
        return ''
    mol = Chem.MolFromSmiles(smiles)
    lfc = rdMolStandardize.fragment.LargestFragmentChooser()
    if mol is not None:
        mol2 = lfc.choose(mol)
        smi2 = Chem.MolToSmiles(mol2, isomericSmiles=True)
        smi, _ = NeutraliseCharges(smi2)
        return smi
    else:
        return ''


def process_smiles(smiles):
    """处理 SMILES"""
    return merge_smiles(smiles)


def kekulize_smiles(smi):
    """将 SMILES 转为 Kekulé 形式"""
    smi = smi.replace(' ', '')
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    smi = Chem.MolToSmiles(mol, canonical=False, kekuleSmiles=True)
    return smi


def remove_leavegroup(smi):
    """去除 SMILES 中的离去基团"""
    frag = smi.split(".")
    result = [i for i in frag if 'Au' not in i]
    return ".".join(result)


def get_similarity(smiles1, smiles2):
    """计算两个 SMILES 的 Tanimoto 相似度"""
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 and mol2:
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        return TanimotoSimilarity(fp1, fp2)
    return 0

# ==================== Data Preparation ====================

n_best = 20
substrate_lis = []
true_tgt_lis = []

# 读取底物文件
with open('test_drug_input.txt') as f:
    for line in f.readlines():
        line = line.replace('\n', '').replace(" ", "").split("|")[0]
        substrate_lis += [line] * n_best
        true_tgt_lis.append(line)

# 读取预测文件
pred_tgt_lis = []
with open('test_drug_nbest20.txt') as f:
    for line in f.readlines():
        line = line.replace('\n', '').replace(" ", "")
        pred_tgt_lis.append(line)

# 构建反应格式
pred = [f"{x}>>>{y}" for x, y in zip(substrate_lis, pred_tgt_lis)]
pred_smi = Parallel(n_jobs=-1)(delayed(process_smiles)(x) for x in pred)

# 去除离去基团
predict_smiles_all = list(map(remove_leavegroup, pred_smi))

# ==================== Canonicalization and Deduplication ====================

# 记录每个 canonical_smiles 对应的所有原始 smiles 的频率
smiles_frequency = defaultdict(Counter)
for smiles in predict_smiles_all:
    canonical_smiles = canonicalize_predict(smiles)
    if canonical_smiles:
        smiles_frequency[canonical_smiles][smiles] += 1

# 构建结果字典，只保留出现频率最高的原始 smiles
predict_smiles_dict = {canonical_smiles: max(counter, key=counter.get)
                       for canonical_smiles, counter in smiles_frequency.items()}

# 标准化
predict_smiles = list(map(canonicalize_predict, predict_smiles_all))
predict_cocan_smiles = list(map(canonicalize_smiles, predict_smiles))

# ==================== Top-K Accuracy ====================

fold = 10
n = n_best  # nbest
score_1 = [1/k**2 for _ in range(fold) for k in range(1, n+1)]

vote_score_lis = [-1] * len(predict_cocan_smiles)
num = int(len(true_tgt_lis) / 10)

# 计算投票得分
for j in range(num):
    smiles_score_dic = {}
    for i in range(fold):
        for k in range(n):
            idx = j * fold * n + i*n + k
            smiles = predict_cocan_smiles[idx]
            s = score_1[idx] if smiles != '' else 0
            smiles_score_dic[smiles] = smiles_score_dic.get(smiles, 0) + s

    for i in range(fold):
        for k in range(n):
            idx = j * fold * n + i*n + k
            smiles = predict_cocan_smiles[idx]
            vote_score_lis[idx] = smiles_score_dic[smiles]

# ==================== Generate Top-50 Predictions ====================

pre_smiles_50_lis = []
for j in range(num):
    smiles_g, vote_score_g = [], []
    for i in range(fold):
        for k in range(n):
            idx = j * fold * n + i*n + k
            smiles_g.append(predict_cocan_smiles[idx])
            vote_score_g.append(vote_score_lis[idx])

    sorted_zip = sorted(zip(smiles_g, vote_score_g), key=lambda x: x[-1], reverse=True)
    smiles_g, _ = zip(*sorted_zip)
    
    pre_smiles_50 = []
    for i in smiles_g:
        if i not in pre_smiles_50 and len(pre_smiles_50) < 50:
            pre_smiles_50.append(i)
    pre_smiles_50 += [''] * (50 - len(pre_smiles_50))  # 补齐长度
    pre_smiles_50_lis.append(pre_smiles_50)

# ==================== Output Results ====================

all = []
for i in range(num):
    data = [canonical_smiles(kekulize_smiles(pre_smiles_50_lis[i][k])) for k in range(20)]
    all.append("|".join(data))

# 保存结果
df = pd.read_csv("test_drugs.csv", encoding='gb18030')
df['predict'] = all
df.to_csv("pred_results.csv", index=False)
