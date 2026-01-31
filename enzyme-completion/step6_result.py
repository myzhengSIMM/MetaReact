
from rdkit import Chem
from neutral import NeutraliseCharges
from multiprocessing import Pool
from rdkit import RDLogger
import pandas as pd
import json
from collections import defaultdict
import pandas as pd
from e_smiles import get_e_smiles, merge_smiles, get_edit_from_e_smiles
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')

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
    with open('test_enzyme_argue10_input_group.txt') as f:
        for line in f.readlines():
            line = line.replace('\n','').replace(" ","")
            substrate_lis += [line] * n_best
            true_tgt_lis.append(line)

    data=pd.read_csv("test_enzyme_argue10_add_smiles.tsv")

    n_best=10
    enzymes = data['enzyme'].fillna('').astype(str).tolist()
    predict_cocan_smiles = data['smiles'].fillna('').tolist()

    fold = 10
    num = int(data.shape[0]/ 100)
    n = n_best  
    enzyme_prefixes = ['CYP', 'UGT', 'SULT', 'AKR', 'GST', 'CES', 'FMO', 'ALDH', 'NAT', 'COMT', 'AOX', 'XDH', 'EPHX', 'NQO', 'ADH']

    # 计算分数
    score_1 = []
    for i in range(fold):
        for j in range(num):
            for k in range(1, n + 1):
                score = 1 / k**2
                score_1.append(score)

    # 初始化
    vote_score_lis = [-1] * len(enzymes)
    id_list = []

    # 计算 vote_score_lis
    for j in range(num):
        enzyme_score_dic = {}
        enzyme_category_map = {}

        for i in range(fold):
            for k in range(n):
                idx = j * fold * n + i * n + k
                id_list.append(idx)
                
                enzyme = enzymes[idx]
                s = score_1[idx] if enzyme != '' else 0
                
                if enzyme == '':
                    continue  # 跳过空值

                # 根据 enzyme 的前缀分类
                enzyme_category = None
                for prefix in enzyme_prefixes:
                    if enzyme.startswith(prefix):
                        enzyme_category = prefix
                        break

                if enzyme_category is None:
                    continue  # 如果没有匹配到大类，跳过
                
                # 大类和亚型分别累积打分
                key = enzyme_category
                subtype_key = enzyme

                # 大类累积打分
                if key not in enzyme_score_dic:
                    enzyme_score_dic[key] = s
                else:
                    enzyme_score_dic[key] += s

                # 亚型累积打分
                if subtype_key not in enzyme_category_map:
                    enzyme_category_map[subtype_key] = s
                else:
                    enzyme_category_map[subtype_key] += s

        # 更新 vote_score_lis
        for i in range(fold):
            for k in range(n):
                idx = j * fold * n + i * n + k
                enzyme = enzymes[idx]

                # 根据大类和亚型分开处理
                enzyme_category = None
                for prefix in enzyme_prefixes:
                    if enzyme.startswith(prefix):
                        enzyme_category = prefix
                        break

                if enzyme_category is None:
                    continue  # 如果没有匹配到大类，跳过

                # 更新总分
                vote_score_lis[idx] = enzyme_category_map.get(enzyme, -1)

    # 创建最终的 pre_enzyme_category_50_lis 和 pre_enzyme_50_lis
    pre_enzyme_category_50_lis = []
    pre_enzyme_50_lis = []
    normalized_values=[]
    for j in range(num):
        enzyme_g = []
        vote_score_g = []
        
        for i in range(fold):
            for k in range(n):
                idx = j * fold * n + i * n + k
                enzyme_g.append(enzymes[idx])
                vote_score_g.append(vote_score_lis[idx])
                
        zip_a_b = zip(enzyme_g, vote_score_g)
        sorted_zip = sorted(zip_a_b, key=lambda x: x[1], reverse=True)
        enzyme_g, vote_score_g = zip(*sorted_zip)

        aggregated_values = defaultdict(list)
        for key, value in zip(enzyme_g, vote_score_g):
            aggregated_values[key].append(value)


        average_values = {
        key: (sum(float(value) for value in values) / len(values) if values else 0)
        for key, values in aggregated_values.items()
        }
        # 计算 min 和 max
        min_value = min(average_values.values())
        max_value = max(average_values.values())

        # 避免 max_value 和 min_value 相等导致除以零的情况
        if max_value == min_value:
            normalized_values[j] = {key: 0.0 for key in average_values}  # 所有值都设置为 0
        else:
            normalized_values.append( {
                key: (value - min_value) / (max_value - min_value)*0.9
                for key, value in average_values.items()
            })

            
        pre_enzyme_category_50 = []
        pre_enzyme_50 = []
        
        seen_enzyme_category = set()
        seen_enzyme = set()
        
        for enzyme in enzyme_g:
            enzyme_category = None
            for prefix in enzyme_prefixes:
                if enzyme.startswith(prefix):
                    enzyme_category = prefix
                    break

            if enzyme_category is None:
                continue  # 如果没有匹配到大类，跳过

            # 先检查是否已经有该大类的条目
            if enzyme_category in seen_enzyme_category:
                # 如果大类已经存在，则仅保留带有亚型的酶
                if enzyme != enzyme_category and enzyme not in seen_enzyme and len(pre_enzyme_50) < 50:
                    pre_enzyme_50.append(enzyme)
                    seen_enzyme.add(enzyme)
            else:
                # 大类处理
                if len(pre_enzyme_category_50) < 50:
                    pre_enzyme_category_50.append(enzyme_category)
                    seen_enzyme_category.add(enzyme_category)

                # 亚型处理
                if enzyme != enzyme_category and len(pre_enzyme_50) < 50:
                    pre_enzyme_50.append(enzyme)
                    seen_enzyme.add(enzyme)
        
        pre_enzyme_category_50.extend([''] * (50 - len(pre_enzyme_category_50)))
        pre_enzyme_50.extend([''] * (50 - len(pre_enzyme_50)))

        pre_enzyme_category_50_lis.append(pre_enzyme_category_50)
        pre_enzyme_50_lis.append(pre_enzyme_50)


    combined_list = []

    # 合并 pre_enzyme_category_50_lis 和 pre_enzyme_50_lis，只保留非空元素
    for category_list, enzyme_list in zip(pre_enzyme_category_50_lis, pre_enzyme_50_lis):
        # 过滤掉空字符串
        filtered_category = [item for item in category_list if item != '']
        filtered_enzyme = [item for item in enzyme_list if item != '']
        
        # 合并非空元素
        combined = filtered_category + filtered_enzyme
        combined_list.append(combined)


    # 保存为 JSON 文件
    with open('combined_enzyme.json', 'w', encoding='utf-8') as f:
        json.dump(combined_list, f, ensure_ascii=False, indent=4)


    # Load JSON data from files
    def load_json_data(filepath):
        with open(filepath, "r") as file:
            return json.load(file)

    # Save data to a new JSON file
    def save_json_data(data, filepath):
        with open(filepath, "w") as file:
            json.dump(data, file)

    # Function to find a metabolite matching by enzyme name or prefix
    def find_metabolite_by_prefix(enzyme, metabolite_list):#, prefixes):
        matches = []
        flag=False
        # Try to find an exact match first
        for metabolite in metabolite_list:
            if enzyme == metabolite[0]:
                matches.append(canonicalize_smiles(metabolite[1]))
                flag=True
        if flag==True:
            return '|'.join(matches) if matches else None
        # If no exact match, try to match by known prefixes

            # Check if any metabolite starts with the same prefix
        for metabolite in metabolite_list:
            if metabolite[0].startswith(enzyme):
                matches.append(canonicalize_smiles(metabolite[1]))
        return '|'.join(list(set(matches))) if matches else None  # Return None if no match found

    # Enzyme prefixes for matching
    enzyme_prefixes = ['CYP', 'UGT', 'SULT', 'AKR', 'GST', 'CES', 'FMO', 'ALDH', 'NAT', 'COMT', 'AOX', 'XDH', 'EPHX', 'NQO', 'ADH']

    # Load both datasets
    enzymes_data = load_json_data("combined_enzyme.json")
    metabolites_data = load_json_data("combined_list.json")

    # Check if the lists are of the correct length and structure
    if len(enzymes_data) == num and len(metabolites_data) == num:
        combined_data = []
        for enzyme_list, metabolite_list,normalized_value in zip(enzymes_data, metabolites_data,normalized_values):
            combined_group = []
            for enzyme in enzyme_list:
                # Find matching metabolite by enzyme or prefix
                metabolite = find_metabolite_by_prefix(enzyme, metabolite_list)#, enzyme_prefixes)
                if enzyme not in enzyme_prefixes:
                    score_enzyme =normalized_value[enzyme]
                    if score_enzyme<0.1:
                        score_enzyme=0.3
                    combined_group.append((enzyme, metabolite,score_enzyme))
                else:
                    combined_group.append((enzyme, metabolite,1))
            combined_data.append(combined_group)

        save_json_data(combined_data, "path_to_updated_combined_data.json")
    else:
        print("Error: Data lists do not match the required structure or length.")


    for i in range(len(combined_list)):
        combined_list[i] = combined_list[i] + [""] * (20 - len(combined_list[i]))

    all=[]
    for i in range(num):
        data=[]
        for k in range(20):     
            data.append(combined_list[i][k])
        all.append("|".join(data))

    def kekulize_smiles(smi):
        smi = smi.replace(' ','')
        mol = Chem.MolFromSmiles(smi)
        Chem.Kekulize(mol, clearAromaticFlags = True)
        smi = Chem.MolToSmiles(mol, canonical = False, kekuleSmiles = True)
        return smi

    df=pd.read_csv("test_drugs.csv")

    formatted_results = []
    for data in combined_data:
        formatted_data = "<separated>".join([f"{enzyme}|{metabolite}|{score:.2f}" for enzyme, metabolite, score in data])
        formatted_results.append(formatted_data)

    df['pred']=formatted_results
    df['substrate']=df['substrate'].apply(canonicalize_smiles)
    df.to_csv("pred_results.csv",index=False)

    # grouped_data = defaultdict(list)

    # for info in combined_data:
    #     enzyme, metabolite, score = info
    #     grouped_data[metabolite].append((enzyme, score))

    # # 计算每个 metabolite 对应的酶和这些酶的平均分数
    # formatted_results = []

    # for metabolite, enzymes_scores in grouped_data.items():
    #     enzymes = [enzyme for enzyme, score in enzymes_scores]
    #     avg_score = sum(score for enzyme, score in enzymes_scores) / len(enzymes_scores)
    #     formatted_data = f"Metabolite: {metabolite}, Enzymes: {', '.join(enzymes)}, Average Score: {avg_score:.2f}"
    #     formatted_results.append(formatted_data)

    # # 输出结果
    # for result in formatted_results:
    #     print(result)


