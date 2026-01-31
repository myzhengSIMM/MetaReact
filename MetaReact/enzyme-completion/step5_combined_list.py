
import pandas as pd
import json
data=pd.read_csv("test_enzyme_argue10_add_smiles.tsv")


data['pred']=data['enzyme']+"|"+data['smiles']
n_best=10
enzymes = data['pred'].str.split('|').str[0].fillna('').astype(str).tolist()
data['pred'].str.split('|').str[1].fillna('').astype(str).tolist()
enzymes = data['pred'].str.split('|').str[0].fillna('').astype(str).tolist()
metabolite = data['pred'].str.split('|').str[1].fillna('').astype(str).tolist()

fold = 10
num = int(data.shape[0]/ 100)
n = n_best  
enzyme_prefixes = ['CYP', 'UGT', 'SULT', 'AKR', 'GST', 'CES', 'FMO', 'ALDH', 'NAT', 'COMT', 'AOX', 'XDH', 'EPHX', 'NQO', 'ADH']

score_1 = []
for i in range(fold):
    for j in range(num):
        for k in range(1, n + 1):
            score = 1 / k**2
            score_1.append(score)

vote_score_lis = [-1] * len(enzymes)
id_list = []

# 用于记录每个 enzyme 对应的 metabolite
enzyme_metabolite_map = [-1] * len(enzymes)

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
                    break  # 找到匹配的前缀后立即退出循环
            
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

    # 更新 vote_score_lis 和 enzyme_metabolite_map
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
            enzyme_metabolite_map[idx] = metabolite[idx]  # 记录对应的 metabolite


pre_enzyme_category_50_lis = []
pre_enzyme_50_lis = []
pre_metabolite_50_lis = []  

for j in range(num):
    enzyme_g = []
    vote_score_g = []
    metabolite_g = []  

    for i in range(fold):
        for k in range(n):
            idx = j * fold * n + i * n + k
            enzyme_g.append(enzymes[idx])
            vote_score_g.append(vote_score_lis[idx])
            metabolite_g.append(metabolite[idx])  # 同时记录 metabolite
    
    # zip 和排序，同时包括 metabolite
    zip_a_b_c = zip(enzyme_g, vote_score_g, metabolite_g)
    sorted_zip = sorted(zip_a_b_c, key=lambda x: x[1], reverse=True)
    enzyme_g, vote_score_g, metabolite_g = zip(*sorted_zip)
    
    pre_enzyme_category_50 = []
    pre_enzyme_50 = []
    pre_metabolite_50 = []  # 新增用于存储 metabolite
    
    seen_enzyme_category = set()
    seen_enzyme = set()
    
    for enzyme, metabolite_val in zip(enzyme_g, metabolite_g):  # 包括 metabolite
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
                pre_metabolite_50.append(metabolite_val)  # 添加对应的 metabolite
                seen_enzyme.add(enzyme)
        else:
            # 大类处理
            if len(pre_enzyme_category_50) < 50:
                pre_enzyme_category_50.append(enzyme_category)
                seen_enzyme_category.add(enzyme_category)

            # 亚型处理
            if enzyme != enzyme_category and len(pre_enzyme_50) < 50:
                pre_enzyme_50.append(enzyme)
                pre_metabolite_50.append(metabolite_val)  # 添加对应的 metabolite
                seen_enzyme.add(enzyme)
    
    # 填充空位
    pre_enzyme_category_50.extend([''] * (50 - len(pre_enzyme_category_50)))
    pre_enzyme_50.extend([''] * (50 - len(pre_enzyme_50)))
    pre_metabolite_50.extend([''] * (50 - len(pre_metabolite_50)))  # 填充 metabolite 列表的空位

    # 添加到最终的列表中
    pre_enzyme_category_50_lis.append(pre_enzyme_category_50)
    pre_enzyme_50_lis.append(pre_enzyme_50)
    pre_metabolite_50_lis.append(pre_metabolite_50)  # 记录对应的 metabolite 列表


combined_list = []

# 合并 pre_enzyme_category_50_lis, pre_enzyme_50_lis 和 pre_metabolite_50_lis
for category_list, enzyme_list, metabolite_list in zip(pre_enzyme_category_50_lis, pre_enzyme_50_lis, pre_metabolite_50_lis):
    combined = []
    
    # 遍历所有 enzyme 和 metabolite，确保一一对应
    for enzyme, metabolite in zip(enzyme_list, metabolite_list):
        # 过滤掉空的 enzyme 和 metabolite
        if enzyme != '' and metabolite != '':
            # 根据 enzyme 的前缀找到对应的 category
            enzyme_category = None
            for prefix in enzyme_prefixes:
                if enzyme.startswith(prefix):
                    enzyme_category = prefix
                    break
            
            # 如果找到了 category，合并结果
            if enzyme_category:
                combined.append((enzyme, metabolite))  # 合并 category, enzyme 和 metabolite

    combined_list.append(combined)

# for combined in combined_list:
#     print(combined)

with open('combined_list.json', 'w', encoding='utf-8') as f:
    json.dump(combined_list, f, ensure_ascii=False, indent=4)




