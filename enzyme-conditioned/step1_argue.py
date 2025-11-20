import os
import pandas as pd
from SmilesEnumerator import SmilesEnumerator

def augment_smiles(smiles, sme, augmentation_factor):
    augmented_smiles = ["".join(list(smiles))]
    for i in range(augmentation_factor):#固定随机化的种子
        randomized_smiles = sme.randomize_smiles(smiles,8079+i)#42.68.8080.69,41
        augmented_smiles.append("".join(list(randomized_smiles)))
    return augmented_smiles

def augment_smiles_file(input_file_path, output_file_path, augmentation_factor=9):
    # 创建 SmilesEnumerator 实例
    sme = SmilesEnumerator()
    
    # 读取原始文件
    try:
        df = pd.read_csv(input_file_path)
    except:
        df = pd.read_csv(input_file_path, encoding="gb18030")
    print(df.columns)
    
    # 检查是否包含 'substrate' 列
    if 'substrate' not in df.columns:
        raise ValueError("输入文件必须包含 'substrate' 列")
    
    # 数据增强
    augmented_data = {'substrate': [], 'Enzyme': []}
    for _, row in df.iterrows():
        substrate = row['substrate'].strip()
        enzyme = row['Enzyme']  # 保留原始 Enzyme 的值
        
        # 增强 substrate，同时扩展 Enzyme 列
        augmented_smiles = augment_smiles(substrate, sme, augmentation_factor)
        augmented_data['substrate'].extend(augmented_smiles)
        augmented_data['Enzyme'].extend([enzyme] * len(augmented_smiles))  # 保持与 SMILES 对应的 Enzyme 不变
    
    # 保存增强后的数据
    augmented_df = pd.DataFrame(augmented_data)
    augmented_df.to_csv(output_file_path, index=False)
    
    print(f"数据增强完成，生成的文件名为：{output_file_path}")

# 指定输入和输出文件路径
input_file_path = "test_drugs.csv"
output_file_path = "predict_argu10.csv"

# 调用函数处理文件
augment_smiles_file(input_file_path, output_file_path)
