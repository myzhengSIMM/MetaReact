import os
import pandas as pd
from SmilesEnumerator import SmilesEnumerator

def augment_smiles(smiles, sme, augmentation_factor):
    augmented_smiles = ["".join(list(smiles))]
    for i in range(augmentation_factor):#固定随机化的种子
        randomized_smiles = sme.randomize_smiles(smiles,i)
        augmented_smiles.append("".join(list(randomized_smiles)))
    return augmented_smiles

def augment_smiles_file(input_file_path, output_file_path, augmentation_factor=9):
    # 创建SmilesEnumerator实例
    sme = SmilesEnumerator()
    
    # 读取原始文件
    try:
        df = pd.read_csv(input_file_path)
    except:
        df=pd.read_csv(input_file_path,encoding="gb18030")
    print(df.columns)
    
    # 检查是否包含'substrate'列
    if 'substrate' not in df.columns:
        raise ValueError("输入文件必须包含'substrate'列")
    
    # 数据增强
    augmented_data = {'substrate': []}
    for _, row in df.iterrows():
        substrate = row['substrate'].strip()

        
        augmented_data['substrate'].extend(augment_smiles(substrate, sme, augmentation_factor))

    
    # 保存增强后的数据
    augmented_df = pd.DataFrame(augmented_data)
    augmented_df.to_csv(output_file_path, index=False)
    
    print(f"数据增强完成，生成的文件名为：{output_file_path}")

# for type in ['train','val','test']:
    # 指定输入和输出文件
input_file_path = "test_drugs.csv"
output_file_path = "test_drugs_argu10.csv"

# 调用函数处理文件
augment_smiles_file(input_file_path, output_file_path)
