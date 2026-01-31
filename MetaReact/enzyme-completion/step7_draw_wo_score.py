
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdDepictor
rdDepictor.SetPreferCoordGen(True)
IPythonConsole.drawOptions.minFontSize=20
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdDepictor
rdDepictor.SetPreferCoordGen(True)
IPythonConsole.drawOptions.minFontSize=20
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from PIL import Image, ImageDraw, ImageFont
from IPython.display import Image as IPImage




def extract_enzyme_and_mol2(pred_str):
    # 分隔字符串，首先按照 <separated> 分开
    parts = pred_str.split('<separated>')[1:]

    # 分割后第一个部分是 enzyme, 第二部分是 mol2 信息
    enzyme_mol2_list = []
    for part in parts:
        # 按照竖线分割 enzyme 和 mol2
        if '|' in part:
            enzyme, mol2 = part.split('|')
            enzyme_mol2_list.append((enzyme, mol2))
        else:
            # 如果没有竖线，直接视为没有 mol2
            enzyme_mol2_list.append((part, None))

    return enzyme_mol2_list




def extract_enzyme_and_mol2(pred_str):
    # 分隔字符串，首先按照 <separated> 分开
    parts = pred_str.split('<separated>')[1:]  # 从第二个部分开始

    # 使用字典存储 mol2 及其对应的 enzyme 列表
    mol2_dict = {}

    # 分割后第一个部分是 enzyme, 第二部分是 mol2 信息
    for part in parts:
        # 按照竖线分割 enzyme 和 mol2
        if '|' in part:
            enzyme, mol2 = part.split('|')
            if mol2 not in mol2_dict:
                mol2_dict[mol2] = []
            mol2_dict[mol2].append(enzyme)

    # 将相同 mol2 的 enzyme 合并为逗号分隔的字符串
    mol2_enzyme_list = [(mol2, ', '.join(enzyme_list)) for mol2, enzyme_list in mol2_dict.items()]

    return mol2_enzyme_list



def kekulize_smiles(smi):
    smi = smi.replace(' ','')
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags = True)
    smi = Chem.MolToSmiles(mol, canonical = False, kekuleSmiles = True)
    return smi


def view_difference_with_arrow(mol1, mol2, label="CYP"):
    """
    比较两个分子的不同之处，并在中间添加箭头和标注。
    
    Args:
        mol1: 第一个分子 (RDKit Mol 对象)
        mol2: 第二个分子 (RDKit Mol 对象)
        label: 箭头上的文本（默认 "CYP"）
    
    Returns:
        Image: 带箭头和标注的图片
    """
    # 找到两个分子的最大公共子结构 (MCS)
    mcs = rdFMCS.FindMCS([mol1, mol2])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

    # 获取 mol1 和 mol2 中未匹配的原子索引
    match1 = mol1.GetSubstructMatch(mcs_mol)
    target_atm1 = [atom.GetIdx() for atom in mol1.GetAtoms() if atom.GetIdx() not in match1]

    match2 = mol2.GetSubstructMatch(mcs_mol)
    target_atm2 = [atom.GetIdx() for atom in mol2.GetAtoms() if atom.GetIdx() not in match2]

    # 创建分子对比图
    mol_img1 = Draw.MolToImage(mol1, highlightAtoms=target_atm1, size=(300, 300))
    mol_img2 = Draw.MolToImage(mol2, highlightAtoms=target_atm2, size=(300, 300))

    # 创建一个空白的图像，作为背景，宽度是两张图像的总和，调整分子间的距离
    total_width = mol_img1.width + mol_img2.width + 100  # 100 是分子间距，可以根据需要调整
    max_height = max(mol_img1.height, mol_img2.height)

    # 创建一个白色背景
    background = Image.new("RGBA", (total_width, max_height), (255, 255, 255))

    # 将两张分子图像粘贴到背景中
    background.paste(mol_img1, (0, 0))
    background.paste(mol_img2, (mol_img1.width + 100, 0))  # 在分子之间增加间距

    # 转换为 RGBA 图像
    img = background.convert("RGBA")
    width, height = img.size

    # 创建绘图对象
    draw = ImageDraw.Draw(img)

    # 在分子间绘制箭头
    arrow_start = (mol_img1.width + 10, height // 2)
    arrow_end = (mol_img1.width + 50 + 50, height // 2)
    draw.line([arrow_start, arrow_end], fill="black", width=3)  # 减小线条宽度
    draw.polygon([(arrow_end[0], arrow_end[1]), 
                  (arrow_end[0] - 6, arrow_end[1] - 6),  # 减小三角形的大小
                  (arrow_end[0] - 6, arrow_end[1] + 6)], 
                 fill="black")  # 箭头三角部分

    # 添加文字标注
    try:
        # 加载默认字体
        font = ImageFont.truetype("arial.ttf", size=30)
    except:
        # 如果加载字体失败，则使用默认字体
        font = ImageFont.load_default(size=15)

    # 使用 textbbox 来计算文本的边界框
    bbox = draw.textbbox((0, 0), label, font=font)
    text_width = bbox[2] - bbox[0]  # 计算文本宽度
    text_height = bbox[3] - bbox[1]  # 计算文本高度

    # 计算文本位置
    text_position = (width // 2 - text_width // 2, height // 2 - 40)
    draw.text(text_position, label, fill="black", font=font)

    return img


if __name__ == "__main__":

    data=pd.read_csv("my_predict_enzyme_metabolite.csv")
    substrate=pd.read_csv("test_drugs.csv")['substrate'].tolist()
    origin=pd.read_csv("test_drugs.csv")
    output = data['pred'].apply(extract_enzyme_and_mol2)
    import os

    output_dir = "enzyme_plots"
    os.makedirs(output_dir, exist_ok=True)

    for i in range(len(substrate)):
        for j in range(len(output[i])):
            reactant_smiles = substrate[i]
            product_smiles = output[i][j][0]
            enzyme = output[i][j][1]
            mol1 = Chem.MolFromSmiles(reactant_smiles)
            mol2 = Chem.MolFromSmiles(product_smiles)
            image = view_difference_with_arrow(mol1, mol2, label=enzyme)
            output_file = f"{origin['name'][i]}.png"
            path=os.path.join(output_dir,output_file)
            os.makedirs(path, exist_ok=True)
            file_name = os.path.join(path,f"plot_enzyme_{j+1}.png")
            image.save(file_name)
            print(f"Saved: {file_name}")


