from rdkit import Chem
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from PIL import Image, ImageDraw, ImageFont
import numpy as np
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import rdDepictor
rdDepictor.SetPreferCoordGen(True)

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
    total_width = mol_img1.width + mol_img2.width + 100  # 100 是分子间距
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
    draw.line([arrow_start, arrow_end], fill="black", width=3)
    draw.polygon([(arrow_end[0], arrow_end[1]), 
                  (arrow_end[0] - 6, arrow_end[1] - 6), 
                  (arrow_end[0] - 6, arrow_end[1] + 6)], 
                 fill="black")  # 箭头三角部分

    # 添加文字标注
    try:
        font = ImageFont.truetype("arial.ttf", size=20)
    except:
        font = ImageFont.load_default()

    text_position = (width // 2 - 10, height // 2 - 40)
    draw.text(text_position, label, fill="black", font=font)

    return img


def create_grid(images, rows=5, cols=2, output_file="output.png"):
    """
    将图像排列成网格并保存到文件。
    
    Args:
        images: 包含所有图片的列表 (Pillow Image 对象)
        rows: 行数
        cols: 列数
        output_file: 输出文件名
    """
    if len(images) > rows * cols:
        raise ValueError("图像数量超过网格容量！")

    # 单张图片的尺寸
    img_width, img_height = images[0].size

    # 创建背景
    grid_width = img_width * cols
    grid_height = img_height * rows
    grid_img = Image.new("RGBA", (grid_width, grid_height), (255, 255, 255))

    # 加载字体
    try:
        font = ImageFont.truetype("arial.ttf", size=20)
    except:
        font = ImageFont.load_default(size=25)

    # 将每张图像粘贴到网格上
    for idx, img in enumerate(images):
        row, col = divmod(idx, cols)
        x_offset = col * img_width
        y_offset = row * img_height
        grid_img.paste(img, (x_offset, y_offset))

        # 在图片上方标注序号
        draw = ImageDraw.Draw(grid_img)
        label_position = (x_offset + 10, y_offset + 10)
        draw.text(label_position, f"{idx + 1}", fill="black", font=font)

    # 保存最终的图片
    grid_img.save(output_file)
    print(f"图片已保存为 {output_file}")
if __name__ == "__main__":
    import os
# reactant_smiles = "COc1n[nH]c2ncc(C#Cc3c(F)ccc(CS(=O)(=O)Nc4ccccc4)c3F)cc12"  
    df=pd.read_csv("pred_results.csv")
    for i in range(df.shape[0]):
        reactant_smiles = df['substrate'][i]
        output= df['our_predict'].iloc[i].split("|")
        images = []
        output_dir = "predict_metabolite"
        output_file=f"molecule{i}.png"
        output_path = os.path.join(output_dir, output_file)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for i in range(10):

            product_smiles = output[i]
            # prob = output[i].split(" ")[1]
            mol1 = Chem.MolFromSmiles(reactant_smiles)
            mol2 = Chem.MolFromSmiles(product_smiles)
            image = view_difference_with_arrow(mol1, mol2,label=" ")
            images.append(image)
        # 创建网格并保存
        create_grid(images, rows=5, cols=2, output_file=output_path)