# MetaReact

## Installation

 Environment 1

#### Create the environment using the `BBAR.yaml` file:

#### `conda env create -f BBAR.yaml`

Environment 2：rdkit2019 (for data processing related to rdkit and indigo)**

```
conda create -n rdkit2019 python==3.7
conda activate rdkit2019
conda install -c rdkit rdkit=2019.03.2 -y
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple epam.indigo
pip install ipykernel --upgrade
rdkit2019 requires:
\- python <=3.7
```

## model

model.pt, which can be downloaded from https://zenodo.org/records/17898172.

## Use

For each of the following directories — `enzyme-conditioned`, `enzyme-agnostic`, and `enzyme-completion` — please follow these steps:

1. Place your corresponding data and model.pt  into the `test_drugs.csv` file within each folder.
2. Execute the corresponding prediction script:

```
predict_enzyme-conditioned.sh
predict_enzyme-agnostic.sh
predict_enzyme-completion.sh
```

The results are saved in `pred_results.csv`. Additionally, visual diagrams are generated to illustrate the structural transformations between the substrate and its predicted metabolites.

`enzyme-conditioned` is suitable for cases where both the substrate and enzyme are provided, which is referred to as "*enzyme-conditioned*".

`enzyme-agnostic` ignores enzyme information; the input is a substrate, and the output is the corresponding metabolites.

`enzyme-completion` takes the substrate as input and predicts both the possible enzymes and the resulting metabolites.

## Running Details

**Step 1: Data Augmentation**
Apply Test-Time Augmentation (TTA) to the input data.

**Step 2: SMILES Preprocessing**
Kekulize the SMILES strings, then tokenize them into separate tokens for model input.

**Step 3: Model Prediction**
Use the trained model to generate predictions based on the processed tokens.

**Step 4: Output Conversion**
Convert the model's output from ReactSeq format into canonical SMILES representations.
