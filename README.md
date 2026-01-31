# MetaReact

MetaReact is a deep learning framework for **metabolic reaction prediction**, designed to flexibly support multiple prediction scenarios with varying levels of enzyme information.  
It enables the prediction of metabolite structures, enzyme–substrate relationships, and complete metabolic pathways in a unified modeling framework.

---

## ⚙️ Installation

To set up the project, you need to configure two separate Conda environments: one for the core model and one for RDKit-based data processing.

### 🐍 Environment 1: Core Environment
This environment is used for the main model inference. Create it using the provided YAML file:

```bash
conda env create -f BBAR.yaml
```

Activate the environment before running predictions:

```bash
conda activate BBAR
```

### 🧪 Environment 2: rdkit2019

(for data processing related to RDKit and Indigo)

This environment is required for molecular structure preprocessing, including SMILES handling, molecule standardization, and visualization.
⚠️ Important
Due to RDKit version constraints, this environment requires Python ≤ 3.7.

```bash
conda create -n rdkit2019 python==3.7
conda activate rdkit2019
conda install -c rdkit rdkit=2019.03.2 -y
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple epam.indigo
pip install ipykernel --upgrade
```
rdkit2019 environment requirements:

Python ≤ 3.7
RDKit 2019.03.2
EPAM Indigo


## 🤖 Model

The pretrained model file ```model.pt``` is required for running MetaReact.

File name: ```model.pt```

Download link:
https://zenodo.org/records/17898172

After downloading, please place model.pt into the corresponding task directory before executing the prediction scripts.

## 🚀 Use

MetaReact provides three prediction modes, corresponding to different application scenarios:

```enzyme-conditioned```

```enzyme-agnostic```

```enzyme-completion```

For each of the following directories:
```enzyme-conditioned```, ```enzyme-agnostic```, and ```enzyme-completion```
please follow the steps below.

### 📁 Step 1: Prepare Input Files

Place your input data into the file:

```test_drugs.csv```

Place the pretrained model file ```model.pt``` in the same directory.

### ▶️ Step 2: Run the Prediction Script

Execute the corresponding shell script:

```
predict_enzyme-conditioned.sh
predict_enzyme-agnostic.sh
predict_enzyme-completion.sh
```
### 📊 Step 3: Check Outputs

Prediction results are saved in:

```pred_results.csv```

In addition, visual diagrams are automatically generated to illustrate the structural transformations between the input substrate and its predicted metabolites 🧩➡️🧩


## 🧠 Prediction Modes Description
### 🧬 Enzyme-Conditioned

This mode is suitable for cases where both the substrate and enzyme are provided.

```Input: substrate + enzyme```

```Output: predicted metabolites```

This setting corresponds to the enzyme-conditioned prediction task.

### 🧪 Enzyme-Agnostic

This mode ignores enzyme information entirely.

```Input: substrate only```

```Output: predicted metabolites```

The model infers metabolic transformations without explicit enzyme constraints.

### 🔍 Enzyme-Completion

This mode predicts both enzymatic context and metabolic outcomes.

```Input: substrate only```

```Output: predicted enzymes and corresponding metabolites```

It is suitable for exploratory analysis when enzyme identity is unknown.

## ⚙️ Running Details

The MetaReact prediction pipeline consists of the following steps:
### Step 1: Data Augmentation

Apply Test-Time Augmentation (TTA) to the input data to improve robustness and reduce prediction variance.

### Step 2: SMILES Preprocessing

Kekulize SMILES strings and tokenize them into separate tokens suitable for model input.

### Step 3: Model Prediction

Use the trained MetaReact model to generate predictions based on the processed token sequences.

### Step 4: Output Conversion

Convert the model outputs from ReactSeq format into canonical SMILES representations, ensuring chemically valid and standardized structures.

