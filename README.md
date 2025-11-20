# MetaReact

## Installation

Environment 1

#### Create the environment using the `BBAR.yaml` file:

####  `conda env create -f BBAR.yaml`

Environment 2：rdkit2019 (for data processing related to rdkit and indigo)**

`conda create -n rdkit2019 python==3.7`

`conda activate rdkit2019`

`conda install -c rdkit rdkit=2019.03.2 -y`

`pip install -i https://pypi.tuna.tsinghua.edu.cn/simple epam.indigo`

`pip install ipykernel --upgrade`

`rdkit2019 requires:`

`\- python <=3.7`

## Use

For each of the folders — `enzyme-conditioned`, `enzyme-agnostic`, and `enzyme-completion` — please execute steps 1, 2, 3 and 4 sequentially.

`enzyme-conditioned` is suitable for cases where both the substrate and enzyme are provided, which is referred as "*enzyme-conditioned*".

`enzyme-agnostic` ignores enzyme information; the input is a substrate, and the output is the corresponding metabolites.

`enzyme-completion` takes the substrate as input and predicts both the possible enzymes and the resulting metabolites,which is referred as "*enzyme-completion*".

## model

model.pt