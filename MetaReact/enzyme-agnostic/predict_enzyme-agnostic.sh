#!/bin/bash


conda activate rdkit2019
python step1_argue.py
python step2_process.py
conda deactivate


conda activate BBAR
bash step3_test_neibu_best.sh
conda deactivate

conda activate rdkit2019
python step4_result.py
conda deactivate
conda activate BBAR
python step5_draw.py
# python step4_result_add_SoM.py

conda deactivate
