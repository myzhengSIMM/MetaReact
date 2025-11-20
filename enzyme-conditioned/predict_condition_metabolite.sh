#!/bin/bash

# Step 1 and Step 2
conda activate rdkit2019
python step1_argue.py
python step2_process_input.py
conda deactivate

# Step 3
conda activate BBAR
bash step3_test_best.sh
conda deactivate

# Step 4
conda activate rdkit2019
python step4_calculate_result.py
conda deactivate

# conda activate BBAR
# python step5_draw.py
