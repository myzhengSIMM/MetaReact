#!/bin/bash



# Step 1 and Step 2
conda activate rdkit2019
python step1_argue_enzyme.py
python step2_process.py
conda deactivate

# Step 3
conda activate BBAR
bash step3_test_neibu_best.sh
conda deactivate

# Step 4
conda activate rdkit2019
python step4_output_process.py
conda deactivate

conda activate BBAR
python step5_combined_list.py
python step6_result.py
python step7_draw_wo_score.py