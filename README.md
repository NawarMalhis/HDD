# HDD
Prepare all figures and tables for the manuscript “The Hidden Disorder Divide: Reconciling Benchmark Inconsistencies in Intrinsically Disordered Protein Binding Site Prediction”.
https://doi.org/10.64898/2026.06.24.733783

# How to Generate Everything
Run these commands in the terminal (they assume your data is in Data/):
```bash
python Figure_01_ROC.py
python Figure_02_Violin_Training.py
python Figure_03_Histogram_Training.py
python Figure_03_Histogram_Test.py
python Figure_04A_Violin_IUPred.py
python Figure_04B_Violin_IUPred.py
python Figure_S01_ROC.py
python Figure_S02_Violin_IUPred.py
python Figure_S02_Violin_LIST.py

python Table_01_AUC.py
python Table_02_Merge_DBs_CAID23.py
python Table_03_Short_Long_AUC.py
python Table_S01_DataCompositions.py
python Table_S02_AUC.py
python Table_S03.py
python Table_S04_AUC_PDBvsIDR.py
python Table_S05_AUC_Merge_DBs_CAID23.py
python Table_S06_AUC_Merge_DBs_CAID1.py
```
All outputs go to Data/results/Figure_* and Data/results/Tables/.
