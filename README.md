# HDD
Prepare all figures and tables for the manuscript “The Hidden Disorder Divide: Reconciling Benchmark Inconsistencies in Intrinsically Disordered Protein Binding Site Prediction”.
https://doi.org/10.64898/2026.06.24.733783

## Install
```bash
# Clone HDD
git clone https://github.com/NawarMalhis/HDD.git
# Clone the "annotated fasta format" library:	
git clone https://github.com/NawarMalhis/AFF.git
# Change directory:	
cd HDD
# Update the path to the AFF (annotated fasta format) folder in param.py.
aff_path = '/xxx/xxx/AFF/'
# Create the hdd environment:
conda env create -f hdd.yml
```

## How to Generate Everything


```bash
# First activate the hdd environment:
conda activate hdd

# Then run these commands in the terminal; they assume your data is in Data/ and
#   all outputs go to Data/results/Figure_* and Data/results/Tables/:
python3 Figure_01_ROC.py
python3 Figure_02_Violin_Training.py
python3 Figure_03_Left.py
python3 Figure_03_Right.py

python3 Figure_S01_ROC.py
python3 Figure_S02_Histogram_Test.py
python3 Figure_S02_Histogram_Training.py
python3 Figure_S03.py

python3 Table_01_AUC.py
python3 Table_02_Merge.py

python3 Table_S01_DataCompositions.py
python3 Table_S02_AUC.py
python3 Table_S03.py
python3 Table_S04_AUC_PDBvsIDR.py
python3 Table_S05_AUC_Merge_DBs_CAID23.py
python3 Table_S06_AUC_Merge_DBs_CAID1.py
python3 Table_S07_AUC_Merge_CAID1_CAID23.py
python3 Table_S08_Short_Long_AUC.py
```
