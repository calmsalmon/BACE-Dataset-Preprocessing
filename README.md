# BACE Dataset Preprocessing
A program the preprocessed the BACE dataset from DeepChem to be used for testing of a Quantum-hybrid model.

## Dataset
The dataset is from DeepChem and is the BACE regression dataset. Due to dependency issues, a downloaded version was used for preprocessing.

## Processing Steps
1. Remove empty columns and rows missing required information
2. Canonicalize SMILES and remove invalid molecules
3. Remove duplicate molecules
4. Convert experiment data into numerical values
5. Remove duplicate descriptor columns
6. Replace missing descriptor values
7. Reorder columns
8. Save data
9. Print logistics 
