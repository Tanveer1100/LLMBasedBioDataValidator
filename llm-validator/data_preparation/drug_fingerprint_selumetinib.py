import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Step 1: Load the PubChem CSV
df = pd.read_csv('s3://xx/PubChem_compound_text_Selumetinib.csv')
# Step 2: Inspect what columns you have
print("Columns:", df.columns.tolist())
print("First few rows:\n", df.head())

# Step 3: Extract SMILES (assuming the field is named 'CanonicalSMILES' or similar)
smiles_column_candidates = ['CanonicalSMILES', 'SMILES', 'IsomericSMILES']

for col in smiles_column_candidates:
    if col in df.columns:
        smiles = df[col].values[0]
        print("Found SMILES:", smiles)
        break
else:
    print("❗No SMILES field found!")

# Step 4: Convert SMILES to ECFP4 fingerprint
mol = Chem.MolFromSmiles(smiles)
fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
fingerprint_array = list(fp)

import numpy as np
fingerprint_df = pd.DataFrame([fingerprint_array], index=['Selumetinib'])

fingerprint_df.to_csv('s3://xx/drug_fingerprint_selumetinib.csv')
