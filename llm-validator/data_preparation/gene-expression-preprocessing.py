import pandas as pd
import numpy as np

### Load Gene Expression Data ---
df = pd.read_csv('s3://sxx/Expression_Public_24Q4_subsetted.csv')
print(df.head())
print(df.shape)
print(df.columns[:20])
print(df.describe())

# --- Step 1: Filter for Lung Cancer Cell Lines Only ---
# Save metadata separately first
metadata_cols = ['depmap_id', 'cell_line_display_name', 
                 'lineage_1', 'lineage_2', 'lineage_3', 'lineage_4', 'lineage_6']

# Filter only lung cancer based on lineage_1
if 'lineage_1' in df.columns:
    df_lung = df[df['lineage_1'] == 'Lung'].copy()
else:
    raise ValueError("Lineage_1 column not found!")

print(f"Shape after filtering lung cancer cell lines: {df_lung.shape}")

# --- Step 2: Drop Metadata Columns ---
cols_to_drop = [col for col in metadata_cols if col in df_lung.columns]
df_lung = df_lung.drop(columns=cols_to_drop)

print("Shape after dropping metadata columns:", df_lung.shape)
print("Remaining columns (sample):", df_lung.columns[:10].tolist())

# --- Step 3: Log2 Transformation after adding Pseudo-count ---
pseudo_count = 1e-6  # very small to avoid changing high values
df_log2 = np.log2(df_lung + pseudo_count)

print("Log2 transformation done.")

# --- Step 4: Z-score Normalization Across Samples (Per Gene) ---
df_zscore = (df_log2 - df_log2.mean()) / df_log2.std()

print("Z-score normalization done.")

# --- Step 5: Select Top 5000 Variable Genes ---
gene_variances = df_zscore.var(axis=0)
top_5000_genes = gene_variances.sort_values(ascending=False).head(5000).index.tolist()

df_final = df_zscore[top_5000_genes]

print("Shape after selecting top 5000 genes:", df_final.shape)

# --- Step 6: Save to S3 ---
df_final.to_csv('s3://xx/preprocessed_gene_expression_top5000_lung_only.csv', index=False)
print("Saved preprocessed lung cancer gene expression data.")
