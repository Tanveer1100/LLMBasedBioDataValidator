import pandas as pd
import numpy as np

### Load Gene Expression Data ---
df = pd.read_csv('s3://sxx/Expression_Public_24Q4_subsetted.csv')
print(df.head())
print(df.shape)
print(df.columns[:20]) 
print(df.describe())


# --- Step 1: Drop Metadata Columns ---
metadata_cols = ['depmap_id', 'cell_line_display_name', 
                 'lineage_1', 'lineage_2', 'lineage_3', 'lineage_4', 'lineage_6']
# Drop only if they exist
cols_to_drop = [col for col in metadata_cols if col in df.columns]
df = df.drop(columns=cols_to_drop)

print("Shape after dropping metadata columns:", df.shape)
print("Remaining columns (sample):", df.columns[:10].tolist())

# --- Step 2: Log2 Transformation after adding Pseudo-count ---
# Gene expression values might have zeros, add small pseudo-count first
pseudo_count = 1e-6  # very small to avoid changing high values
df_log2 = np.log2(df + pseudo_count)

print("Log2 transformation done.")

# --- Step 3: Z-score Normalization Across Samples (Per Gene) ---
df_zscore = (df_log2 - df_log2.mean()) / df_log2.std()

print("Z-score normalization done.")

gene_variances = df_zscore.var(axis=0)
top_5000_genes = gene_variances.sort_values(ascending=False).head(5000).index.tolist()

df_final = df_zscore[top_5000_genes]

print("Shape after selecting top 5000 genes:", df_final.shape)

df_final.to_csv('s3://xx/preprocessed_gene_expression_top5000.csv', index=False)
