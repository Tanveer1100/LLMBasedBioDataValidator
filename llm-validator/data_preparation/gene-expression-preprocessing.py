import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

# Example: Input gene expression DataFrame
# Rows = samples (cell lines), Columns = genes
np.random.seed(42)
gene_expression_df = pd.DataFrame(
    np.random.rand(100, 20000),  # 100 samples, 20,000 genes
    columns=[f"Gene_{i}" for i in range(20000)]
)

# Step 1: Log2 transformation with pseudo-count addition
pseudo_count = 1
log2_transformed = np.log2(gene_expression_df + pseudo_count)

# Step 2: Z-score normalization per gene
scaler = StandardScaler(with_mean=True, with_std=True)
zscore_normalized = pd.DataFrame(
    scaler.fit_transform(log2_transformed.T).T,  # Normalize across samples, per gene
    index=log2_transformed.index,
    columns=log2_transformed.columns
)

# Step 3: Selection of top 5000 most variable genes
gene_variances = zscore_normalized.var(axis=0)
top_5000_genes = gene_variances.sort_values(ascending=False).head(5000).index

# Final preprocessed DataFrame
filtered_gene_expression_df = zscore_normalized[top_5000_genes]

# Preview
filtered_gene_expression_df.head()
