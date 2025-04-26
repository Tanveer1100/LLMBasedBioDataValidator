```python
import pandas as pd
import numpy as np


# Load Gene Expression Data ---

df = pd.read_csv('s3://sxx/Expression_Public_24Q4_subsetted.csv')
print(df.head())
print(df.shape)
print(df.columns[:20]) 
print(df.describe())
```

    /tmp/ipykernel_29252/3084145054.py:7: DtypeWarning: Columns (5) have mixed types. Specify dtype option on import or set low_memory=False.
      df = pd.read_csv('s3://shelfspace-alpha-sandbox/users/tyshaikh/llm-validator/Expression_Public_24Q4_subsetted.csv')


        depmap_id cell_line_display_name    lineage_1  \
    0  ACH-001435                1156QE8       Testis   
    1  ACH-001270                 127399  Soft Tissue   
    2  ACH-002680                170MGBA    CNS/Brain   
    3  ACH-001438             1777NRPMET       Testis   
    4  ACH-002401                  21MT2       Breast   
    
                              lineage_2                         lineage_3  \
    0  Non-Seminomatous Germ Cell Tumor               Embryonal Carcinoma   
    1                  Synovial Sarcoma                  Synovial Sarcoma   
    2                    Diffuse Glioma                      Glioblastoma   
    3  Non-Seminomatous Germ Cell Tumor               Embryonal Carcinoma   
    4         Invasive Breast Carcinoma  Breast Invasive Ductal Carcinoma   
    
      lineage_6  lineage_4    TSPAN6      TNMD      DPM1  ...   SPDYE11      H3C2  \
    0       NaN        NaN  5.255123  0.000000  6.061776  ...  0.000000  0.536053   
    1       NaN        NaN  3.243364  0.000000  6.863691  ...  0.014355  1.646163   
    2       NaN        NaN  4.155425  0.124328  6.681730  ...  0.000000  0.659925   
    3       NaN        NaN  4.438293  0.000000  7.098243  ...  0.000000  1.304511   
    4     HER2+        NaN  4.363872  0.028569  6.704042  ...  0.000000  0.659925   
    
           H3C3  DUS4L-BCAP29  C8orf44-SGK3  ELOA3BP    NPBWR1  ELOA3DP  ELOA3P  \
    0  0.000000      1.599318      0.321928      0.0  0.014355      0.0     0.0   
    1  1.480265      1.978196      0.526069      0.0  0.000000      0.0     0.0   
    2  0.150560      2.432959      0.731183      0.0  0.000000      0.0     0.0   
    3  0.000000      2.321928      0.238787      0.0  0.000000      0.0     0.0   
    4  0.545968      1.310340      1.673556      0.0  0.765535      0.0     0.0   
    
           CDR1  
    0  0.000000  
    1  0.000000  
    2  0.495695  
    3  0.000000  
    4  0.000000  
    
    [5 rows x 19160 columns]
    (1673, 19160)
    Index(['depmap_id', 'cell_line_display_name', 'lineage_1', 'lineage_2',
           'lineage_3', 'lineage_6', 'lineage_4', 'TSPAN6', 'TNMD', 'DPM1',
           'SCYL3', 'FIRRM', 'FGR', 'CFH', 'FUCA2', 'GCLC', 'NFYA', 'STPG1',
           'NIPAL3', 'LAS1L'],
          dtype='object')
           lineage_4       TSPAN6         TNMD         DPM1        SCYL3  \
    count        0.0  1673.000000  1673.000000  1673.000000  1673.000000   
    mean         NaN     3.364755     0.084035     6.511230     2.351832   
    std          NaN     1.646095     0.421914     0.650400     0.544043   
    min          NaN     0.000000     0.000000     3.117695     0.594549   
    25%          NaN     2.873813     0.000000     6.104127     1.992768   
    50%          NaN     3.808385     0.000000     6.503190     2.327687   
    75%          NaN     4.426936     0.000000     6.928370     2.653060   
    max          NaN     8.132680     7.814422     9.175250     4.747387   
    
                 FIRRM          FGR          CFH        FUCA2         GCLC  ...  \
    count  1673.000000  1673.000000  1673.000000  1673.000000  1673.000000  ...   
    mean      3.621073     0.431686     2.231446     5.104079     4.561802  ...   
    std       0.834547     1.216362     2.303087     1.840251     1.149210  ...   
    min       0.056584     0.000000     0.000000     0.000000     1.157044  ...   
    25%       3.158660     0.014355     0.124328     4.657640     3.855990  ...   
    50%       3.704872     0.042644     1.321928     5.653920     4.463361  ...   
    75%       4.165108     0.150560     4.099295     6.275380     5.190615  ...   
    max       5.972463     8.014299     9.693016     8.633068     9.581991  ...   
    
               SPDYE11         H3C2         H3C3  DUS4L-BCAP29  C8orf44-SGK3  \
    count  1673.000000  1673.000000  1673.000000   1673.000000   1673.000000   
    mean      0.032003     1.015787     0.750993      2.059848      0.405048   
    std       0.079981     0.732034     0.777924      0.686985      0.367234   
    min       0.000000     0.000000     0.000000      0.263034      0.000000   
    25%       0.000000     0.475085     0.000000      1.580145      0.124328   
    50%       0.014355     0.887525     0.555816      2.017922      0.333424   
    75%       0.042644     1.400538     1.195348      2.503349      0.594549   
    max       2.066950     5.173127     4.894818      4.433627      2.247928   
    
               ELOA3BP       NPBWR1      ELOA3DP       ELOA3P         CDR1  
    count  1673.000000  1673.000000  1673.000000  1673.000000  1673.000000  
    mean      0.016839     0.264674     0.017209     0.023767     0.063297  
    std       0.052869     0.592130     0.060571     0.078146     0.184181  
    min       0.000000     0.000000     0.000000     0.000000     0.000000  
    25%       0.000000     0.000000     0.000000     0.000000     0.000000  
    50%       0.000000     0.014355     0.000000     0.000000     0.000000  
    75%       0.000000     0.189034     0.000000     0.014355     0.000000  
    max       0.622930     4.254745     0.992768     1.056584     1.867896  
    
    [8 rows x 19154 columns]



```python
# --- Step 1: Drop Metadata Columns ---
metadata_cols = ['depmap_id', 'cell_line_display_name', 
                 'lineage_1', 'lineage_2', 'lineage_3', 'lineage_4', 'lineage_6']
# Drop only if they exist
cols_to_drop = [col for col in metadata_cols if col in df.columns]
df = df.drop(columns=cols_to_drop)


print("Shape after dropping metadata columns:", df.shape)
print("Remaining columns (sample):", df.columns[:10].tolist())
```

    Shape after dropping metadata columns: (1673, 19153)
    Remaining columns (sample): ['TSPAN6', 'TNMD', 'DPM1', 'SCYL3', 'FIRRM', 'FGR', 'CFH', 'FUCA2', 'GCLC', 'NFYA']



```python
# --- Step 2: Log2 Transformation after adding Pseudo-count ---
# Gene expression values might have zeros, add small pseudo-count first
pseudo_count = 1e-6  # very small to avoid changing high values
df_log2 = np.log2(df + pseudo_count)

print("Log2 transformation done.")
```

    Log2 transformation done.



```python
# --- Step 3: Z-score Normalization Across Samples (Per Gene) ---
df_zscore = (df_log2 - df_log2.mean()) / df_log2.std()

print("Z-score normalization done.")
```

    Z-score normalization done.



```python
gene_variances = df_zscore.var(axis=0)
top_5000_genes = gene_variances.sort_values(ascending=False).head(5000).index.tolist()

df_final = df_zscore[top_5000_genes]

print("Shape after selecting top 5000 genes:", df_final.shape)
```

    Shape after selecting top 5000 genes: (1673, 5000)



```python
df_final.to_csv('s3://xx/preprocessed_gene_expression_top5000.csv', index=False)
```


```python

```

print("Remaining columns (sample):", df.columns[:10].tolist())
