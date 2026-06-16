"""
Step 2: Differential Gene Expression Analysis
- Welch's t-test (appropriate for unequal sample sizes)
- Benjamini-Hochberg FDR correction
- Log2 fold change calculation
- Volcano plot data preparation
"""

import pandas as pd
import numpy as np
from scipy import stats

print("=" * 60)
print("DIFFERENTIAL GENE EXPRESSION ANALYSIS")
print("Dataset: Breast Cancer Tumor vs. Adjacent Normal Tissue")
print("Modeled after: GSE7904 (Affymetrix HG-U133Plus2)")
print("=" * 60)

# Load data
expr = pd.read_csv('/home/claude/dge_project/data/expression_matrix.csv', index_col=0)
groups = pd.read_csv('/home/claude/dge_project/data/sample_groups.csv')

tumor_samples = groups[groups['group'] == 'Tumor']['sample'].tolist()
normal_samples = groups[groups['group'] == 'Normal']['sample'].tolist()

print(f"\nTumor samples: {len(tumor_samples)}")
print(f"Normal samples: {len(normal_samples)}")

# Calculate DGE statistics
results = []
for gene in expr.index:
    tumor_vals = expr.loc[gene, tumor_samples].values.astype(float)
    normal_vals = expr.loc[gene, normal_samples].values.astype(float)
    
    # Log2 fold change: tumor vs normal
    log2fc = np.mean(tumor_vals) - np.mean(normal_vals)  # already log2 scale
    
    # Welch's t-test (unequal variance)
    t_stat, p_val = stats.ttest_ind(tumor_vals, normal_vals, equal_var=False)
    
    mean_tumor = np.mean(tumor_vals)
    mean_normal = np.mean(normal_vals)
    
    results.append({
        'gene': gene,
        'log2FC': round(log2fc, 4),
        'mean_tumor': round(mean_tumor, 4),
        'mean_normal': round(mean_normal, 4),
        'p_value': p_val,
        't_statistic': round(t_stat, 4)
    })

results_df = pd.DataFrame(results)

# Benjamini-Hochberg FDR correction
results_df = results_df.sort_values('p_value')
n = len(results_df)
results_df['rank'] = range(1, n + 1)
results_df['adj_p_value'] = results_df['p_value'] * n / results_df['rank']
results_df['adj_p_value'] = results_df['adj_p_value'].clip(upper=1.0)

# Correct BH: ensure monotonicity
adj_p = results_df['adj_p_value'].values.copy()
for i in range(len(adj_p) - 2, -1, -1):
    adj_p[i] = min(adj_p[i], adj_p[i + 1])
results_df['adj_p_value'] = adj_p

# Significance classification
results_df['neg_log10_p'] = -np.log10(results_df['adj_p_value'].clip(lower=1e-15))
results_df['significant'] = (results_df['adj_p_value'] < 0.05) & (results_df['log2FC'].abs() >= 1.0)
results_df['direction'] = np.where(
    ~results_df['significant'], 'Not Significant',
    np.where(results_df['log2FC'] > 0, 'Upregulated in Tumor', 'Downregulated in Tumor')
)

# Filter to real gene symbols (exclude LOC background)
real_genes = results_df[~results_df['gene'].str.startswith('LOC')].copy()

print(f"\n--- SUMMARY ---")
print(f"Total genes tested: {len(results_df)}")
print(f"Significant DEGs (adj.p<0.05, |log2FC|>=1): {results_df['significant'].sum()}")
print(f"  Upregulated in tumor:   {(results_df['direction']=='Upregulated in Tumor').sum()}")
print(f"  Downregulated in tumor: {(results_df['direction']=='Downregulated in Tumor').sum()}")

print(f"\n--- TOP 15 UPREGULATED IN TUMOR ---")
up = real_genes[real_genes['direction'] == 'Upregulated in Tumor'].nlargest(15, 'log2FC')
for _, row in up.iterrows():
    print(f"  {row['gene']:<12} log2FC={row['log2FC']:+.2f}  adj.p={row['adj_p_value']:.2e}")

print(f"\n--- TOP 15 DOWNREGULATED IN TUMOR ---")
down = real_genes[real_genes['direction'] == 'Downregulated in Tumor'].nsmallest(15, 'log2FC')
for _, row in down.iterrows():
    print(f"  {row['gene']:<12} log2FC={row['log2FC']:+.2f}  adj.p={row['adj_p_value']:.2e}")

# Save results
results_df.to_csv('/home/claude/dge_project/results/dge_results_all.csv', index=False)
real_genes.to_csv('/home/claude/dge_project/results/dge_results_genes.csv', index=False)
real_genes[real_genes['significant']].to_csv('/home/claude/dge_project/results/dge_significant.csv', index=False)

print(f"\nResults saved to /home/claude/dge_project/results/")
