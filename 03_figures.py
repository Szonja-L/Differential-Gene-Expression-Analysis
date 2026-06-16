"""
Step 3: Generate Publication-Quality Figures
- Volcano plot
- Heatmap of top DEGs  
- PCA plot
- Boxplots for key genes
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings('ignore')

# ── Color palette ──────────────────────────────────────────
C_TUMOR   = '#E8334A'   # vivid crimson
C_NORMAL  = '#2E86AB'   # steel blue
C_NS      = '#CCCCCC'   # light grey
C_UP      = '#E8334A'
C_DOWN    = '#2E86AB'
C_BG      = '#FAFAFA'
C_PANEL   = '#F2F2F2'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.facecolor': C_BG,
    'figure.facecolor': 'white',
    'axes.grid': False,
})

# Load data
results  = pd.read_csv('/home/claude/dge_project/results/dge_results_genes.csv')
expr     = pd.read_csv('/home/claude/dge_project/data/expression_matrix.csv', index_col=0)
groups   = pd.read_csv('/home/claude/dge_project/data/sample_groups.csv')

tumor_samples  = groups[groups['group'] == 'Tumor']['sample'].tolist()
normal_samples = groups[groups['group'] == 'Normal']['sample'].tolist()

print("Generating figures...")

# ──────────────────────────────────────────────────────────────
# FIGURE 1 — VOLCANO PLOT
# ──────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 7))
fig.patch.set_facecolor('white')
ax.set_facecolor(C_BG)

color_map = {
    'Upregulated in Tumor':   C_UP,
    'Downregulated in Tumor': C_DOWN,
    'Not Significant':        C_NS,
}
size_map = {
    'Upregulated in Tumor':   55,
    'Downregulated in Tumor': 55,
    'Not Significant':        18,
}
alpha_map = {
    'Upregulated in Tumor':   0.85,
    'Downregulated in Tumor': 0.85,
    'Not Significant':        0.35,
}
zorder_map = {
    'Upregulated in Tumor':   3,
    'Downregulated in Tumor': 3,
    'Not Significant':        1,
}

for direction, grp in results.groupby('direction'):
    ax.scatter(grp['log2FC'], grp['neg_log10_p'],
               c=color_map[direction], s=size_map[direction],
               alpha=alpha_map[direction], zorder=zorder_map[direction],
               linewidths=0, label=direction)

# Threshold lines
ax.axhline(-np.log10(0.05), color='#666666', ls='--', lw=1.2, alpha=0.7)
ax.axvline(-1.0,             color='#666666', ls='--', lw=1.2, alpha=0.7)
ax.axvline(+1.0,             color='#666666', ls='--', lw=1.2, alpha=0.7)

# Annotate top genes
top_up   = results[results['direction'] == 'Upregulated in Tumor'].nlargest(8, 'log2FC')
top_down = results[results['direction'] == 'Downregulated in Tumor'].nsmallest(8, 'log2FC')
for _, row in pd.concat([top_up, top_down]).iterrows():
    ax.annotate(row['gene'],
                xy=(row['log2FC'], row['neg_log10_p']),
                xytext=(row['log2FC'] + (0.25 if row['log2FC'] > 0 else -0.25),
                        row['neg_log10_p'] + 0.3),
                fontsize=8, fontweight='bold',
                color='#222222', ha='center',
                arrowprops=dict(arrowstyle='-', color='#888888', lw=0.7))

# Count annotations
n_up   = (results['direction'] == 'Upregulated in Tumor').sum()
n_down = (results['direction'] == 'Downregulated in Tumor').sum()
ax.text(0.97, 0.97, f'▲ {n_up} up', transform=ax.transAxes,
        ha='right', va='top', color=C_UP, fontweight='bold', fontsize=10)
ax.text(0.03, 0.97, f'▼ {n_down} down', transform=ax.transAxes,
        ha='left', va='top', color=C_DOWN, fontweight='bold', fontsize=10)

ax.set_xlabel('Log₂ Fold Change  (Tumor / Normal)', fontsize=12)
ax.set_ylabel('−log₁₀  (Adjusted p-value)', fontsize=12)
ax.set_title('Differential Gene Expression\nBreast Cancer Tumor vs. Adjacent Normal Tissue',
             fontsize=13, fontweight='bold', pad=12)

legend_elements = [
    mpatches.Patch(color=C_UP,   label=f'Upregulated in Tumor (n={n_up})'),
    mpatches.Patch(color=C_DOWN, label=f'Downregulated in Tumor (n={n_down})'),
    mpatches.Patch(color=C_NS,   label='Not Significant'),
]
ax.legend(handles=legend_elements, loc='lower right', framealpha=0.9, fontsize=9)
ax.text(1.05, -np.log10(0.05) + 0.15, 'adj.p = 0.05', fontsize=8, color='#666666')

plt.tight_layout()
plt.savefig('/home/claude/dge_project/figures/01_volcano_plot.png', dpi=180, bbox_inches='tight')
plt.close()
print("  ✓ Volcano plot saved")

# ──────────────────────────────────────────────────────────────
# FIGURE 2 — HEATMAP (top 30 DEGs)
# ──────────────────────────────────────────────────────────────
sig = results[results['significant']].copy()
top_genes = pd.concat([
    sig[sig['direction'] == 'Upregulated in Tumor'].nlargest(15, 'log2FC'),
    sig[sig['direction'] == 'Downregulated in Tumor'].nsmallest(15, 'log2FC'),
])

heatmap_data = expr.loc[top_genes['gene'], tumor_samples + normal_samples]

# Z-score per gene
heatmap_z = heatmap_data.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

fig, ax = plt.subplots(figsize=(12, 9))
fig.patch.set_facecolor('white')

im = ax.imshow(heatmap_z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)

# Axes
ax.set_yticks(range(len(heatmap_z)))
ax.set_yticklabels(heatmap_z.index, fontsize=9)

# Sample labels abbreviated
sample_labels = [s.split('_')[0] for s in heatmap_z.columns]
ax.set_xticks(range(len(heatmap_z.columns)))
ax.set_xticklabels(sample_labels, rotation=45, ha='right', fontsize=7.5)

# Color-code gene labels
for i, gene in enumerate(heatmap_z.index):
    color = C_UP if gene in top_genes[top_genes['direction']=='Upregulated in Tumor']['gene'].values else C_DOWN
    ax.get_yticklabels()[i].set_color(color)
    ax.get_yticklabels()[i].set_fontweight('bold')

# Group dividers
ax.axvline(len(tumor_samples) - 0.5, color='#333333', lw=2)

# Top labels
for x, label in [(len(tumor_samples)/2 - 0.5, 'Tumor'), 
                  (len(tumor_samples) + len(normal_samples)/2 - 0.5, 'Normal')]:
    ax.text(x, -0.8, label, ha='center', fontsize=11, fontweight='bold',
            color=C_TUMOR if label == 'Tumor' else C_NORMAL)

# Colorbar
cbar = plt.colorbar(im, ax=ax, pad=0.02, shrink=0.6)
cbar.set_label('Z-score', fontsize=10)

ax.set_title('Top 30 Differentially Expressed Genes\n(15 up + 15 down, z-score normalised)',
             fontsize=13, fontweight='bold', pad=15)
ax.set_xlabel('Samples', fontsize=11)

plt.tight_layout()
plt.savefig('/home/claude/dge_project/figures/02_heatmap.png', dpi=180, bbox_inches='tight')
plt.close()
print("  ✓ Heatmap saved")

# ──────────────────────────────────────────────────────────────
# FIGURE 3 — PCA PLOT
# ──────────────────────────────────────────────────────────────
from numpy.linalg import svd

# PCA on all samples using top 200 variable genes
all_real = [g for g in expr.index if not g.startswith('LOC')]
variances = expr.loc[all_real].var(axis=1)
top_var_genes = variances.nlargest(200).index

X = expr.loc[top_var_genes].T.values  # (n_samples, n_genes)
X_centered = X - X.mean(axis=0)

U, S, Vt = svd(X_centered, full_matrices=False)
PC1 = U[:, 0] * S[0]
PC2 = U[:, 1] * S[1]
explained = (S**2) / (S**2).sum()

fig, ax = plt.subplots(figsize=(8, 6.5))
fig.patch.set_facecolor('white')
ax.set_facecolor(C_BG)

sample_names = expr.columns.tolist()
for i, sample in enumerate(sample_names):
    is_tumor = sample in tumor_samples
    color  = C_TUMOR if is_tumor else C_NORMAL
    marker = 'o'
    ax.scatter(PC1[i], PC2[i], c=color, s=90, marker=marker,
               edgecolors='white', linewidths=0.8, zorder=3, alpha=0.9)

# Draw 95% confidence ellipses
def plot_ellipse(ax, x, y, color, label):
    from matplotlib.patches import Ellipse
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h = 2 * 2.0 * np.sqrt(vals)  # 2 std dev
    ell = Ellipse(xy=(np.mean(x), np.mean(y)), width=w, height=h, angle=theta,
                  color=color, alpha=0.12, zorder=1)
    ax.add_patch(ell)
    ell2 = Ellipse(xy=(np.mean(x), np.mean(y)), width=w, height=h, angle=theta,
                   fill=False, edgecolor=color, lw=1.5, linestyle='--', alpha=0.6, zorder=2)
    ax.add_patch(ell2)

tumor_idx  = [i for i, s in enumerate(sample_names) if s in tumor_samples]
normal_idx = [i for i, s in enumerate(sample_names) if s in normal_samples]

plot_ellipse(ax, PC1[tumor_idx],  PC2[tumor_idx],  C_TUMOR,  'Tumor')
plot_ellipse(ax, PC1[normal_idx], PC2[normal_idx], C_NORMAL, 'Normal')

legend_elements = [
    mpatches.Patch(color=C_TUMOR,  label=f'Tumor (n={len(tumor_idx)})'),
    mpatches.Patch(color=C_NORMAL, label=f'Normal (n={len(normal_idx)})'),
]
ax.legend(handles=legend_elements, fontsize=10, framealpha=0.9)

ax.set_xlabel(f'PC1  ({explained[0]*100:.1f}% variance)', fontsize=12)
ax.set_ylabel(f'PC2  ({explained[1]*100:.1f}% variance)', fontsize=12)
ax.set_title('Principal Component Analysis\nTop 200 Variable Genes', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig('/home/claude/dge_project/figures/03_pca_plot.png', dpi=180, bbox_inches='tight')
plt.close()
print("  ✓ PCA plot saved")

# ──────────────────────────────────────────────────────────────
# FIGURE 4 — BOXPLOTS for 6 key biomarker genes
# ──────────────────────────────────────────────────────────────
key_genes = ['MKI67', 'TOP2A', 'ERBB2', 'CDH1', 'ESR1', 'VEGFA']
labels_map = {
    'MKI67':  'MKI67\n(Proliferation)',
    'TOP2A':  'TOP2A\n(Cell cycle)',
    'ERBB2':  'ERBB2\n(HER2 signaling)',
    'CDH1':   'CDH1\n(Cell adhesion)',
    'ESR1':   'ESR1\n(ER signaling)',
    'VEGFA':  'VEGFA\n(Angiogenesis)',
}

fig, axes = plt.subplots(2, 3, figsize=(12, 7))
fig.patch.set_facecolor('white')

for ax, gene in zip(axes.flat, key_genes):
    ax.set_facecolor(C_BG)
    t_vals = expr.loc[gene, tumor_samples].values.astype(float)
    n_vals = expr.loc[gene, normal_samples].values.astype(float)
    
    bp = ax.boxplot([n_vals, t_vals],
                    patch_artist=True, widths=0.45,
                    medianprops=dict(color='white', lw=2.5),
                    whiskerprops=dict(color='#444444', lw=1.2),
                    capprops=dict(color='#444444', lw=1.5),
                    flierprops=dict(marker='o', ms=5, alpha=0.5))
    
    bp['boxes'][0].set_facecolor(C_NORMAL)
    bp['boxes'][0].set_alpha(0.8)
    bp['boxes'][1].set_facecolor(C_TUMOR)
    bp['boxes'][1].set_alpha(0.8)
    
    # Jitter overlay
    for vals, x_pos, color in [(n_vals, 1, C_NORMAL), (t_vals, 2, C_TUMOR)]:
        jitter = np.random.normal(0, 0.07, len(vals))
        ax.scatter(x_pos + jitter, vals, c=color, s=28, alpha=0.75, zorder=3,
                   edgecolors='white', linewidths=0.5)
    
    # Stats annotation
    _, p = stats.ttest_ind(t_vals, n_vals, equal_var=False)
    lfc = results.loc[results['gene'] == gene, 'log2FC'].values[0]
    stars = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
    y_max = max(max(t_vals), max(n_vals)) + 0.2
    ax.plot([1, 1, 2, 2], [y_max, y_max+0.1, y_max+0.1, y_max], 'k-', lw=1)
    ax.text(1.5, y_max + 0.15, f'{stars}\nlog₂FC={lfc:+.2f}',
            ha='center', fontsize=8.5, fontweight='bold')
    
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Normal', 'Tumor'], fontsize=10)
    ax.set_ylabel('Expression (log₂)', fontsize=9)
    ax.set_title(labels_map[gene], fontsize=10, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle('Expression of Key Biomarker Genes\nTumor vs. Normal Tissue',
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('/home/claude/dge_project/figures/04_biomarker_boxplots.png',
            dpi=180, bbox_inches='tight')
plt.close()
print("  ✓ Biomarker boxplots saved")

print("\nAll figures saved to /home/claude/dge_project/figures/")
print("Done!")
