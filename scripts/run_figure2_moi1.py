#!/usr/bin/env python
"""
Run Figure 2 Analysis for MOI=1 data
Generates outputs with _moi1 suffix
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import ElasticNet
from sklearn.cluster import KMeans
import scipy.sparse as sp_sparse
import scipy.stats as sp_stats
import os

print("=" * 60)
print("FIGURE 2 ANALYSIS - MOI=1 DATA")
print("=" * 60)

# Set plotting style
plt.style.use('default')
sc.settings.verbosity = 3

# Paths
full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
adata_arr_pn = os.path.join(full_exp_pn, 'adata_arrs')
sequencing_info_pn = os.path.join(full_exp_pn, 'sequencing_info')
CROP_pn = os.path.join(full_exp_pn, 'CROP')
design_mat_pn = os.path.join(CROP_pn, 'design_mats')
lm_pn = os.path.join(CROP_pn, 'linear_model')

# Output directories
output_dir = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me/results/figure2_moi1'
fig_dir = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me/figures/figure2_moi1'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {fig_dir}")

# Load data
print("\n" + "=" * 60)
print("Loading data...")
print("=" * 60)

adata_full = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))
print(f"Full dataset: {adata_full.shape[0]} cells, {adata_full.shape[1]} features")

# Filter to MOI=1
moi_numeric = adata_full.obs['MOI'].astype(float)
moi1_mask = moi_numeric == 1
adata = adata_full[moi1_mask].copy()
print(f"MOI=1 dataset: {adata.shape[0]} cells")
print(f"\nConditions:")
print(adata.obs['Condition'].value_counts())

# Load target names
target_names = np.load(os.path.join(design_mat_pn, 'target_names.npy'))
print(f"\nTotal targets in library: {len(target_names)}")

# ============================================================
# Figure 2A: Cells per Target
# ============================================================
print("\n" + "=" * 60)
print("Figure 2A: Cells per Target Analysis")
print("=" * 60)

# Load sgRNA information
gw_guides_pn = os.path.join(sequencing_info_pn, 'gw_lib_guides.csv')
target_sgRNA_df = pd.read_csv(gw_guides_pn)
target_sgRNA_names = target_sgRNA_df['sgRNA_name'].to_numpy().astype(str)

control_sgRNA_df = pd.read_csv(os.path.join(sequencing_info_pn, 'human_TF_CTL_guides.txt'),
                               sep='\t', header=None, names=['sgRNA_name', 'sgRNA_barcode'])
control_sgRNA_names = control_sgRNA_df['sgRNA_name'].to_numpy().astype(str)
combined_sgRNA_names = np.append(target_sgRNA_names, control_sgRNA_names)

# Get sgRNA assignments
sgRNA_assignments = adata.obs['sgRNA'].to_numpy().astype(str)
sample_arr = adata.obs['Sample'].to_numpy().astype(str)

high_mask = np.char.startswith(sample_arr, 'High_')
low_mask = np.char.startswith(sample_arr, 'Low_')

print(f"HLA_High cells: {high_mask.sum()}")
print(f"HLA_Low cells: {low_mask.sum()}")

def count_cells_per_target(sgRNA_assignments, target_names, combined_sgRNA_names):
    cells_per_target = np.zeros(len(target_names))
    for target_i, target in enumerate(target_names):
        target_sgRNAs = [s for s in combined_sgRNA_names if s.startswith(target + '_')]
        cells_per_target[target_i] = np.isin(sgRNA_assignments, target_sgRNAs).sum()
    return cells_per_target

print("Counting cells per target...")
cells_per_target_high = count_cells_per_target(sgRNA_assignments[high_mask], target_names, combined_sgRNA_names)
cells_per_target_low = count_cells_per_target(sgRNA_assignments[low_mask], target_names, combined_sgRNA_names)

print(f"Total cells in HLA_High targets: {cells_per_target_high.sum():.0f}")
print(f"Total cells in HLA_Low targets: {cells_per_target_low.sum():.0f}")

# Create DataFrame
cells_per_target_df = pd.DataFrame({
    'target': target_names,
    'cells_high': cells_per_target_high,
    'cells_low': cells_per_target_low,
    'cells_total': cells_per_target_high + cells_per_target_low
})

cells_df_filt = cells_per_target_df[
    ~cells_per_target_df['target'].str.contains('NO_SITE|NON-GENE', na=False)
].copy()

print(f"\nTop 10 targets by HLA_High cells:")
print(cells_df_filt.nlargest(10, 'cells_high')[['target', 'cells_high', 'cells_low']])

print(f"\nTop 10 targets by HLA_Low cells:")
print(cells_df_filt.nlargest(10, 'cells_low')[['target', 'cells_low', 'cells_high']])

# Plot Figure 2A
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
top_low = cells_df_filt.nlargest(30, 'cells_low')
ax.barh(range(len(top_low)), top_low['cells_low'].values, color='steelblue')
ax.set_yticks(range(len(top_low)))
ax.set_yticklabels(top_low['target'].values)
ax.invert_yaxis()
ax.set_xlabel('Number of Cells')
ax.set_title(f'HLA_Low: Top 30 Targets\n(n={int(cells_per_target_low.sum())} total cells)')

ax = axes[1]
top_high = cells_df_filt.nlargest(30, 'cells_high')
ax.barh(range(len(top_high)), top_high['cells_high'].values, color='coral')
ax.set_yticks(range(len(top_high)))
ax.set_yticklabels(top_high['target'].values)
ax.invert_yaxis()
ax.set_xlabel('Number of Cells')
ax.set_title(f'HLA_High: Top 30 Targets\n(n={int(cells_per_target_high.sum())} total cells)')

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Fig2A_cells_per_target_moi1.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, 'Fig2A_cells_per_target_moi1.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved Figure 2A")

# Save cells per target
cells_per_target_df.to_csv(os.path.join(output_dir, 'cells_per_target_moi1.csv'), index=False)
np.save(os.path.join(output_dir, 'cells_per_target_high_moi1.npy'), cells_per_target_high)
np.save(os.path.join(output_dir, 'cells_per_target_low_moi1.npy'), cells_per_target_low)

# ============================================================
# Figure 2D: HLA Program Score
# ============================================================
print("\n" + "=" * 60)
print("Figure 2D: HLA Program Score")
print("=" * 60)

hla_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1', 'TAP2', 'TAPBP', 'PSMB8', 'PSMB9']
available_hla_genes = [g for g in hla_genes if g in adata.var_names]
print(f"Available HLA program genes: {available_hla_genes}")

hla_expr = adata[:, available_hla_genes].X
if sp_sparse.issparse(hla_expr):
    hla_expr = hla_expr.toarray()

hla_zscore = sp_stats.zscore(hla_expr, axis=0)
hla_program_score = np.nanmean(hla_zscore, axis=1)
adata.obs['HLA_program_score'] = hla_program_score

# Plot Figure 2D
fig, ax = plt.subplots(figsize=(8, 5))

# Map condition names (the data uses 'High'/'Low' not 'HLA_High'/'HLA_Low')
condition_map = {'HLA_Low': 'Low', 'Control': 'Control', 'HLA_High': 'High'}
conditions = ['HLA_Low', 'Control', 'HLA_High']
colors = ['steelblue', 'gray', 'coral']

for i, cond in enumerate(conditions):
    data_cond = condition_map[cond]
    mask = adata.obs['Condition'] == data_cond
    scores = adata.obs.loc[mask, 'HLA_program_score'].dropna()

    if len(scores) > 0:
        parts = ax.violinplot([scores], positions=[i], showmeans=True, showmedians=False)
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i])
            pc.set_alpha(0.7)

        mean_score = scores.mean()
        n_cells = len(scores)
        ax.text(i, scores.max() + 0.3, f'n={n_cells}\nmean={mean_score:.2f}', ha='center', fontsize=9)

ax.set_xticks(range(len(conditions)))
ax.set_xticklabels(conditions)
ax.set_ylabel('HLA Program Score (z-scored)')
ax.set_title('Figure 2D: HLA Program Score by Condition (MOI=1)')
ax.axhline(0, color='black', linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Fig2D_HLA_program_score_moi1.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, 'Fig2D_HLA_program_score_moi1.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved Figure 2D")

# Statistical test
high_scores = adata.obs.loc[adata.obs['Condition'] == 'High', 'HLA_program_score'].dropna()
low_scores = adata.obs.loc[adata.obs['Condition'] == 'Low', 'HLA_program_score'].dropna()

stat, pval = sp_stats.mannwhitneyu(high_scores, low_scores, alternative='two-sided')
pooled_std = np.sqrt((high_scores.std()**2 + low_scores.std()**2) / 2)
cohens_d = (high_scores.mean() - low_scores.mean()) / pooled_std

print(f"\nHLA_High vs HLA_Low comparison:")
print(f"  HLA_High: mean={high_scores.mean():.3f}, n={len(high_scores)}")
print(f"  HLA_Low:  mean={low_scores.mean():.3f}, n={len(low_scores)}")
print(f"  Mann-Whitney U p-value: {pval:.2e}")
print(f"  Cohen's d: {cohens_d:.3f}")

# Save HLA program scores
hla_scores_df = adata.obs[['Condition', 'HLA_program_score']].copy()
hla_scores_df.to_csv(os.path.join(output_dir, 'HLA_program_scores_moi1.csv'))

# ============================================================
# Figure 2E: Beta Matrix and Clustering
# ============================================================
print("\n" + "=" * 60)
print("Figure 2E: Beta Matrix and Clustering")
print("=" * 60)

cells_per_target_thresh = 17
num_clusters_covs = 7
num_clusters_features = 9

# Load existing beta matrix
existing_lm_path = os.path.join(lm_pn, f'{cells_per_target_thresh}_cells_per_target/all_features')
print(f"Loading from: {existing_lm_path}")

B_mat = np.load(os.path.join(existing_lm_path, 'EN_B_EM.npy'))
cov_names = np.load(os.path.join(existing_lm_path, 'cov_names.npy'))
feature_names = np.load(os.path.join(existing_lm_path, 'feature_names.npy'))

print(f"Beta matrix shape: {B_mat.shape}")
print(f"  Features (rows): {B_mat.shape[0]}")
print(f"  Covariates (columns): {B_mat.shape[1]}")

# Helper functions
def calc_corr_matrix(mat):
    corr_mat = pd.DataFrame(mat).corr().to_numpy()
    corr_mat[np.isnan(corr_mat)] = 0
    return corr_mat

def cluster_by_kmeans(mat, num_clusters):
    km = KMeans(n_clusters=num_clusters, random_state=3, n_init=10)
    km.fit(mat)

    ret_order = np.array([], dtype=int)
    for clust_i in range(num_clusters):
        curr_clust_inds = np.where(clust_i == km.labels_)[0]
        ret_order = np.append(ret_order, curr_clust_inds.astype(int))

    return ret_order, km.labels_[ret_order], km.labels_

def remove_sparse_rows_and_columns(mat, row_names, col_names,
                                    approx_zero=0.02, row_sparse_thresh=0.7, col_sparse_thresh=0.7):
    row_sparsity = np.array([(np.abs(mat[r, :]) <= approx_zero).sum() / len(col_names)
                             for r in range(len(row_names))])
    good_rows = np.where(row_sparsity < row_sparse_thresh)[0]
    filtered_rows = row_names[good_rows]

    col_sparsity = np.array([(np.abs(mat[:, c]) <= approx_zero).sum() / len(row_names)
                             for c in range(len(col_names))])
    good_cols = np.where(col_sparsity < col_sparse_thresh)[0]
    filtered_cols = col_names[good_cols]

    filtered_mat = mat[good_rows, :][:, good_cols]

    return filtered_mat, filtered_rows, filtered_cols, good_rows, good_cols

# Filter sparse rows/columns
filtered_B, filtered_features, filtered_covs, row_idx, col_idx = remove_sparse_rows_and_columns(
    B_mat, feature_names, cov_names,
    approx_zero=0.05, col_sparse_thresh=0.8, row_sparse_thresh=0.25
)

print(f"\nAfter filtering sparse rows/columns:")
print(f"  Features: {len(feature_names)} -> {len(filtered_features)}")
print(f"  Covariates: {len(cov_names)} -> {len(filtered_covs)}")

# Calculate correlation matrices
print("\nCalculating correlation matrices...")
B_correlation_covs = calc_corr_matrix(filtered_B)
B_correlation_features = calc_corr_matrix(filtered_B.T)

# K-means clustering
print(f"Clustering covariates into {num_clusters_covs} modules...")
covs_order, covs_cluster_labels_ordered, covs_cluster_labels = cluster_by_kmeans(
    B_correlation_covs, num_clusters_covs
)

print(f"Clustering features into {num_clusters_features} modules...")
features_order, features_cluster_labels_ordered, features_cluster_labels = cluster_by_kmeans(
    B_correlation_features, num_clusters_features
)

print(f"\nTarget module sizes:")
for i in range(num_clusters_covs):
    n = (covs_cluster_labels == i).sum()
    print(f"  Module {i}: {n} targets")

print(f"\nFeature module sizes:")
for i in range(num_clusters_features):
    n = (features_cluster_labels == i).sum()
    print(f"  Module {i}: {n} features")

# Plot beta matrix heatmap
fig, ax = plt.subplots(figsize=(12, 10))
filtered_plot_mat = filtered_B[features_order, :][:, covs_order]
im = ax.imshow(filtered_plot_mat,
               aspect=(filtered_plot_mat.shape[1]/filtered_plot_mat.shape[0]),
               cmap='bwr', vmin=-1, vmax=1)
cbar = fig.colorbar(im, ax=ax, extend='both', shrink=0.6)
cbar.set_label('Beta coefficient')
ax.set_xlabel('Target Perturbations (clustered)')
ax.set_ylabel('Features (clustered)')
ax.set_title(f'Figure 2E: Beta Matrix\n({num_clusters_features} feature modules x {num_clusters_covs} target modules)')
ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Fig2E_beta_matrix_moi1.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, 'Fig2E_beta_matrix_moi1.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved Figure 2E beta matrix")

# Plot target correlation heatmap
fig, ax = plt.subplots(figsize=(10, 10))
corr_plot = B_correlation_covs[covs_order, :][:, covs_order]
im = ax.imshow(corr_plot, cmap='PRGn', vmin=-0.2, vmax=0.2)
cbar = fig.colorbar(im, ax=ax, extend='both', shrink=0.6)
cbar.set_label('Correlation')
ax.set_xlabel('Target Perturbations')
ax.set_ylabel('Target Perturbations')
ax.set_title(f'Target Correlation Matrix ({num_clusters_covs} modules)')
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Fig2E_target_correlation_moi1.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, 'Fig2E_target_correlation_moi1.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved target correlation heatmap")

# Plot feature correlation heatmap
fig, ax = plt.subplots(figsize=(10, 10))
corr_plot = B_correlation_features[features_order, :][:, features_order]
im = ax.imshow(corr_plot, cmap='PRGn', vmin=-0.2, vmax=0.2)
cbar = fig.colorbar(im, ax=ax, extend='both', shrink=0.6)
cbar.set_label('Correlation')
ax.set_xlabel('Features')
ax.set_ylabel('Features')
ax.set_title(f'Feature Correlation Matrix ({num_clusters_features} modules)')
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Fig2E_feature_correlation_moi1.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, 'Fig2E_feature_correlation_moi1.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved feature correlation heatmap")

# Create cluster DataFrames
target_clusters_df = pd.DataFrame({
    'target': filtered_covs,
    'cluster': covs_cluster_labels
})

feature_clusters_df = pd.DataFrame({
    'feature': filtered_features,
    'cluster': features_cluster_labels
})

# Print targets in each cluster
print("\n" + "=" * 60)
print("TARGET MODULES")
print("=" * 60)
for i in range(num_clusters_covs):
    targets_in_cluster = target_clusters_df[target_clusters_df['cluster'] == i]['target'].tolist()
    gene_targets = [t for t in targets_in_cluster
                    if not any(x in t for x in ['TENX_', 'G1', 'G2M', 'S', 'UMI', 'cells_per_target'])]
    print(f"\nModule {i} ({len(gene_targets)} genes):")
    print(f"  {', '.join(gene_targets[:20])}{'...' if len(gene_targets) > 20 else ''}")

# Save clustering results
target_clusters_df.to_csv(os.path.join(output_dir, 'target_clusters_moi1.csv'), index=False)
feature_clusters_df.to_csv(os.path.join(output_dir, 'feature_clusters_moi1.csv'), index=False)

# Save beta matrix
beta_df = pd.DataFrame(filtered_B, index=filtered_features, columns=filtered_covs)
beta_df.to_csv(os.path.join(output_dir, 'beta_matrix_moi1.csv'))

# Save numpy arrays
np.save(os.path.join(output_dir, 'beta_matrix_filtered_moi1.npy'), filtered_B)
np.save(os.path.join(output_dir, 'filtered_feature_names_moi1.npy'), filtered_features)
np.save(os.path.join(output_dir, 'filtered_cov_names_moi1.npy'), filtered_covs)
np.save(os.path.join(output_dir, 'feature_cluster_labels_moi1.npy'), features_cluster_labels)
np.save(os.path.join(output_dir, 'target_cluster_labels_moi1.npy'), covs_cluster_labels)

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("FIGURE 2 ANALYSIS SUMMARY (MOI=1 DATA)")
print("=" * 60)
print(f"\nDataset:")
print(f"  Total cells (MOI=1): {adata.shape[0]:,}")
print(f"  HLA_High cells: {(adata.obs['Condition'] == 'High').sum():,}")
print(f"  HLA_Low cells: {(adata.obs['Condition'] == 'Low').sum():,}")
print(f"  Control cells: {(adata.obs['Condition'] == 'Control').sum():,}")

print(f"\nFigure 2A - Cells per Target:")
print(f"  Targets with >{cells_per_target_thresh} cells: {((cells_per_target_high > cells_per_target_thresh) | (cells_per_target_low > cells_per_target_thresh)).sum()}")

print(f"\nFigure 2D - HLA Program Score:")
print(f"  Cohen's d (High vs Low): {cohens_d:.3f}")

print(f"\nFigure 2E - Beta Matrix Clustering:")
print(f"  Filtered features: {len(filtered_features)}")
print(f"  Filtered targets/covariates: {len(filtered_covs)}")
print(f"  Target modules: {num_clusters_covs}")
print(f"  Feature modules: {num_clusters_features}")

print(f"\nOutput files saved to:")
print(f"  {output_dir}")
for f in os.listdir(output_dir):
    print(f"    - {f}")

print(f"\nFigure files saved to:")
print(f"  {fig_dir}")
for f in os.listdir(fig_dir):
    print(f"    - {f}")

print("\n" + "=" * 60)
print("Figure 2 analysis complete!")
print("=" * 60)
