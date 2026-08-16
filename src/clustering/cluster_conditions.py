import numpy as np
import pandas as pd
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt

full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
fig_pn = os.path.join(full_exp_pn, 'figures')

adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))

high_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'High')[0]
low_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Low')[0]
control_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Control')[0]

bp()

cell_cycle_arr = np.load(os.path.join(full_exp_pn, 'adata_arrs/cell_cycle_Tirosh.npy'))
UMI_arr = np.load(os.path.join(full_exp_pn, 'adata_arrs/UMI_count.npy'))

adata_RNA_CITE.obs['Cell Cycle'] = cell_cycle_arr
adata_RNA_CITE.obs['UMI'] = UMI_arr

fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE, color=['Cell Cycle'], ax=ax)
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_cell_cycle.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE, color=['UMI'], ax=ax, color_map='hsv')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_UMI.png'), bbox_inches='tight', dpi=800)
plt.close()

fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE, color=['Condition'], ax=ax)
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_conditions.png'), bbox_inches='tight', dpi=800)
plt.close()

condition_arr = adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str)

low_inds = np.where(condition_arr == 'Low')[0]
high_inds = np.where(condition_arr == 'High')[0]
control_inds = np.where(condition_arr == 'Control')[0]
dot_size = 1
dot_alpha = 1
label_fontsize = 20
color_arr = ['#e1812c', '#3a923a', '#3274a1']
fig, ax = plt.subplots(figsize=[2*6.4, 2*4.8])
ax.scatter(adata_RNA_CITE.obsm['X_umap'][low_inds, 0], adata_RNA_CITE.obsm['X_umap'][low_inds, 1], edgecolor=None, s=dot_size, alpha=dot_alpha, c=color_arr[0], label='Low')
ax.scatter(adata_RNA_CITE.obsm['X_umap'][high_inds, 0], adata_RNA_CITE.obsm['X_umap'][high_inds, 1], edgecolor=None, s=dot_size, alpha=dot_alpha, c=color_arr[1], label='High')
ax.scatter(adata_RNA_CITE.obsm['X_umap'][control_inds, 0], adata_RNA_CITE.obsm['X_umap'][control_inds, 1], edgecolor=None, s=dot_size, alpha=dot_alpha, c=color_arr[2], label='Control')
ax.set_xlabel('UMAP1', fontsize=label_fontsize, fontname='Helvetica'), ax.set_ylabel('UMAP2', fontsize=label_fontsize, fontname='Helvetica')
ax.set_xticks([]), ax.set_yticks([])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale=10)
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_conditions.png'), bbox_inches='tight', dpi=200)
plt.close()

fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[high_inds, :], color=['Condition'], ax=ax)
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_high_condition.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[low_inds, :], color=['Condition'], ax=ax)
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_low_condition.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[control_inds, :], color=['Condition'], ax=ax)
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_control_condition.png'), bbox_inches='tight', dpi=800)
plt.close()

fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[high_inds,:], color=['Sample'], ax=ax)
ax.set_title('Channel - High')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_channel_high.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[low_inds,:], color=['Sample'], ax=ax)
ax.set_title('Channel - Low')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_channel_low.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[control_inds,:], color=['Sample'], ax=ax)
ax.set_title('Channel - CTL')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_channel_CTL.png'), bbox_inches='tight', dpi=800)
plt.close()

# Shade individual features
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE, color=['CITE-HLA_A'], ax=ax, color_map='hsv', vmin=2, vmax=6)
ax.set_title('HLA-A,B,C Protein')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_HLA_A.png'), bbox_inches='tight', dpi=800)
plt.close()

# Annotate 10X Run
channel_arr = adata_RNA_CITE.obs['Sample'].to_numpy().astype(np.str)
run1_inds = np.where(np.char.find(channel_arr, '1') != -1)[0]
run2_inds = np.where(np.char.find(channel_arr, '2') != -1)[0]
run12_inds = np.intersect1d(run1_inds, run2_inds)
run3_inds = np.where(np.char.find(channel_arr, '3') != -1)[0]

run_arr = np.full(channel_arr.size, '00')
run_arr[run1_inds] = '1'
run_arr[run2_inds] = '2'
run_arr[run12_inds] = '12'
run_arr[run3_inds] = '3'
assert(np.where(run_arr == '00')[0].size==0)
adata_RNA_CITE.obs['10X Run'] = pd.Categorical(run_arr)

fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[high_inds,:], color=['10X Run'], ax=ax)
ax.set_title('10X Run - High')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_10X_run_high.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[low_inds,:], color=['10X Run'], ax=ax)
ax.set_title('10X Run - Low')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_10X_run_low.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[control_inds,:], color=['10X Run'], ax=ax)
ax.set_title('10X Run - CTL')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_10X_run_control.png'), bbox_inches='tight', dpi=800)
plt.close()

# Annotate sorter
sorterA_inds = np.where(np.char.find(channel_arr, 'A') != -1)[0]
sorterB_inds = np.where(np.char.find(channel_arr, 'B') != -1)[0]
sorter_arr = np.full(channel_arr.size, '0')
sorter_arr[sorterA_inds] = 'A'
sorter_arr[sorterB_inds] = 'B'
assert(np.where(sorter_arr == '0')[0].size==0)
adata_RNA_CITE.obs['Sorter'] = pd.Categorical(sorter_arr)

fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[high_inds,:], color=['Sorter'], ax=ax)
ax.set_title('Sorter - High')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_sorter_high.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[low_inds,:], color=['Sorter'], ax=ax)
ax.set_title('Sorter - Low')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_sorter_low.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[control_inds,:], color=['Sorter'], ax=ax)
ax.set_title('Sorter - CTL')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_sorter_control.png'), bbox_inches='tight', dpi=800)
plt.close()

# Annotate UMI counts
UMI_counts = np.load(os.path.join(full_exp_pn, 'adata_arrs/UMI_count.npy'))
adata_RNA_CITE.obs['UMI'] = UMI_counts

fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[high_inds,:], color=['UMI'], ax=ax, color_map='hsv')
ax.set_title('UMI - High')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_UMI_high.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[low_inds,:], color=['UMI'], ax=ax, color_map='hsv')
ax.set_title('UMI - Low')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_UMI_low.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[control_inds,:], color=['UMI'], ax=ax, color_map='hsv')
ax.set_title('UMI - CTL')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_UMI_control.png'), bbox_inches='tight', dpi=800)
plt.close()
