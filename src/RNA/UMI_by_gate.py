import os
from pdb import set_trace as bp
import numpy as np
import pandas as pd
import glob
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as sp_stats
import seaborn as sns
import anndata
import scanpy as sc

full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
figures_pn = os.path.join(full_exp_pn, 'figures')

adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))

bp()

UMI_arr = np.load(os.path.join(exp_pn, 'adata_arrs/UMI_count.npy'))

high_inds = np.where(adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str) == 'High')[0]
low_inds = np.where(adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str) == 'Low')[0]
ctl_inds = np.where(adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str) == 'Control')[0]

plot_df = pd.DataFrame({'UMI' : UMI_arr, 'Condition' : adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str)})

fig, ax = plt.subplots()
sns.violinplot(data=plot_df, x='Condition', y='UMI')
ax.set_xlabel(''), ax.set_ylabel('UMI Count')
fig.savefig(os.path.join(figures_pn, 'UMI_count_by_condition.png'), bbox_inches='tight', dpi=800)
plt.close()

num_bins = 200
xlims = [0, 18e3]
fig, ax = plt.subplots(nrows=3)
ax[0].hist(UMI_arr[low_inds], bins=num_bins)
ax[0].set_xlim(xlims), ax[0].set_xticks([])
ax[0].set_ylabel('Low')
ax[1].hist(UMI_arr[ctl_inds], bins=num_bins)
ax[1].set_xlim(xlims), ax[1].set_xticks([])
ax[1].set_ylabel('Control')
ax[2].hist(UMI_arr[high_inds], bins=num_bins)
ax[2].set_xlim(xlims), ax[2].set_xlabel('UMI Count')
ax[2].set_ylabel('High')
fig.savefig(os.path.join(figures_pn, 'UMI_count_by_condition_hist.png'), bbox_inches='tight', dpi=800)
plt.close()
