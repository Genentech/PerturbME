import numpy as np
import pandas as pd
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import h5py
import seaborn as sns

full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
fig_pn = os.path.join(full_exp_pn, 'figures/manuscript_figures/Fig2')
save_pn = os.path.join(full_exp_pn, 'adata_arrs')
DE_pn = os.path.join(full_exp_pn, 'DE/DE_csvs')

# Use pre-filtered list of genes from raw h5 to account for genes filtered out due to minimum cell threshold
test_h5_pn = os.path.join(full_exp_pn, 'terra/terra_h5s/CTL_A1a.h5')
gene_names = h5py.File(test_h5_pn , 'r')['matrix']['features']['name'][:].astype(str)

adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))

high_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'High')[0]
low_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Low')[0]
control_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Control')[0]

# Load HLA program from DE
FDR_val=0.05
highE_lowE_df = pd.read_csv(os.path.join(DE_pn, 'RNA_DE_results_MAST_highE_lowE.csv'))
highE_lowE_features = highE_lowE_df['Unnamed: 0'].loc[highE_lowE_df['p_val_adj'] < FDR_val].to_numpy().astype(np.str)
highE_lowE_LFC = highE_lowE_df['avg_logFC'].loc[highE_lowE_df['p_val_adj'] < FDR_val].to_numpy()

upreg_inds = np.where(highE_lowE_LFC > 0)[0]
downreg_inds = np.where(highE_lowE_LFC < 0)[0]

assert(np.all(np.equal(np.union1d(upreg_inds,downreg_inds), np.arange(highE_lowE_LFC.size)))) # Make sure everything up or down regulated
assert(np.intersect1d(upreg_inds, downreg_inds).size==0) # But not both

num_bins=50
HLA_composite_score = np.load(os.path.join(save_pn, 'HLA_composite_score_%d_bins.npy' % num_bins))

# Histogram

fig, ax = plt.subplots(nrows=3)
xlims = [-2, 0.5]
nbins = 200
ax[0].hist(HLA_composite_score[low_inds], bins=nbins, color='#e1812c')
ax[0].set_xlim(xlims), ax[0].set_xticks([])
ax[0].set_ylabel('Low')
ax[1].hist(HLA_composite_score[control_inds], bins=nbins, color='#3274a1')
ax[1].set_ylabel('Control')
ax[1].set_xlim(xlims), ax[1].set_xticks([])
ax[2].hist(HLA_composite_score[high_inds], bins=nbins, color='#3a923a')
ax[2].set_ylabel('High')
ax[2].set_xlim(xlims)
ax[0].set_title('HLA Composite')
fig.savefig(os.path.join(fig_pn, '2D_HLA_composite_by_cond.pdf'), bbox_inches='tight', dpi=800)
plt.close()
