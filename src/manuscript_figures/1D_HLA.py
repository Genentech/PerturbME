import numpy as np
import pandas as pd
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt

full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
fig_pn = os.path.join(full_exp_pn, 'figures/manuscript_figures/Fig1')

adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))

high_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'High')[0]
low_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Low')[0]
control_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Control')[0]

condition_arr = adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str)

low_inds = np.where(condition_arr == 'Low')[0]
high_inds = np.where(condition_arr == 'High')[0]
control_inds = np.where(condition_arr == 'Control')[0]

# Shade individual features
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE, color=['CITE-HLA_A'], ax=ax, color_map='hsv', vmin=2, vmax=6)
ax.set_title('HLA-A,B,C Protein')
fig.savefig(os.path.join(fig_pn, '1D_UMAP_RNA_CITE_HLA_A.pdf'), bbox_inches='tight', dpi=800)
plt.close()
