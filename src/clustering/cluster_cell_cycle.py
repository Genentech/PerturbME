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
save_pn = os.path.join(full_exp_pn, 'adata_arrs')

adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))

high_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'High')[0]
low_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Low')[0]
control_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Control')[0]

cycle_df = pd.read_excel('/ahg/regevdata/users/cfrangie/Projects/01_Perturb_Melanoma/full_exp/signatures/Tirosh_cycles.xlsx')

G1S_genes = cycle_df['G1/S'].to_numpy().astype(np.str)
G1S_genes = G1S_genes[~(G1S_genes == 'nan')]
G1S_genes_filt = [gene.strip() for gene in G1S_genes if gene.strip() in adata_RNA_CITE.var_names]

G2M_genes = cycle_df['G2/M'].to_numpy().astype(np.str)
G2M_genes = G2M_genes[~(G2M_genes == 'nan')]
G2M_genes_filt = [gene.strip() for gene in G2M_genes if gene.strip() in adata_RNA_CITE.var_names]

# Score cell cycle genes
sc.tl.score_genes_cell_cycle(adata_RNA_CITE, s_genes=G1S_genes_filt, g2m_genes=G2M_genes_filt)
np.save(os.path.join(save_pn, 'cell_cycle_Tirosh.npy'), adata_RNA_CITE.obs['phase'].to_numpy().astype(np.str))


fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[high_inds,:], color=['phase'], ax=ax)
ax.set_title('Cell Cycle - High')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_cell_cycle_high.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[low_inds,:], color=['phase'], ax=ax)
ax.set_title('Cell Cycle - Low')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_cell_cycle_low.png'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
sc.pl.umap(adata_RNA_CITE[control_inds,:], color=['phase'], ax=ax)
ax.set_title('Cell Cycle - CTL')
fig.savefig(os.path.join(fig_pn, 'Clustering/UMAP_RNA_CITE_cell_cycle_control.png'), bbox_inches='tight', dpi=800)
plt.close()
