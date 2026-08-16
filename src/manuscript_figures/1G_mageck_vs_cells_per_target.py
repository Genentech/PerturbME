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
from tqdm import tqdm
from adjustText import adjust_text

full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
figures_pn = os.path.join(full_exp_pn, 'figures/manuscript_figures/Fig1')
gDNA_pn = os.path.join(exp_pn, 'guide_seq')
mageck_out_pn = os.path.join(gDNA_pn, 'mageck_out')

adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))
sgRNA_arr = adata_RNA_CITE.obs['sgRNA'].to_numpy().astype(np.str)
high_inds = np.where(adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str) == 'High')[0]
low_inds = np.where(adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str) == 'Low')[0]

# all_samples_count = pd.read_csv(os.path.join(mageck_out_pn, 'all_samples.count.txt'), sep='\t')
all_samples_count = pd.read_csv(os.path.join(mageck_out_pn, 'all_samples.count.txt'), sep='\t')
all_sgRNAs = all_samples_count['sgRNA'].to_numpy().astype(np.str)
all_genes = np.unique(all_samples_count['Gene'].to_numpy().astype(np.str))

high_counts = all_samples_count['High_A'].to_numpy()+all_samples_count['High_B'].to_numpy()
low_counts = all_samples_count['Low_A'].to_numpy()+all_samples_count['Low_B'].to_numpy()
assert(high_counts.size == low_counts.size)

reads_per_gene_high = np.load(os.path.join(mageck_out_pn, 'reads_per_gene_high.npy'))
reads_per_gene_low = np.load(os.path.join(mageck_out_pn, 'reads_per_gene_low.npy'))
cells_per_gene_high = np.load(os.path.join(mageck_out_pn, 'cells_per_gene_high.npy'))
cells_per_gene_low = np.load(os.path.join(mageck_out_pn, 'cells_per_gene_low.npy'))
all_genes = np.load(os.path.join(mageck_out_pn, 'all_genes.npy'))

# save_mat = np.column_stack((reads_per_gene_high, cells_per_gene_high, reads_per_gene_low, cells_per_gene_low))
# mat_labels = ['Reads / Gene High', 'Cells / Gene High', 'Reads / Gene Low', 'Cells / Gene Low']
# pd.DataFrame(data=save_mat, index=all_genes, columns=mat_labels).to_csv(os.path.join(figures_pn, 'gDNA/reads_vs_cells.csv'))

control_inds = np.where( (all_genes=='NO_SITE') | (all_genes=='ONE_NON-GENE_SITE'))[0]

fig, ax = plt.subplots()
ax.scatter(reads_per_gene_high, cells_per_gene_high, s=1)
ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
ax.set_title('HLA High')
ax.set_xlim([0, 150e3]), ax.set_ylim([0, 80])
label_inds = np.where(cells_per_gene_high > 25)[0]
label_inds = np.union1d(label_inds, np.where( (reads_per_gene_high > 40e3) & (cells_per_gene_high > 25) )[0])
label_inds = np.delete(label_inds, np.where(np.isin(label_inds, control_inds))[0])
text_objs = [ax.text(reads_per_gene_high[ind], cells_per_gene_high[ind], all_genes[ind], fontsize=10, fontweight='bold') for ind in label_inds]
adjust_text(text_objs, arrowprops=dict(arrowstyle='-', color='red'))
fig.savefig(os.path.join(figures_pn, '1G_reads_vs_cells_per_target_high.pdf'), bbox_inches='tight', dpi=800)
plt.close()
fig, ax = plt.subplots()
ax.scatter(reads_per_gene_low, cells_per_gene_low, s=1)
ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
ax.set_title('HLA Low')
ax.set_xlim([0, 3e5]), ax.set_ylim([0, 200])
label_inds = np.where(cells_per_gene_low > 40)[0]
label_inds = np.union1d(label_inds, np.where( (reads_per_gene_low > 75e3) & (cells_per_gene_low > 40) )[0])
label_inds = np.delete(label_inds, np.where(np.isin(label_inds, control_inds))[0])
text_objs = [ax.text(reads_per_gene_low[ind], cells_per_gene_low[ind], all_genes[ind], fontsize=10, fontweight='bold') for ind in label_inds]
adjust_text(text_objs, arrowprops=dict(arrowstyle='-', color='red'))
fig.savefig(os.path.join(figures_pn, '1G_reads_vs_cells_per_target_low.pdf'), bbox_inches='tight', dpi=800)
plt.close()
