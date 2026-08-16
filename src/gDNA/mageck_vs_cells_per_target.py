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
figures_pn = os.path.join(full_exp_pn, 'figures')
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

if False:

    reads_per_gene_high = np.zeros(all_genes.size)
    reads_per_gene_low = np.zeros(all_genes.size)
    cells_per_gene_high = np.zeros(all_genes.size)
    cells_per_gene_low = np.zeros(all_genes.size)
    for gene_i, curr_gene in tqdm(enumerate(all_genes)):

        mageck_inds = np.where(np.char.startswith(all_sgRNAs, curr_gene+'_'))[0]
        if mageck_inds.size != 3:
            print(curr_gene)

        reads_per_gene_high[gene_i] = high_counts[mageck_inds].sum()
        reads_per_gene_low[gene_i] = low_counts[mageck_inds].sum()

        cells_with_target = np.where(np.char.startswith(sgRNA_arr, curr_gene+'_'))[0]

        cells_per_gene_high[gene_i] = np.intersect1d(cells_with_target, high_inds).size
        cells_per_gene_low[gene_i] = np.intersect1d(cells_with_target, low_inds).size

    np.save(os.path.join(mageck_out_pn, 'reads_per_gene_high.npy'), reads_per_gene_high)
    np.save(os.path.join(mageck_out_pn, 'reads_per_gene_low.npy'), reads_per_gene_low)
    np.save(os.path.join(mageck_out_pn, 'cells_per_gene_high.npy'), cells_per_gene_high)
    np.save(os.path.join(mageck_out_pn, 'cells_per_gene_low.npy'), cells_per_gene_low)
    np.save(os.path.join(mageck_out_pn, 'all_genes.npy'), all_genes)

else:

    reads_per_gene_high = np.load(os.path.join(mageck_out_pn, 'reads_per_gene_high.npy'))
    reads_per_gene_low = np.load(os.path.join(mageck_out_pn, 'reads_per_gene_low.npy'))
    cells_per_gene_high = np.load(os.path.join(mageck_out_pn, 'cells_per_gene_high.npy'))
    cells_per_gene_low = np.load(os.path.join(mageck_out_pn, 'cells_per_gene_low.npy'))
    all_genes = np.load(os.path.join(mageck_out_pn, 'all_genes.npy'))

bp()
save_mat = np.column_stack((reads_per_gene_high, cells_per_gene_high, reads_per_gene_low, cells_per_gene_low))
mat_labels = ['Reads / Gene High', 'Cells / Gene High', 'Reads / Gene Low', 'Cells / Gene Low']
pd.DataFrame(data=save_mat, index=all_genes, columns=mat_labels).to_csv(os.path.join(figures_pn, 'gDNA/reads_vs_cells.csv'))

control_inds = np.where( (all_genes=='NO_SITE') | (all_genes=='ONE_NON-GENE_SITE'))[0]

cmap =  matplotlib.cm.get_cmap('Reds')
cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=20)
color_arr = cmap(cmap_norm(cells_per_gene_high))
fig, ax = plt.subplots()
ax.scatter(reads_per_gene_high, reads_per_gene_low, s=1, c=color_arr)
ax.set_xlabel('Reads / Target High'), ax.set_ylabel('Read / Target Low')
ax.set_title('HLA High')
ax.set_xlim([-1000, 150e3]), ax.set_ylim([-2000, 3e5])
fig.savefig(os.path.join(figures_pn, 'reads_per_target_high_vs_low_high.png'), bbox_inches='tight', dpi=800)
plt.close()
color_arr = cmap(cmap_norm(cells_per_gene_low))
fig, ax = plt.subplots()
ax.scatter(reads_per_gene_high, reads_per_gene_low, s=1, c=color_arr)
ax.set_xlabel('Reads / Target High'), ax.set_ylabel('Reads / Target Low')
ax.set_title('HLA Low')
ax.set_xlim([-1000, 150e3]), ax.set_ylim([-2000, 3e5])
fig.savefig(os.path.join(figures_pn, 'reads_per_target_high_vs_low_low.png'), bbox_inches='tight', dpi=800)
plt.close()

fig, ax = plt.subplots(figsize=(1, 12))
fig.subplots_adjust(bottom=0.5)
norm = matplotlib.colors.Normalize(vmin=0, vmax=20)
cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=plt.cm.Reds, norm=norm, orientation='vertical')
cb1.set_label('Cells / Target', fontsize=16, fontweight='bold')
fig.savefig(os.path.join(figures_pn, 'reads_per_target_high_vs_low_cmap.png'), bbox_inches='tight', dpi=800)
plt.close()

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
fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_high.png'), bbox_inches='tight', dpi=800)
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
fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_low.png'), bbox_inches='tight', dpi=800)
plt.close()
