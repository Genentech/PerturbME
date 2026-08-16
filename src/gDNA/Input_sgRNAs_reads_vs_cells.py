import os
from pdb import set_trace as bp
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import scipy.stats as sp_stats
import seaborn as sns
import anndata
import scanpy as sc

full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/20200809_exp/'
figures_pn = os.path.join(full_exp_pn, 'figures')
gDNA_pn = os.path.join(exp_pn, 'guide_seq')
mageck_out_pn = os.path.join(gDNA_pn, 'mageck_out')

adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))
sgRNA_arr = adata_RNA_CITE.obs['sgRNA'].to_numpy().astype(np.str)

# all_samples_count = pd.read_csv(os.path.join(mageck_out_pn, 'all_samples.count.txt'), sep='\t')
all_samples_count = pd.read_csv(os.path.join(mageck_out_pn, 'Input_samples.count.txt'), sep='\t')

all_sgRNAs = all_samples_count['sgRNA'].to_numpy().astype(np.str)
reads_per_sgRNA = np.zeros(all_sgRNAs.size)
cells_per_sgRNA = np.zeros(all_sgRNAs.size)
for sgRNA_i, curr_sgRNA in enumerate(all_sgRNAs):

    reads_per_sgRNA[sgRNA_i] = all_samples_count.loc[all_samples_count['sgRNA']==curr_sgRNA].iloc[:,2:].to_numpy().mean()
    cells_per_sgRNA[sgRNA_i] = np.where(sgRNA_arr == curr_sgRNA)[0].size

bp()

fig, ax = plt.subplots()
ax.scatter(reads_per_sgRNA, cells_per_sgRNA, s=1)
ax.set_xlabel('Reads / sgRNA'), ax.set_ylabel('Cells / sgRNA')
fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_sgRNA.png'), bbox_inches='tight', dpi=800)
plt.close()
