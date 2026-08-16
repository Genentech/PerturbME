import numpy as np
import pandas as pd
import seaborn as sns
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt

READS_THRESH = 1
READS_PERC_THRESH = 200
# PERC_THRESH = 0.2
PERC_THRESH = 2 # PERC THRESH OFF
UMI_COUNT_THRESH = 10e3
UMI_READS_THRESH = 10e3
FILTERED_UMI_THRESH = 1

exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
figures_pn = os.path.join(exp_pn, 'figures')
sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
CROP_dir = os.path.join(exp_pn, 'CROP')
design_dir = os.path.join(CROP_dir, 'design_mats/')
filtered_RNA_pn = os.path.join(exp_pn, 'filtered_RNA')

# Load sgRNA information
gw_guides_pn = os.path.join(sequencing_info_pn, 'gw_lib_guides.csv')
target_sgRNA_df = pd.read_csv(gw_guides_pn)
target_sgRNA_barcodes = target_sgRNA_df['sgRNA_barcode'].to_numpy().astype(np.str)
target_sgRNA_names =  target_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
assert(target_sgRNA_barcodes.size == target_sgRNA_names.size)
control_sgRNA_df = pd.read_csv(os.path.join(sequencing_info_pn, 'human_TF_CTL_guides.txt'), sep='\t', header=None, names=['sgRNA_name', 'sgRNA_barcode'])
control_sgRNA_barcodes = control_sgRNA_df['sgRNA_barcode'].to_numpy().astype(np.str)
control_sgRNA_names =  control_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
assert(control_sgRNA_barcodes.size == control_sgRNA_names.size)

all_design_mats = glob.glob(os.path.join(design_dir, '*.sgRNA_design_mat.%d_%d_%.2f_%d_%d_%d.npy' % (READS_THRESH,READS_PERC_THRESH,PERC_THRESH,UMI_COUNT_THRESH,UMI_READS_THRESH,FILTERED_UMI_THRESH)))

HLA_low_dict = {}
HLA_high_dict = {}
for curr_mat in all_design_mats:

	cond = curr_mat.split('/')[-1].split('.')[0]
	print('Processing', cond)

	if 'CTL' not in cond:

		design_mat = np.load(curr_mat)

		MOI = design_mat.sum(1)
		MOI2_cells = np.where(MOI == 2)[0]

		for MOI2_cell in MOI2_cells:
			doublet_sgRNA_inds = np.where(design_mat[MOI2_cell, :].flatten() != 0)[0]
			assert(doublet_sgRNA_inds.size == 2)

			sgRNA_names = np.sort(target_sgRNA_names[doublet_sgRNA_inds]) # Sort so dictionary key is consistent
			assert(sgRNA_names.size == 2)

			dict_key = sgRNA_names[0].split('_')[0]+'-'+sgRNA_names[1].split('_')[0]

			if 'High' in cond:
				if dict_key in HLA_high_dict: # Increment count for this combination
					HLA_high_dict[dict_key] += 1
					print('High found %d combinations for %s' % (HLA_high_dict[dict_key], dict_key))
				else: # Add to dictionary
					HLA_high_dict[dict_key] = 1
			elif 'Low' in cond:
				if dict_key in HLA_low_dict: # Increment count for this combination
					HLA_low_dict[dict_key] += 1
					print('Low found %d combinations for %s' % (HLA_low_dict[dict_key], dict_key))
				else: # Add to dictionary
					HLA_low_dict[dict_key] = 1
			else:
				bp()

bp()

cells_per_combo_low = np.fromiter(HLA_low_dict.values(), dtype=np.int)
cells_per_combo_high = np.fromiter(HLA_high_dict.values(), dtype=np.int)

sgRNAs_to_plot = 100
xtick_fontsize = 4
save_fig_pn = os.path.join(figures_pn, 'CROP/Low_cells_per_sgRNA_sorted_%d_sgRNAs.png' % sgRNAs_to_plot)
fig, ax = plt.subplots()
sgRNA_plot_inds = np.flip(np.argsort(HLA_low_cpg))[:sgRNAs_to_plot]
ax.bar(np.arange(sgRNAs_to_plot), HLA_low_cpg[sgRNA_plot_inds])
ax.set_xticks(np.arange(sgRNAs_to_plot)), ax.set_xticklabels(target_sgRNA_names[sgRNA_plot_inds])
ax.set_ylabel('# Cells'), plt.xticks(rotation=90, fontsize=xtick_fontsize)
ax.set_title('HLA Low')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()
save_fig_pn = os.path.join(figures_pn, 'CROP/High_cells_per_sgRNA_sorted_%d_sgRNAs.png' % sgRNAs_to_plot)
fig, ax = plt.subplots()
sgRNA_plot_inds = np.flip(np.argsort(HLA_high_cpg))[:sgRNAs_to_plot]
ax.bar(np.arange(sgRNAs_to_plot), HLA_high_cpg[sgRNA_plot_inds])
ax.set_xticks(np.arange(sgRNAs_to_plot)), ax.set_xticklabels(target_sgRNA_names[sgRNA_plot_inds])
ax.set_ylabel('# Cells'), plt.xticks(rotation=90, fontsize=xtick_fontsize)
ax.set_title('HLA High')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()
save_fig_pn = os.path.join(figures_pn, 'CROP/CTL_cells_per_sgRNA_sorted_%d_sgRNAs.png' % sgRNAs_to_plot)
fig, ax = plt.subplots()
sgRNA_plot_inds = np.flip(np.argsort(CTL_cpg))[:sgRNAs_to_plot]
ax.bar(np.arange(sgRNAs_to_plot), CTL_cpg[sgRNA_plot_inds])
ax.set_xticks(np.arange(sgRNAs_to_plot)), ax.set_xticklabels(control_sgRNA_names[sgRNA_plot_inds])
ax.set_ylabel('# Cells'), plt.xticks(rotation=90, fontsize=xtick_fontsize)
ax.set_title('CTL')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()

# Targets with > 1 sgRNA with > 10 cells
np.sort(target_sgRNA_names[np.where(HLA_low_cpg > 10)[0]])
np.sort(target_sgRNA_names[np.where(HLA_high_cpg > 10)[0]])
