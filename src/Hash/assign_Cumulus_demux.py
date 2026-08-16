import os,glob
from pdb import set_trace as bp
import numpy as np
import pandas as pd
import seaborn as sns
import anndata
import h5py
import matplotlib.pyplot as plt
import scipy.stats as sp_stats
# Conda environment anndat

exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
fig_pn = os.path.join(exp_pn, 'figures')
terra_pn = os.path.join(exp_pn, 'terra')
barcode_txt_pn = os.path.join(terra_pn, 'Barcodes')
terra_Hash_path = os.path.join(exp_pn, 'Hash/cumulus_demux')

demux_FPs = glob.glob(os.path.join(terra_Hash_path, '*_demux.h5ad'))
plot_df = pd.DataFrame([])
plot_df_perc = pd.DataFrame([])
for demux_FP in demux_FPs:

	channel = demux_FP.split('/')[-1].split('demux')[0][:-1]

	print('Processing channel %s' % channel)

	curr_adata = anndata.read_h5ad(demux_FP)

	RNA_barcodes = np.loadtxt(os.path.join(barcode_txt_pn, channel+'.txt'), dtype=str)
	adata_barcodes = curr_adata.obs.index.to_numpy().astype(np.str)

	good_adata_barcodes = np.where(np.isin(adata_barcodes, RNA_barcodes))[0]
	# assert(good_adata_barcodes.size == RNA_barcodes.size)

	HTO_assignment_type = curr_adata.obs['demux_type'].to_numpy().astype(np.str)[good_adata_barcodes]
	HTO_assignment = curr_adata.obs['assignment'].to_numpy().astype(np.str)[good_adata_barcodes]
	Cumulus_barcodes = adata_barcodes[good_adata_barcodes]

	np.save(os.path.join(terra_Hash_path, '%s.Cumulus_barcodes.npy' % channel), Cumulus_barcodes.astype(np.str))
	np.save(os.path.join(terra_Hash_path, '%s.Cumulus_assignment_type.npy' % channel), HTO_assignment_type.astype(np.str))
	np.save(os.path.join(terra_Hash_path, '%s.Cumulus_assignment.npy' % channel), HTO_assignment.astype(np.str))

	singlet_count = np.where(HTO_assignment_type == 'singlet')[0].size
	doublet_count = np.where(HTO_assignment_type == 'doublet')[0].size
	unknown_count = np.where(HTO_assignment_type == 'unknown')[0].size

	plot_df = plot_df.append(pd.DataFrame({'Channel' : channel, 'Count' : [singlet_count], 'Type' : 'Singlet'}), ignore_index=True)
	plot_df = plot_df.append(pd.DataFrame({'Channel' : channel, 'Count' : [doublet_count], 'Type' : 'Doublet'}), ignore_index=True)
	plot_df = plot_df.append(pd.DataFrame({'Channel' : channel, 'Count' : [unknown_count], 'Type' : 'Unkown'}), ignore_index=True)

	plot_df_perc = plot_df_perc.append(pd.DataFrame({'Channel' : channel, 'Count' : [100*(singlet_count/good_adata_barcodes.size)], 'Type' : 'Singlet'}), ignore_index=True)
	plot_df_perc = plot_df_perc.append(pd.DataFrame({'Channel' : channel, 'Count' : [100*(doublet_count/good_adata_barcodes.size)], 'Type' : 'Doublet'}), ignore_index=True)
	plot_df_perc = plot_df_perc.append(pd.DataFrame({'Channel' : channel, 'Count' : [100*(unknown_count/good_adata_barcodes.size)], 'Type' : 'Unkown'}), ignore_index=True)

bp()
plot_df = plot_df.sort_values(by='Channel')
save_fig_pn = os.path.join(fig_pn, 'Cumulus_assignments.png')
fig, ax = plt.subplots()
sns.barplot(data=plot_df, x='Channel', y='Count', hue='Type', ax=ax)
ax.set_xlabel(''), ax.set_ylabel('Number Cells')
plt.xticks(rotation=90)
fig.savefig(save_fig_pn, dpi=800, bbox_inches='tight')
plt.close()

plot_df_perc = plot_df_perc.sort_values(by='Channel')
save_fig_pn = os.path.join(fig_pn, 'Cumulus_assignments_percent.png')
fig, ax = plt.subplots()
sns.barplot(data=plot_df_perc, x='Channel', y='Count', hue='Type', ax=ax)
ax.set_xlabel(''), ax.set_ylabel('Percent Total Cells')
plt.xticks(rotation=90)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig.savefig(save_fig_pn, dpi=800, bbox_inches='tight')
plt.close()
