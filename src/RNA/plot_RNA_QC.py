import numpy as np
import pandas as pd
import scanpy as sc
import os,glob
import h5py
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import seaborn as sns
import matplotlib.pyplot as plt

# MIN_CELLS_PER_GENE = 3
# MIN_GENES_PER_CELL = 200
# MAX_PERCENT_MITO = 0.18

# Parse h5 files
def parse_h5(fp):
	f = h5py.File(fp,'r')
	barcodes = f['matrix']['barcodes'][:].astype(str)
	genes = f['matrix']['features']['name'][:].astype(str)
	data = f['matrix']['data']
	indices = f['matrix']['indices']
	indptr = f['matrix']['indptr']
	shape = f['matrix']['shape']
	X = sp_sparse.csr_matrix((data[:],indices[:],indptr[:]), shape=(shape[1], shape[0]))

	return X,barcodes,genes

def main():
	exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	fig_pn = os.path.join(exp_pn, 'figures')
	terra_pn = os.path.join(exp_pn, 'terra')
	h5_dir = os.path.join(terra_pn, 'terra_h5s')

	h5_pns = glob.glob(os.path.join(h5_dir, '*.h5'))
	plot_df_cells = pd.DataFrame()
	# plot_df_perc  = pd.DataFrame()
	plot_df_MT  = pd.DataFrame()
	plot_df_gpc = pd.DataFrame()
	plot_df_cpg = pd.DataFrame()
	plot_df_MT_by_run = pd.DataFrame()
	plot_df_gpc_by_run = pd.DataFrame()
	for pn_i, h5_pn in enumerate(h5_pns):

		channel = h5_pn.split('/')[-1].split('.')[0]
		print('Processing %s' % channel)

		# adata = sc.read_10x_h5(h5_pn)
		X, B, G = parse_h5(h5_pn)
		# Full expression matrix
		curr_expr = np.array(X.todense())

		# Find bad cells and genes
		gene_counts_per_cell = np.diff(X.tocsr().indptr)
		cell_counts_per_gene = np.diff(X.tocsc().indptr)

		# Find mitochondrial %
		MT_gene_inds = np.where(np.char.find(G, 'MT-') == 0)[0] # Gene names starting with 'MT-'
		MT_perc = curr_expr[:, MT_gene_inds].sum(1) / curr_expr.sum(1)

		# Inds to delete
		# delete_gene_inds = np.where(cell_counts_per_gene < MIN_CELLS_PER_GENE)[0]
		# delete_cell_inds = np.where(gene_counts_per_cell < MIN_GENES_PER_CELL)[0]
		# delete_cell_inds = np.union1d(delete_cell_inds, np.where(MT_perc > MAX_PERCENT_MITO)[0])

		# Calculate percent of lost cells
		# lost_cells_perc = 100*(delete_cell_inds.size / B.size)

		# Delete cells and genes
		# B = np.delete(B, delete_cell_inds)
		# G = np.delete(G, delete_gene_inds)
		# curr_expr = np.delete(curr_expr, delete_cell_inds, axis=0)
		# curr_expr = np.delete(curr_expr, delete_gene_inds, axis=1)

		plot_df_cells = plot_df_cells.append(pd.DataFrame({'Sample' : channel, 'Number of Cells' : [B.size]}), ignore_index=True)
		plot_df_MT = plot_df_MT.append(pd.DataFrame({'Sample' : channel, 'MT Percent' : MT_perc}), ignore_index=True)
		# plot_df_perc = plot_df_perc.append(pd.DataFrame({'Sample' : channel, 'Percent Filtered' : [lost_cells_perc]}), ignore_index=True)
		plot_df_gpc = plot_df_gpc.append(pd.DataFrame({'Sample' : channel, 'Genes/Cell' : gene_counts_per_cell}), ignore_index=True)
		plot_df_cpg = plot_df_cpg.append(pd.DataFrame({'Sample' : channel, 'Cells/Gene' : cell_counts_per_gene}), ignore_index=True)

		if '2' in channel:
			plot_df_MT_by_run = plot_df_MT_by_run.append(pd.DataFrame({'Run' : '2', 'MT Percent' : MT_perc}), ignore_index=True)
			plot_df_gpc_by_run = plot_df_gpc_by_run.append(pd.DataFrame({'Run' : '2', 'Genes/Cell' : gene_counts_per_cell}), ignore_index=True)
		elif '1' in channel:
			plot_df_MT_by_run = plot_df_MT_by_run.append(pd.DataFrame({'Run' : '1', 'MT Percent' : MT_perc}), ignore_index=True)
			plot_df_gpc_by_run = plot_df_gpc_by_run.append(pd.DataFrame({'Run' : '1', 'Genes/Cell' : gene_counts_per_cell}), ignore_index=True)
		elif '3' in channel:
			plot_df_MT_by_run = plot_df_MT_by_run.append(pd.DataFrame({'Run' : '3', 'MT Percent' : MT_perc}), ignore_index=True)
			plot_df_gpc_by_run = plot_df_gpc_by_run.append(pd.DataFrame({'Run' : '3', 'Genes/Cell' : gene_counts_per_cell}), ignore_index=True)
		else:
			bp()

	bp()

	plot_df_MT_by_run = plot_df_MT_by_run.sort_values(by='Run')
	fig_save_pn = os.path.join(fig_pn, 'MT_perc_by_run.png')
	fig, ax = plt.subplots(figsize = [1*6.4, 1*4.8])
	sns.violinplot(data=plot_df_MT_by_run, y='MT Percent', x='Run', ax=ax)
	sns.stripplot(data=plot_df_MT_by_run, y='MT Percent', x='Run', ax=ax, size=0.5, color='k')
	ax.axhline(y=0.18, color='r')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	fig_save_pn = os.path.join(fig_pn, 'MT_perc_by_run_hist.png')
	fig, ax = plt.subplots(nrows=3)
	xlims = [0, 0.2]
	nbins=100
	ax[0].hist(plot_df_MT_by_run['MT Percent'].loc[plot_df_MT_by_run['Run'] == '1'], bins=nbins)
	ax[0].set_ylabel('Run #1'), ax[0].set_xlim(xlims), ax[0].set_xticks([])
	ax[1].hist(plot_df_MT_by_run['MT Percent'].loc[plot_df_MT_by_run['Run'] == '2'], bins=nbins)
	ax[1].set_ylabel('Run #2'), ax[1].set_xlim(xlims), ax[1].set_xticks([])
	ax[2].hist(plot_df_MT_by_run['MT Percent'].loc[plot_df_MT_by_run['Run'] == '3'], bins=nbins)
	ax[2].set_ylabel('Run #3'), ax[2].set_xlim(xlims), ax[2].set_xlabel('MT Percent')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	plot_df_gpc_by_run = plot_df_gpc_by_run.sort_values(by='Run')
	fig_save_pn = os.path.join(fig_pn, 'genes_per_cell_by_run.png')
	fig, ax = plt.subplots(figsize = [1*6.4, 1*4.8])
	sns.violinplot(data=plot_df_gpc_by_run, y='Genes/Cell', x='Run', ax=ax)
	sns.stripplot(data=plot_df_gpc_by_run, y='Genes/Cell', x='Run', ax=ax, size=0.5, color='k')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	fig_save_pn = os.path.join(fig_pn, 'genes_per_cell_by_run_hist.png')
	fig, ax = plt.subplots(nrows=3)
	xlims = [0, 7e3]
	nbins=100
	gpc_thresh = 800
	ax[0].hist(plot_df_gpc_by_run['Genes/Cell'].loc[plot_df_gpc_by_run['Run'] == '1'], bins=nbins)
	ax[0].set_ylabel('Run #1'), ax[0].set_xlim(xlims), ax[0].set_xticks([])
	ax[0].axvline(x=gpc_thresh, color='r')
	ax[1].hist(plot_df_gpc_by_run['Genes/Cell'].loc[plot_df_gpc_by_run['Run'] == '2'], bins=nbins)
	ax[1].set_ylabel('Run #2'), ax[1].set_xlim(xlims), ax[1].set_xticks([])
	ax[1].axvline(x=gpc_thresh, color='r')
	ax[2].hist(plot_df_gpc_by_run['Genes/Cell'].loc[plot_df_gpc_by_run['Run'] == '3'], bins=nbins)
	ax[2].set_ylabel('Run #3'), ax[2].set_xlim(xlims), ax[2].set_xlabel('Genes / Cell')
	ax[2].axvline(x=gpc_thresh, color='r')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()


	plot_df_cells = plot_df_cells.sort_values(by='Sample')
	fig_save_pn = os.path.join(fig_pn, 'unfiltered_cell_count.png')
	fig, ax = plt.subplots(figsize = [1*6.4, 1*4.8])
	sns.barplot(data=plot_df_cells, x='Sample', y='Number of Cells', ax=ax, color='b')
	plt.xticks(rotation=90)
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	plot_df_MT = plot_df_MT.sort_values(by='Sample')
	fig_save_pn = os.path.join(fig_pn, 'MT_perc.png')
	fig, ax = plt.subplots(figsize = [1*6.4, 1*4.8])
	sns.violinplot(data=plot_df_MT, x='MT Percent', y='Sample', ax=ax)
	# sns.stripplot(data=plot_df_MT, x='MT Percent', y='Sample', ax=ax, size=1)
	ax.axvline(x=0.15, color='r')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	#
	# plot_df_perc = plot_df_perc.sort_values(by='Sample')
	# fig_save_pn = os.path.join(fig_pn, 'filtered_cell_percent.png')
	#
	# fig, ax = plt.subplots(figsize = [2*6.4, 2*4.8])
	# sns.barplot(data=plot_df_perc, x='Sample', y='Percent Filtered', ax=ax)
	# fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	# plt.close()

	plot_df_gpc = plot_df_gpc.sort_values(by='Sample')
	fig_save_pn = os.path.join(fig_pn, 'genes_per_cell.png')
	fig, ax = plt.subplots(figsize = [1*6.4, 1*4.8])
	sns.violinplot(data=plot_df_gpc, x='Genes/Cell', y='Sample', ax=ax)
	# sns.stripplot(data=plot_df_gpc, x='Genes/Cell', y='Sample', ax=ax, size=2)
	genes_per_cell_thresh = 800
	ax.axvline(x=genes_per_cell_thresh, color='r')
	ax.axvline(x=6000, color='r')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	HLA_high_inds = np.where(plot_df_gpc['Sample'] == 'HLA_high')[0]
	print(100*(np.where(plot_df_gpc['Genes/Cell'].to_numpy()[HLA_high_inds] < genes_per_cell_thresh)[0].size/HLA_high_inds.size))
	HLA_low_inds = np.where(plot_df_gpc['Sample'] == 'HLA_low')[0]
	print(100*(np.where(plot_df_gpc['Genes/Cell'].to_numpy()[HLA_low_inds] < genes_per_cell_thresh)[0].size/HLA_low_inds.size))
	Input_inds = np.where(plot_df_gpc['Sample'] == 'Input')[0]
	print(100*(np.where(plot_df_gpc['Genes/Cell'].to_numpy()[Input_inds] < genes_per_cell_thresh)[0].size/Input_inds.size))

	# fig_save_pn = os.path.join(fig_pn, 'genes_per_cell.png')
	# fig, ax = plt.subplots()
	# ax.hist(genes_per_cell, bins=3000)
	# ax.set_xlabel('Genes / Cell'), ax.set_ylabel('Count')
	# # ax.axvline(MIN_GENES_PER_CELL, color='r', linewidth=1)
	# fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	# plt.close()

	plot_df_cpg = plot_df_cpg.sort_values(by='Sample')
	plot_df_cpg_log = plot_df_cpg.copy()
	plot_df_cpg_log['Cells/Gene'] = np.log10(plot_df_cpg['Cells/Gene'] + 1)
	fig_save_pn = os.path.join(fig_pn, 'cells_per_gene_log.png')
	fig, ax = plt.subplots(figsize = [1*6.4, 1*4.8])
	sns.violinplot(data=plot_df_cpg_log, x='Cells/Gene', y='Sample', ax=ax)
	# sns.stripplot(data=plot_df_cpg, x='Cells/Gene', y='Sample', ax=ax, size=2)
	cells_per_gene_thresh = np.log10(10)
	ax.axvline(x=cells_per_gene_thresh, color='r')
	ax.set_xlabel('Log(Cells/Gene)')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	# fig_save_pn = os.path.join(fig_pn, 'cells_per_gene.png')
	# fig, ax = plt.subplots()
	# ax.hist(cells_per_gene, bins=6000)
	# ax.set_xlabel('Cells / Gene'), ax.set_ylabel('Count')
	# # ax.axvline(MIN_GENES_PER_CELL, color='r', linewidth=1)
	# fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	# plt.close()

if __name__ == "__main__":
    main()
