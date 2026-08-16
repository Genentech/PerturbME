import numpy as np
import pandas as pd
import scanpy as sc
import os,glob
import h5py
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import seaborn as sns
import matplotlib.pyplot as plt

# MIN_CELLS_PER_GENE = 10
MIN_GENES_PER_CELL = 800
MAX_GENES_PER_CELL = 6000
MAX_PERCENT_MITO = 0.15

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
	terra_Hash_path = os.path.join(exp_pn, 'Hash/cumulus_demux')
	save_pn = os.path.join(exp_pn, 'filtered_RNA')

	# Load filtered gene list
	# filtered_genes = np.load(os.path.join(save_pn, 'filtered_genes_%d_cells.npy' % MIN_CELLS_PER_GENE))
	# assert(filtered_genes.size == np.unique(filtered_genes).size)

	h5_pns = glob.glob(os.path.join(h5_dir, '*.h5'))

	plot_df_cells = pd.DataFrame()
	plot_df_perc  = pd.DataFrame()
	for pn_i, h5_pn in enumerate(h5_pns):

		channel = h5_pn.split('/')[-1].split('.')[0]
		print('Processing %s' % channel)

		# adata = sc.read_10x_h5(h5_pn)
		X, B, G = parse_h5(h5_pn)

		if pn_i == 0:
			filtered_genes = np.unique(G)
			np.save(os.path.join(save_pn, 'all_genes.npy'), filtered_genes.astype(np.str))
		else:
			assert(np.all(np.char.equal(filtered_genes, np.unique(G))))

		# Full expression matrix
		curr_expr = np.array(X.todense())

		# Find # genes in each cell
		gene_counts_per_cell = np.diff(X.tocsr().indptr)

		# Find mitochondrial % in each cell
		MT_gene_inds = np.where(np.char.find(G, 'MT-') == 0)[0] # Gene names starting with 'MT-'
		MT_perc = curr_expr[:, MT_gene_inds].sum(1) / curr_expr.sum(1)

		# Get cumulus demux information
		HTO_assignment = np.load(os.path.join(terra_Hash_path, '%s.Cumulus_assignment_type.npy' % channel))
		HTO_barcodes = np.char.add(np.load(os.path.join(terra_Hash_path, '%s.Cumulus_barcodes.npy' % channel)), '-1')

		# Loop through barcodes and find hashing assignment
		aligned_HTO_assignment = np.empty(B.size, dtype='<U32')
		for B_i, curr_B in enumerate(B):
			cell_ind_in_HTO = np.where(curr_B == HTO_barcodes)[0]
			if cell_ind_in_HTO.size == 1:
				aligned_HTO_assignment[B_i] = HTO_assignment[cell_ind_in_HTO][0]
			elif cell_ind_in_HTO.size == 0:
				aligned_HTO_assignment[B_i] = 'unknown'
			else:
				bp()

		assert(aligned_HTO_assignment.size == B.size)
		assert(gene_counts_per_cell.size == B.size)
		assert(MT_perc.size == B.size)

		# Inds to delete
		delete_cell_inds = np.where((gene_counts_per_cell < MIN_GENES_PER_CELL) | (gene_counts_per_cell > MAX_GENES_PER_CELL))[0] # Filter based on genes / cell
		delete_cell_inds = np.union1d(delete_cell_inds, np.where(MT_perc > MAX_PERCENT_MITO)[0]) # Filter based on mitochondrial percent
		delete_cell_inds = np.union1d(delete_cell_inds, np.where(aligned_HTO_assignment != 'singlet')[0]) # Filter based on hashing

		# Calculate percent of lost cells
		lost_cells_perc = 100*(delete_cell_inds.size / B.size)

		# Delete bad cells
		B = np.delete(B, delete_cell_inds)
		curr_expr = np.delete(curr_expr, delete_cell_inds, axis=0)

		# Build save matrix gene-wise (gene names non-unique so sum)
		save_expr = np.zeros((B.size, filtered_genes.size))
		for gene_i, filtered_gene in enumerate(filtered_genes):
			G_ind = np.where(G == filtered_gene)[0]
			assert(G_ind.size > 0)
			save_expr[:, gene_i] = np.sum(curr_expr[:, G_ind], axis=1)

		clean_barcodes = np.array([b.split('-')[0] for b in B])

		assert(save_expr.shape[0] == clean_barcodes.size)
		assert(save_expr.shape[1] == filtered_genes.size)

		np.save(os.path.join(save_pn, '%s.expr_mat.npy' % channel), save_expr)
		np.save(os.path.join(save_pn, '%s.barcodes.npy' % channel), clean_barcodes.astype(np.str))

		plot_df_cells = plot_df_cells.append(pd.DataFrame({'Sample' : channel, 'Number of Cells' : [B.size]}), ignore_index=True)
		plot_df_perc = plot_df_perc.append(pd.DataFrame({'Sample' : channel, 'Percent Filtered' : [lost_cells_perc]}), ignore_index=True)

	bp()

	plot_df_cells = plot_df_cells.sort_values(by='Sample')
	fig_save_pn = os.path.join(fig_pn, 'filtered_cell_count.png')
	fig, ax = plt.subplots(figsize = [1*6.4, 1*4.8])
	sns.barplot(data=plot_df_cells, x='Sample', y='Number of Cells', ax=ax, color='b')
	plt.xticks(rotation=90), ax.set_xlabel('')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	plot_df_perc = plot_df_perc.sort_values(by='Sample')
	fig_save_pn = os.path.join(fig_pn, 'filtered_cell_percent.png')
	fig, ax = plt.subplots(figsize = [1*6.4, 1*4.8])
	sns.barplot(data=plot_df_perc, x='Sample', y='Percent Filtered', ax=ax, color='b')
	plt.xticks(rotation=90), ax.set_xlabel('')
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

if __name__ == "__main__":
    main()
