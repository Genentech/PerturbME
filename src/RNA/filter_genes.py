import numpy as np
import pandas as pd
import scanpy as sc
import os,glob
import h5py
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import seaborn as sns
import matplotlib.pyplot as plt

MIN_CELLS_PER_GENE = 100

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

	full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	fig_pn = os.path.join(full_exp_pn, 'figures')
	load_pn = os.path.join(full_exp_pn, 'filtered_RNA')

	all_genes = np.load(os.path.join(load_pn, 'all_genes.npy'))
	cell_counts_per_gene = np.zeros(all_genes.size)

	mat_pns = glob.glob(os.path.join(load_pn, '*.npy'))
	mat_pns = [mat_pn for mat_pn in mat_pns if 'expr_mat' in mat_pn]
	for mat_i, mat_pn in enumerate(mat_pns):

		channel = mat_pn.split('/')[-1].split('.')[0]
		print('Processing %s' % channel)

		curr_expr_mat = np.load(mat_pn)

		assert(all_genes.size == curr_expr_mat.shape[1])

		cell_counts_per_gene += np.count_nonzero(curr_expr_mat, axis=0)

	fig_save_pn = os.path.join(fig_pn, 'cells_per_gene.png')
	fig, ax = plt.subplots()
	ax.hist(cell_counts_per_gene, bins=200)
	ax.set_xlabel('Cells / Gene'), ax.set_ylabel('Count')
	ax.set_ylim([0, 2e3]), ax.set_xlim([0, 50e3])
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	# Gene indices to keep
	good_gene_inds = np.where(cell_counts_per_gene >= MIN_CELLS_PER_GENE)[0]

	# Save filtered gene list
	print('Saving filtered gene list with %d genes' % good_gene_inds.size)
	np.save(os.path.join(load_pn, 'filtered_gene_inds_%d_cells.npy' % MIN_CELLS_PER_GENE), good_gene_inds)
	np.save(os.path.join(load_pn, 'filtered_gene_names_%d_cells.npy' % MIN_CELLS_PER_GENE), all_genes[good_gene_inds])

if __name__ == "__main__":
    main()
