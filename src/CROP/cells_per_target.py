import numpy as np
import pandas as pd
import seaborn as sns
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt

SGRNAS_PER_TARGET = 3

# Assumes guides are labeled as "target_#" (e.g. JAK1_1, JAK1_2, JAK1_3) in sgRNA_names array
def collapse_sgRNAs_to_targets(sgRNA_names, sgRNA_design_mat, target_names):

	assert(sgRNA_design_mat.shape[1] == sgRNA_names.size)
	assert(target_names.size == np.unique(target_names).size)

	target_design_mat = np.zeros((sgRNA_design_mat.shape[0], target_names.size))
	for target_i, target in enumerate(target_names):

		sgRNA_inds = np.where(np.char.find(sgRNA_names, target+'_') == 0)[0]
		if sgRNA_inds.size != SGRNAS_PER_TARGET:
			assert( (target=='NO_SITE') or (target=='ONE_NON-GENE_SITE') )

		target_design_mat[:, target_i] = sgRNA_design_mat[:, sgRNA_inds].sum(1)

	return target_design_mat

def get_target_names(sgRNA_names):

	target_names = np.array([])
	for curr_sgRNA in sgRNA_names:
		if 'NO_SITE' in curr_sgRNA:
			target_names = np.append(target_names, 'NO_SITE')
		elif 'ONE_NON-GENE_SITE' in curr_sgRNA:
			target_names = np.append(target_names, 'ONE_NON-GENE_SITE')
		else:
			target_names = np.append(target_names, curr_sgRNA.split('_')[0])

	target_names = np.unique(target_names).astype(np.object).astype(np.str)

	return target_names

def main():

	exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	figures_pn = os.path.join(exp_pn, 'figures')
	sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
	CROP_dir = os.path.join(exp_pn, 'CROP')
	design_dir = os.path.join(CROP_dir, 'design_mats/')

	# Load sgRNA information
	gw_guides_pn = os.path.join(sequencing_info_pn, 'gw_lib_guides.csv')
	target_sgRNA_df = pd.read_csv(gw_guides_pn)
	target_sgRNA_names =  target_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
	control_sgRNA_df = pd.read_csv(os.path.join(sequencing_info_pn, 'human_TF_CTL_guides.txt'), sep='\t', header=None, names=['sgRNA_name', 'sgRNA_barcode'])
	control_sgRNA_names =  control_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)

	combined_sgRNA_names = np.append(target_sgRNA_names, control_sgRNA_names) # Targeting channels have both control and targeting sgRNAs

	all_design_mats = glob.glob(os.path.join(design_dir, '*.sgRNA_design_mat.*.npz'))
	assert(len(all_design_mats) == 20)

	target_names_pn = os.path.join(design_dir, 'target_names.npy')
	if os.path.exists(target_names_pn):
		target_names = np.load(target_names_pn)
	else:
		target_names = get_target_names(combined_sgRNA_names)
		np.save(target_names_pn, target_names)

	cells_per_target_high = np.zeros(target_names.size)
	cells_per_target_low = np.zeros(target_names.size)
	for curr_mat in all_design_mats:

		cond = curr_mat.split('/')[-1].split('.')[0]

		if 'High' in cond or 'Low' in cond:

			print('Processing', cond)

			target_design_mat_save_pn = curr_mat.replace('sgRNA_design_mat', 'target_design_mat')

			if not os.path.exists(target_design_mat_save_pn):

				sgRNA_design_mat = np.array(sp_sparse.load_npz(curr_mat).todense())
				assert(sgRNA_design_mat.shape[1] == combined_sgRNA_names.size)

				target_design_mat = collapse_sgRNAs_to_targets(combined_sgRNA_names, sgRNA_design_mat, target_names)

				# MOI should be unaffected when collapsing sgRNAs to targets
				assert(np.all(np.equal(sgRNA_design_mat.sum(1), target_design_mat.sum(1))))

				# Save target-wise design matrix
				sp_sparse.save_npz(target_design_mat_save_pn, sp_sparse.csr_matrix(target_design_mat))

			else:

				target_design_mat = np.array(sp_sparse.load_npz(target_design_mat_save_pn).todense())
				assert(target_design_mat.shape[1] == target_names.size)


			MOI = target_design_mat.sum(1)
			MOI1_cells = np.where(MOI == 1)[0]

			if 'High' in cond:
				cells_per_target_high += target_design_mat[MOI1_cells, :].sum(0)
			elif 'Low' in cond:
				cells_per_target_low += target_design_mat[MOI1_cells, :].sum(0)

	bp()
	np.save(os.path.join(CROP_dir, 'cells_per_target_high.npy'), cells_per_target_high)
	np.save(os.path.join(CROP_dir, 'cells_per_target_low.npy'), cells_per_target_low)



if __name__ == "__main__":
    main()
