import numpy as np
import pandas as pd
import seaborn as sns
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from adjustText import adjust_text

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
	gDNA_pn = os.path.join(exp_pn, 'guide_seq')
	mageck_out_pn = os.path.join(gDNA_pn, 'mageck_out')

	# adata_RNA_CITE = anndata.read(os.path.join(exp_pn, 'adata_RNA_CITE.h5ad'))

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

	cells_per_target_high_MOI1 = np.zeros(target_names.size)
	cells_per_target_high_MOI2 = np.zeros(target_names.size)
	cells_per_target_high_MOI3 = np.zeros(target_names.size)
	cells_per_target_high_MOI4 = np.zeros(target_names.size)
	cells_per_target_high_MOI5 = np.zeros(target_names.size)
	cells_per_target_low_MOI1 = np.zeros(target_names.size)
	cells_per_target_low_MOI2 = np.zeros(target_names.size)
	cells_per_target_low_MOI3 = np.zeros(target_names.size)
	cells_per_target_low_MOI4 = np.zeros(target_names.size)
	cells_per_target_low_MOI5 = np.zeros(target_names.size)
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
			MOI2_cells = np.where(MOI == 2)[0]
			MOI3_cells = np.where(MOI == 3)[0]
			MOI4_cells = np.where(MOI == 4)[0]
			MOI5_cells = np.where(MOI == 5)[0]

			if 'High' in cond:
				cells_per_target_high_MOI1 += target_design_mat[MOI1_cells, :].sum(0)
				cells_per_target_high_MOI2 += target_design_mat[MOI2_cells, :].sum(0)
				cells_per_target_high_MOI3 += target_design_mat[MOI3_cells, :].sum(0)
				cells_per_target_high_MOI4 += target_design_mat[MOI4_cells, :].sum(0)
				cells_per_target_high_MOI5 += target_design_mat[MOI5_cells, :].sum(0)
			elif 'Low' in cond:
				cells_per_target_low_MOI1 += target_design_mat[MOI1_cells, :].sum(0)
				cells_per_target_low_MOI2 += target_design_mat[MOI2_cells, :].sum(0)
				cells_per_target_low_MOI3 += target_design_mat[MOI3_cells, :].sum(0)
				cells_per_target_low_MOI4 += target_design_mat[MOI4_cells, :].sum(0)
				cells_per_target_low_MOI5 += target_design_mat[MOI5_cells, :].sum(0)

	bp()

	reads_per_gene_high = np.load(os.path.join(mageck_out_pn, 'reads_per_gene_high.npy'))
	reads_per_gene_low = np.load(os.path.join(mageck_out_pn, 'reads_per_gene_low.npy'))
	all_genes = np.load(os.path.join(mageck_out_pn, 'all_genes.npy'))

	control_inds = np.where( (all_genes=='NO_SITE') | (all_genes=='ONE_NON-GENE_SITE'))[0]

	target_to_gene_map = np.concatenate([np.where(curr_gene==target_names)[0] for curr_gene in all_genes])
	assert(np.all(np.char.equal(all_genes, target_names[target_to_gene_map])))

	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_low, cells_per_target_low_MOI1[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA Low MOI = 1')
	ax.set_xlim([0, 3e5]), ax.set_ylim([-5, 200])
	label_inds = np.where(cells_per_target_low_MOI1[target_to_gene_map] > 40)[0]
	label_inds = np.union1d(label_inds, np.where( (reads_per_gene_low > 75e3) & (cells_per_target_low_MOI1[target_to_gene_map] > 40) )[0])
	label_inds = np.delete(label_inds, np.where(np.isin(label_inds, control_inds))[0])
	text_objs = [ax.text(reads_per_gene_low[ind], cells_per_target_low_MOI1[target_to_gene_map][ind], all_genes[ind], fontsize=6, fontweight='bold') for ind in label_inds]
	adjust_text(text_objs, arrowprops=dict(arrowstyle='-', color='red'))
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_low_MOI1.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_low, cells_per_target_low_MOI2[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA Low MOI = 2')
	ax.set_xlim([0, 3e5]), ax.set_ylim([-5, 200])
	label_inds = np.where(reads_per_gene_low > 75e3)[0]
	label_inds = np.delete(label_inds, np.where(np.isin(label_inds, control_inds))[0])
	text_objs = [ax.text(reads_per_gene_low[ind], cells_per_target_low_MOI2[target_to_gene_map][ind], all_genes[ind], fontsize=6, fontweight='bold') for ind in label_inds]
	adjust_text(text_objs, arrowprops=dict(arrowstyle='-', color='red'))
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_low_MOI2.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_low, cells_per_target_low_MOI3[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA Low MOI = 3')
	ax.set_xlim([0, 3e5]), ax.set_ylim([-5, 200])
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_low_MOI3.png'), bbox_inches='tight', dpi=800)
	plt.close()
	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_low, cells_per_target_low_MOI4[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA Low MOI = 4')
	ax.set_xlim([0, 3e5]), ax.set_ylim([-5, 200])
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_low_MOI4.png'), bbox_inches='tight', dpi=800)
	plt.close()
	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_low, cells_per_target_low_MOI5[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA Low MOI = 5')
	ax.set_xlim([0, 3e5]), ax.set_ylim([-5, 200])
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_low_MOI5.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_high, cells_per_target_high_MOI1[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA High MOI = 1')
	ax.set_xlim([0, 150e3]), ax.set_ylim([-5, 80])
	label_inds = np.where(cells_per_target_high_MOI1[target_to_gene_map] > 20)[0]
	label_inds = np.union1d(label_inds, np.where( (reads_per_gene_high > 40e3) & (cells_per_target_high_MOI1[target_to_gene_map] > 20) )[0])
	label_inds = np.delete(label_inds, np.where(np.isin(label_inds, control_inds))[0])
	text_objs = [ax.text(reads_per_gene_high[ind], cells_per_target_high_MOI1[target_to_gene_map][ind], all_genes[ind], fontsize=6, fontweight='bold') for ind in label_inds]
	adjust_text(text_objs, arrowprops=dict(arrowstyle='-', color='red'))
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_high_MOI1.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_high, cells_per_target_high_MOI2[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA High MOI = 2')
	ax.set_xlim([0, 150e3]), ax.set_ylim([-5, 80])
	label_inds = np.where(reads_per_gene_high > 50e3)[0]
	label_inds = np.delete(label_inds, np.where(np.isin(label_inds, control_inds))[0])
	text_objs = [ax.text(reads_per_gene_high[ind], cells_per_target_high_MOI2[target_to_gene_map][ind], all_genes[ind], fontsize=6, fontweight='bold') for ind in label_inds]
	adjust_text(text_objs, arrowprops=dict(arrowstyle='-', color='red'))
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_high_MOI2.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_high, cells_per_target_high_MOI3[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA High MOI = 3')
	ax.set_xlim([0, 150e3]), ax.set_ylim([-5, 80])
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_high_MOI3.png'), bbox_inches='tight', dpi=800)
	plt.close()
	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_high, cells_per_target_high_MOI4[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA High MOI = 4')
	ax.set_xlim([0, 150e3]), ax.set_ylim([-5, 80])
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_high_MOI4.png'), bbox_inches='tight', dpi=800)
	plt.close()
	fig, ax = plt.subplots()
	ax.scatter(reads_per_gene_high, cells_per_target_high_MOI5[target_to_gene_map], s=1)
	ax.set_xlabel('Reads / Target'), ax.set_ylabel('Cells / Target')
	ax.set_title('HLA High MOI = 5')
	ax.set_xlim([0, 150e3]), ax.set_ylim([-5, 80])
	fig.savefig(os.path.join(figures_pn, 'reads_vs_cells_per_target_high_MOI5.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(cells_per_target_low_MOI1, cells_per_target_low_MOI2, s=1)
	ax.set_xlabel('MOI = 1'), ax.set_ylabel('MOI = 2')
	ax.set_title('Low Cells / Target')
	ax.set_xlim([-5, 200]), ax.set_ylim([-5, 200])
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_low_MOI1_vs_MOI2.png'), bbox_inches='tight', dpi=800)
	ax.set_xlim([-5, 50]), ax.set_ylim([-5, 50])
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_low_MOI1_vs_MOI2_zoom.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(cells_per_target_high_MOI1, cells_per_target_high_MOI2, s=1)
	ax.set_xlabel('MOI = 1'), ax.set_ylabel('MOI = 2')
	ax.set_title('High Cells / Target')
	ax.set_xlim([-5, 80]), ax.set_ylim([-5, 80])
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_high_MOI1_vs_MOI2.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(cells_per_target_low_MOI1, cells_per_target_low_MOI2, s=1, label='X = 2')
	ax.scatter(cells_per_target_low_MOI1, cells_per_target_low_MOI3, s=1, label='X = 3')
	ax.scatter(cells_per_target_low_MOI1, cells_per_target_low_MOI4, s=1, label='X = 4')
	ax.scatter(cells_per_target_low_MOI1, cells_per_target_low_MOI5, s=1, label='X = 5')
	ax.set_xlabel('MOI = 1'), ax.set_ylabel('MOI = X')
	ax.set_title('Low Cells / Target')
	ax.set_xlim([-5, 200]), ax.set_ylim([-5, 200])
	ax.legend()
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_low_MOI1_vs_MOIX.png'), bbox_inches='tight', dpi=800)
	ax.set_xlim([-5, 50]), ax.set_ylim([-5, 50])
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_low_MOI1_vs_MOIX_zoom.png'), bbox_inches='tight', dpi=800)
	plt.close()

	fig, ax = plt.subplots()
	ax.scatter(cells_per_target_high_MOI1, cells_per_target_high_MOI2, s=1, label='X = 2')
	ax.scatter(cells_per_target_high_MOI1, cells_per_target_high_MOI3, s=1, label='X = 3')
	ax.scatter(cells_per_target_high_MOI1, cells_per_target_high_MOI4, s=1, label='X = 4')
	ax.scatter(cells_per_target_high_MOI1, cells_per_target_high_MOI5, s=1, label='X = 5')
	ax.set_xlabel('MOI = 1'), ax.set_ylabel('MOI = X')
	ax.set_title('High Cells / Target')
	ax.set_xlim([-5, 80]), ax.set_ylim([-5, 80])
	ax.legend()
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_high_MOI1_vs_MOIX.png'), bbox_inches='tight', dpi=800)
	plt.close()


if __name__ == "__main__":
    main()
