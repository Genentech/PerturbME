import numpy as np
import pandas as pd
import seaborn as sns
import os,glob
import scipy.sparse as sp_sparse
import scipy.stats as sp_stats
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
import pickle
from adjustText import adjust_text
from tqdm import tqdm

def plot_volcano(eff_size, p_vals, dot_labels, fig_save_pn, title='',
				 lfc_thresh=1, p_thresh=0.01, DPI=800, label_fontsize=2,
				 dot_size=1, FDR_correction=False, alpha_val=0.05,
				 FDR_method='fdr_bh', fig_size_mult=1, adjust_labels=True,
				 draw_arrows=False, arrow_style='-', arrow_color='red',
				 shade_program=False, program_inds=None, program_color='r',
				 program_name='', legend=True, arrow_width=0.5):

	fig, ax = plt.subplots(figsize=[fig_size_mult*6.4, fig_size_mult*4.8])

	if shade_program:
		non_prog_inds = np.delete(np.arange(dot_labels.size), program_inds)
		ax.scatter(eff_size[non_prog_inds], -np.log10(p_vals[non_prog_inds]), s=dot_size, c='b')
		ax.scatter(eff_size[program_inds], -np.log10(p_vals[program_inds]), s=dot_size, c=program_color, label=program_name)
	else:
		ax.scatter(eff_size, -np.log10(p_vals), s=dot_size)

	ax.set_xlabel('LFC'), ax.set_ylabel(u'-log\u2081\u2080(p)')

	if FDR_correction:
		label_inds = np.where(multipletests(p_vals, alpha=alpha_val, method=FDR_method)[0])[0]
	else:
		label_inds = np.where( (np.abs(eff_size) > lfc_thresh) & (p_vals < p_thresh) )[0]

	xlabels = eff_size[label_inds]
	ylabels = -np.log10(p_vals)[label_inds]
	label_text = dot_labels[label_inds]
	text_objs = [ax.text(xlabels[ind], ylabels[ind], label_text[ind], fontsize=label_fontsize, fontweight='bold') for ind in np.arange(label_inds.size)]

	if adjust_labels:
		if draw_arrows:
			adjust_text(text_objs, arrowprops=dict(arrowstyle=arrow_style, color=arrow_color, lw=arrow_width))
		else:
			adjust_text(text_objs)

	ax.set_title(title)
	if shade_program and legend:
		ax.legend()
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=DPI)
	plt.close()

def LFC_by_target(target_names, target_to_cell_dict, score_arr):
	print('Scoring')

	all_mean = np.mean(score_arr)

	LFC = np.zeros(target_names.size)
	for target_i, target in tqdm(enumerate(target_names)):

		if target != 'NO_SITE' and target != 'ONE_NON-GENE_SITE':

			curr_cells = target_to_cell_dict[target_i]

			if curr_cells.size > 5:
				curr_scores = score_arr[curr_cells]

				LFC[target_i] = np.log2( np.mean(curr_scores) / all_mean )

	return LFC

def make_target_to_cell_dict(target_cov_matrix):

	target_to_cell_dict = {}
	for target_i in tqdm(np.arange(target_cov_matrix.shape[1])):
		cells_for_this_target = np.where(target_cov_matrix[:,target_i] != 0)[0]

		target_to_cell_dict[target_i] = cells_for_this_target

	return target_to_cell_dict

def main():

	exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	figures_pn = os.path.join(exp_pn, 'figures/CROP')
	sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
	CROP_dir = os.path.join(exp_pn, 'CROP')
	design_dir = os.path.join(CROP_dir, 'design_mats/')
	load_pn = os.path.join(exp_pn, 'adata_arrs')
	lm_pn = os.path.join(exp_pn, 'CROP/linear_model')

	# adata_RNA_CITE = anndata.read(os.path.join(exp_pn, 'adata_RNA_CITE.h5ad'))

	target_names = np.load(os.path.join(design_dir, 'target_names.npy'))
	MOI_arr = np.load(os.path.join(load_pn, 'adata_parts/MOI_num.npy'))
	MOI1_inds = np.where(MOI_arr == 1)[0]
	cond_names = np.load(os.path.join(exp_pn, 'adata_arrs/adata_parts/cond_categories.npy'))[MOI1_inds]

	target_to_cell_dict_save_pn = os.path.join(load_pn, 'target_to_cell_dict.pkl')
	if os.path.exists(target_to_cell_dict_save_pn):
		target_to_cell_dict = pickle.load(open(target_to_cell_dict_save_pn, 'rb'))
	else:
		target_cov_matrix = np.array(sp_sparse.load_npz(os.path.join(lm_pn, 'MOI1_target_design_mat.npz')).todense())
		target_to_cell_dict = make_target_to_cell_dict(target_cov_matrix)
		pickle.dump(target_to_cell_dict, open(target_to_cell_dict_save_pn, 'wb'))

	num_bins = 50
	HLA_upregulation_score = np.load(os.path.join(load_pn, 'HLA_upregulation_score_%d_bins.npy' % num_bins))[MOI1_inds]
	HLA_downregulation_score = np.load(os.path.join(load_pn, 'HLA_downregulation_score_%d_bins.npy' % num_bins))[MOI1_inds]
	HLA_composite_score = np.load(os.path.join(load_pn, 'HLA_composite_score_%d_bins.npy' % num_bins))[MOI1_inds]

	upregulation_LFC = LFC_by_target(target_names, target_to_cell_dict, HLA_upregulation_score)
	downregulation_LFC = LFC_by_target(target_names, target_to_cell_dict, HLA_downregulation_score)
	composite_LFC = LFC_by_target(target_names, target_to_cell_dict, HLA_composite_score)

	bp()


if __name__ == "__main__":
    main()
