import numpy as np
import pandas as pd
import glob,os
from scipy.spatial import distance
import scipy.stats as sp_stats
import fastcluster
from scipy.cluster.hierarchy import dendrogram
import anndata
import scanpy as sc
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import ElasticNet
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text
import pickle
import scipy.sparse as sp_sparse

def design_from_arrs(sgRNA_names, sgRNA_assignments):

	design_mat = np.zeros((sgRNA_assignments.size, sgRNA_names.size))

	for cell_i in range(sgRNA_assignments.size):
		sgRNA_ind = np.where(sgRNA_assignments[cell_i] == sgRNA_names)[0]
		assert(sgRNA_ind.size == 1)
		design_mat[cell_i, sgRNA_ind] = 1

	return design_mat

def plot_clustermap(plot_arr, row_inds, col_inds, row_names, col_names, save_pn,
					remove_sparse_rows=False, row_sparse_thresh=0.8, clim=1,
					cmap_name='seismic', DPI=800, xlabel='', ylabel='',
					label_ticks=False, fig_size_mult=1):

	fig, ax = plt.subplots(figsize=[fig_size_mult*6.4, fig_size_mult*4.8])
	filtered_plot_mat = plot_arr[row_inds, :][:, col_inds]

	if remove_sparse_rows:
		row_sparsity = np.array([ (np.where(filtered_plot_mat[curr_row, :] == 0)[0].size / col_inds.size) for curr_row in range(row_inds.size)])
		good_rows = np.where(row_sparsity < row_sparse_thresh)[0]
		filtered_plot_mat = filtered_plot_mat[good_rows, :]
		row_inds = row_inds[good_rows]

	ax_ret = ax.imshow(filtered_plot_mat, aspect=(filtered_plot_mat.shape[1]/filtered_plot_mat.shape[0]), cmap=cmap_name, clim=[-clim, clim])
	cbar = fig.colorbar(ax_ret, ax=ax, extend='both')
	cbar.minorticks_on()

	if label_ticks:
		ax.set_yticks(np.arange(row_inds.size))
		ax.set_yticklabels(row_names[row_inds], fontsize=1)

		ax.set_xticks(np.arange(col_inds.size))
		ax.set_xticklabels(col_names[col_inds], fontsize=1, rotation=90)
	else:
		ax.set_yticks([])
		ax.set_xticks([])

	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)

	ax.tick_params(width=0.1)
	ax.set_ylim([filtered_plot_mat.shape[0]-0.5, -0.5])

	fig.savefig(save_pn, bbox_inches='tight', dpi=DPI)
	plt.close()

def fit_and_cluster(X, y, l1_ratio=0.5, alpha=0.0005, max_iter=10000, z_score=False):

	lmfit = ElasticNet(precompute=True, l1_ratio=l1_ratio, alpha=alpha, max_iter=max_iter)

	if z_score:
		y = sp_stats.zscore(y, axis=0)

	lmfit.fit(X, y)

	# cgrid = sns.clustermap(lmfit.coef_)
	# row_inds = np.array(cgrid.dendrogram_row.reordered_ind)
	# col_inds = np.array(cgrid.dendrogram_col.reordered_ind)

	row_linkage = fastcluster.linkage(lmfit.coef_, method='complete')
	col_linkage = fastcluster.linkage(lmfit.coef_.T, method='complete')
	row_dend = dendrogram(row_linkage, no_plot=True)
	col_dend = dendrogram(col_linkage, no_plot=True)
	row_inds = np.array(row_dend['leaves'])
	col_inds = np.array(col_dend['leaves'])

	return lmfit.coef_, row_inds, col_inds, lmfit

# From https://github.com/asncd/MIMOSCA/blob/master/mimosca.py
def bayes_cov_col(Y,X,cols,lm,verbose=False):
	"""
	@Y    = Expression matrix, cells x x genes, expecting pandas dataframe
	@X    = Covariate matrix, cells x covariates, expecting pandas dataframe
	@cols = The subset of columns that the EM should be performed over, expecting list
	@lm   = linear model object
	"""

	#EM iterateit
	Yhat=pd.DataFrame(lm.predict(X))
	Yhat.index=Y.index
	Yhat.columns=Y.columns
	SSE_all=np.square(Y.subtract(Yhat))
	X_adjust=X.copy()


	df_SSE   = []
	df_logit = []

	for cov_i, curcov in enumerate(cols):

		if verbose:
			print('Column',cov_i+1, 'of', len(cols))

		curcells=X[X[curcov]>0].index

		if len(curcells)>2:

		    X_notcur=X.copy()
		    X_notcur[curcov]=[0]*len(X_notcur)

		    X_sub=X_notcur.loc[curcells]

		    Y_sub=Y.loc[curcells]

		    GENE_var=2.0*Y_sub.var(axis=0)
		    vargenes=GENE_var[GENE_var>0].index

		    Yhat_notcur=pd.DataFrame(lm.predict(X_sub))
		    Yhat_notcur.index=Y_sub.index
		    Yhat_notcur.columns=Y_sub.columns

		    SSE_notcur=np.square(Y_sub.subtract(Yhat_notcur))
		    SSE=SSE_all.loc[curcells].subtract(SSE_notcur)
		    SSE_sum=SSE.sum(axis=1)

		    SSE_transform=SSE.div(GENE_var+0.5)[vargenes].sum(axis=1)
		    logitify=np.divide(1.0,1.0+np.exp(SSE_transform))#sum))

		    df_SSE.append(SSE_sum)
		    df_logit.append(logitify)

		    X_adjust[curcov].loc[curcells]=logitify

	return X_adjust

def calc_corr_matrix(coefs):

	corr_mat = pd.DataFrame(coefs).corr().to_numpy()
	corr_mat[np.isnan(corr_mat)] = 0

	# cgrid = sns.clustermap(corr_mat)
	# row_inds = np.array(cgrid.dendrogram_row.reordered_ind)
	# col_inds = np.array(cgrid.dendrogram_col.reordered_ind)

	row_linkage = fastcluster.linkage(corr_mat, method='complete')
	col_linkage = fastcluster.linkage(corr_mat.T, method='complete')
	row_dend = dendrogram(row_linkage, no_plot=True)
	col_dend = dendrogram(col_linkage, no_plot=True)
	row_inds = np.array(row_dend['leaves'])
	col_inds = np.array(col_dend['leaves'])

	assert(np.all(row_inds==col_inds)) # Symmetric matrix

	return corr_mat, row_inds, col_inds

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

def main():

	plot = True
	save = True
	cells_per_target_thresh = 13

	full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	fig_pn = os.path.join(full_exp_pn, 'figures/CROP/linear_model')
	adata_arr_pn = os.path.join(full_exp_pn, 'adata_arrs')
	sequencing_info_pn = os.path.join(full_exp_pn, 'sequencing_info')
	design_mat_pn = os.path.join(full_exp_pn, 'CROP/design_mats')
	target_names = np.load(os.path.join(design_mat_pn, 'target_names.npy'))
	save_pn = os.path.join(full_exp_pn, 'CROP/linear_model')

	# Load sgRNA information
	gw_guides_pn = os.path.join(sequencing_info_pn, 'gw_lib_guides.csv')
	target_sgRNA_df = pd.read_csv(gw_guides_pn)
	target_sgRNA_names =  target_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
	control_sgRNA_df = pd.read_csv(os.path.join(sequencing_info_pn, 'human_TF_CTL_guides.txt'), sep='\t', header=None, names=['sgRNA_name', 'sgRNA_barcode'])
	control_sgRNA_names =  control_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
	combined_sgRNA_names = np.append(target_sgRNA_names, control_sgRNA_names) # Targeting channels have both control and targeting sgRNAs

	print('Loading AnnData')
	adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))
	feature_names = adata_RNA_CITE.var.index.to_numpy().astype(np.str)
	UMI_counts_per_cell = np.load(os.path.join(adata_arr_pn, 'UMI_count.npy'))
	cell_cycle = np.load(os.path.join(adata_arr_pn, 'cell_cycle_Tirosh.npy'))
	assert(UMI_counts_per_cell.size == adata_RNA_CITE.shape[0])
	assert(cell_cycle.size == adata_RNA_CITE.shape[0])

	# For now only run model with MOI = 1  cells
	MOI1_inds = np.where(adata_RNA_CITE.obs['MOI'] == 1)[0]

	# Covariates - cell cycle and UMI
	print('Building covariate matrix with %s cells / target threshold' % cells_per_target_thresh)
	cell_cycle_arr = cell_cycle[MOI1_inds]
	UMI_arr = UMI_counts_per_cell[MOI1_inds].astype(np.int)

	# Covariate - 10X run (there is clear bias in UMAP)
	sample_arr = adata_RNA_CITE.obs['Sample'].to_numpy().astype(np.str)[MOI1_inds]
	run1_inds = np.where(np.char.find(sample_arr, '1') != -1)[0]
	run2_inds = np.where(np.char.find(sample_arr, '2') != -1)[0]
	run3_inds = np.where(np.char.find(sample_arr, '3') != -1)[0]
	tenx_arr = np.full(sample_arr.size, '0')
	tenx_arr[run2_inds] = '2' # Assign 2 first so '12' condition is assigned to 10X run #1 (based on UMAP)
	tenx_arr[run1_inds] = '1'
	tenx_arr[run3_inds] = '3'
	assert(np.where(tenx_arr == '0')[0].size==0)

	# Put cell cycle and 10X covariates into matrix form
	cell_cycle_matrix = np.array([cell_cycle_arr == curr_cycle for curr_cycle in np.unique(cell_cycle)]).T.astype(np.int)
	assert(np.all(cell_cycle_matrix.sum(1) == 1)) # Ensure every cell assigned to one phase of cell cycle

	tenx_matrix = np.array([tenx_arr == curr_tenx for curr_tenx in np.unique(tenx_arr)]).T.astype(np.int)
	assert(np.all(tenx_matrix.sum(1) == 1)) # Ensure every cell assigned to one 10X channel

	target_cov_matrix_save_pn = os.path.join(save_pn, 'MOI1_target_design_mat.npz')
	if os.path.exists(target_cov_matrix_save_pn):
		target_cov_matrix = np.array(sp_sparse.load_npz(target_cov_matrix_save_pn).todense())
	else:
		curr_sgRNAs = adata_RNA_CITE.obs['sgRNA'].to_numpy().astype(np.str)[MOI1_inds]
		sgRNA_cov_matrix = design_from_arrs(combined_sgRNA_names, curr_sgRNAs)
		target_cov_matrix = collapse_sgRNAs_to_targets(combined_sgRNA_names, sgRNA_cov_matrix, target_names)
		sp_sparse.save_npz(target_cov_matrix_save_pn, sp_sparse.csr_matrix(target_cov_matrix))

	assert(target_cov_matrix.shape[1]==target_names.size)
	assert(target_cov_matrix.shape[0]==cell_cycle_matrix.shape[0])
	assert(target_cov_matrix.shape[0]==tenx_matrix.shape[0])
	assert(target_cov_matrix.shape[0]==UMI_arr.size)

	# Only keep targets / cells exceeding minimum cells / target threshold (only count cells from one channel)
	high_inds = np.where(np.char.startswith(sample_arr, 'High_'))[0]
	low_inds = np.where(np.char.startswith(sample_arr, 'Low_'))[0]

	cells_per_target_high = target_cov_matrix[high_inds, :].sum(0)
	cells_per_target_low = target_cov_matrix[low_inds, :].sum(0)
	assert(cells_per_target_high.size == target_names.size)
	assert(cells_per_target_low.size == target_names.size)

	target_inds_to_keep = np.where( (cells_per_target_high > cells_per_target_thresh) | (cells_per_target_low > cells_per_target_thresh) )[0]
	target_names_to_keep = target_names[target_inds_to_keep]

	# targets_we_know_we_want = np.array(['NO_SITE', 'ONE_NON-GENE_SITE', 'JAK1', 'JAK2', 'STAT1', 'IFNGR1', 'IFNGR2'])
	# assert(np.all(np.isin(targets_we_know_we_want, target_names_to_keep)))

	target_cov_matrix_filtered = target_cov_matrix[:, target_inds_to_keep] # Only keep relevant targets
	cells_per_target_filtered = target_cov_matrix_filtered.sum(0)

	# Delete cells that no longer have a sgRNA
	targets_per_cell = target_cov_matrix_filtered.sum(1)
	assert(targets_per_cell.size == target_cov_matrix.shape[0])
	assert(np.where( (targets_per_cell!=0) & (targets_per_cell!=1) )[0].size==0) # Every cell should have 0 or 1 sgRNA
	cells_with_remaining_sgRNA = np.where(targets_per_cell==1)[0]

	# Get info for cells / target covariate
	MOI1_cov_mat = target_cov_matrix_filtered[cells_with_remaining_sgRNA, :]
	cells_per_target_target_ind = np.array([np.where(MOI1_cov_mat[curr_cell,:]==1)[0] for curr_cell in np.arange(cells_with_remaining_sgRNA.size)]).flatten()
	cells_per_target_cov_arr = cells_per_target_filtered[cells_per_target_target_ind]

	# Build covariate matrix
	assert(target_cov_matrix_filtered.shape[0]==cell_cycle_matrix.shape[0])
	assert(target_cov_matrix_filtered.shape[0]==tenx_matrix.shape[0])
	assert(target_cov_matrix_filtered.shape[0]==UMI_arr.size)
	assert(cells_with_remaining_sgRNA.size==cells_per_target_cov_arr.size)
	full_cov_mat = np.column_stack((target_cov_matrix_filtered[cells_with_remaining_sgRNA, :], cell_cycle_matrix[cells_with_remaining_sgRNA, :]))
	full_cov_mat = np.column_stack((full_cov_mat, tenx_matrix[cells_with_remaining_sgRNA, :]))
	full_cov_mat = np.column_stack((full_cov_mat, UMI_arr[cells_with_remaining_sgRNA]))
	full_cov_mat = np.column_stack((full_cov_mat, cells_per_target_cov_arr))
	all_cov_names = np.append(target_names_to_keep, np.append(np.unique(cell_cycle), np.append(np.char.add('TENX_', np.unique(tenx_arr)), 'UMI')))
	all_cov_names = np.append(all_cov_names, 'cells_per_target')
	assert(full_cov_mat.shape[1]==all_cov_names.size)
	assert(full_cov_mat.shape[0]==cells_with_remaining_sgRNA.size)

	# Build expression matrix
	print('Building expression matrix')
	raw_expr_mat = adata_RNA_CITE.X[MOI1_inds, :]
	assert(raw_expr_mat.shape[0]==target_cov_matrix.shape[0])
	assert(raw_expr_mat.shape[1]==feature_names.size)

	# Do feature selection on MOI=1 cells
	DE_features = np.load(os.path.join(save_pn, 'lm_features_from_DE.npy')) # Features from differential expression analysis
	assert(np.all(np.isin(DE_features, feature_names)))
	DE_inds = np.where(np.isin(feature_names, DE_features))[0]

	adata_temp = anndata.AnnData(X=raw_expr_mat)
	sc.pp.highly_variable_genes(adata_temp)
	assert(adata_temp.var['highly_variable'].size == feature_names.size)
	highly_variable_inds = np.where(adata_temp.var['highly_variable'])[0]
	del adata_temp

	features_that_are_targeted = np.where(np.isin(feature_names, target_names_to_keep))[0]
	CITE_inds = np.where(np.char.find(feature_names, 'CITE-') == 0)[0]

	feature_inds_to_keep = np.union1d(highly_variable_inds, np.union1d(features_that_are_targeted, CITE_inds))
	feature_inds_to_keep = np.union1d(feature_inds_to_keep, DE_inds)
	feature_names_to_keep = feature_names[feature_inds_to_keep]

	assert(np.all(np.isin(DE_features, feature_names_to_keep)))

	print('Selected %s features' % str(feature_inds_to_keep.size))

	full_expr_mat = raw_expr_mat[cells_with_remaining_sgRNA, :][:, feature_inds_to_keep]
	assert(full_expr_mat.shape[0]==full_cov_mat.shape[0])
	assert(full_expr_mat.shape[1]==feature_names_to_keep.size)

	if save:
		sp_sparse.save_npz(os.path.join(save_pn, '%s_cells_per_target/full_cov_mat.npz' % cells_per_target_thresh), sp_sparse.csr_matrix(full_cov_mat))
		sp_sparse.save_npz(os.path.join(save_pn, '%s_cells_per_target/full_expr_mat.npz' % cells_per_target_thresh), sp_sparse.csr_matrix(full_expr_mat))
		np.save(os.path.join(save_pn, '%s_cells_per_target/cov_names.npy' % cells_per_target_thresh), all_cov_names)
		np.save(os.path.join(save_pn, '%s_cells_per_target/feature_names.npy' % cells_per_target_thresh), feature_names_to_keep)

	# Run elastic net
	print('Elastic Net #1')
	EN_B, EN_B_row_inds, EN_B_col_inds, EN_lm = fit_and_cluster(full_cov_mat, full_expr_mat)
	if save:
		np.save(os.path.join(save_pn, '%s_cells_per_target/EN_B.npy' % cells_per_target_thresh), EN_B)
		np.save(os.path.join(save_pn, '%s_cells_per_target/EN_B_row_inds.npy' % cells_per_target_thresh), EN_B_row_inds)
		np.save(os.path.join(save_pn, '%s_cells_per_target/EN_B_col_inds.npy' % cells_per_target_thresh), EN_B_col_inds)
		pickle.dump(EN_lm, open(os.path.join(save_pn, '%s_cells_per_target/EN_lm.pkl' % cells_per_target_thresh), 'wb'))

	if plot:
		plot_clustermap(EN_B, EN_B_row_inds, EN_B_col_inds, feature_names_to_keep, all_cov_names, os.path.join(fig_pn, 'EN_B_%s_cells_per_target.png' % cells_per_target_thresh),
						xlabel='Covariates', ylabel='Features', clim=0.3, label_ticks=True, fig_size_mult=3, remove_sparse_rows=True, row_sparse_thresh=0.9)

	# EM-like update
	print('EM-like update')
	targets_for_EM = [target for target in target_names_to_keep if ('NO_SITE' not in target and 'ONE_NON-GENE_SITE' not in target)]

	Y_df = pd.DataFrame(full_expr_mat, columns=feature_names_to_keep)
	X_df = pd.DataFrame(full_cov_mat, columns=all_cov_names)
	full_cov_mat_EM = bayes_cov_col(Y_df,X_df,targets_for_EM,EN_lm,verbose=False)

	if save:
		sp_sparse.save_npz(os.path.join(save_pn, '%s_cells_per_target/full_cov_mat_EM.npz' % cells_per_target_thresh), sp_sparse.csr_matrix(full_cov_mat_EM))

	# Second elastic net
	print('Elastic Net #2')
	EN_B_EM, EN_B_row_inds_EM, EN_B_col_inds_EM, EN_lm_EM = fit_and_cluster(full_cov_mat_EM.to_numpy(), full_expr_mat)
	if save:
		np.save(os.path.join(save_pn, '%s_cells_per_target/EN_B_EM.npy' % cells_per_target_thresh), EN_B_EM)
		np.save(os.path.join(save_pn, '%s_cells_per_target/EN_B_row_inds_EM.npy' % cells_per_target_thresh), EN_B_row_inds_EM)
		np.save(os.path.join(save_pn, '%s_cells_per_target/EN_B_col_inds_EM.npy' % cells_per_target_thresh), EN_B_col_inds_EM)
		pickle.dump(EN_lm_EM, open(os.path.join(save_pn, '%s_cells_per_target/EN_lm_EM.pkl' % cells_per_target_thresh), 'wb'))

	if plot:
		plot_clustermap(EN_B_EM, EN_B_row_inds_EM, EN_B_col_inds_EM, feature_names_to_keep, all_cov_names, os.path.join(fig_pn, 'EN_B_EM_%s_cells_per_target.png' % cells_per_target_thresh),
						xlabel='Covariates', ylabel='Features', clim=0.3, label_ticks=True, fig_size_mult=3, remove_sparse_rows=True, row_sparse_thresh=0.9)

	bp()

if __name__ == "__main__":
    main()
