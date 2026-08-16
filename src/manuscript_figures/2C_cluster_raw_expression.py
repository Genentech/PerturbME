import numpy as np
import pandas as pd
import glob,os
from scipy.sparse import csr_matrix
from scipy.spatial import distance
import scipy.stats as sp_stats
import fastcluster
from scipy.cluster.hierarchy import dendrogram
import anndata
import scanpy as sc
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from sklearn.linear_model import ElasticNet
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.stats.multitest import multipletests
from sklearn.cluster import KMeans
from adjustText import adjust_text
import pickle

def plot_clustermap(plot_arr, row_inds, col_inds, row_names, col_names, save_pn,
					remove_sparse_rows=False, row_sparse_thresh=0.8, clim=1,
					cmap_name='seismic', DPI=800, xlabel='', ylabel='',
					label_ticks=False, fig_size_mult=1, title='', symmetric=True):

    fig, ax = plt.subplots(figsize=[fig_size_mult*6.4, fig_size_mult*4.8])
    filtered_plot_mat = plot_arr[row_inds, :][:, col_inds]

    if remove_sparse_rows:
    	row_sparsity = np.array([ (np.where(filtered_plot_mat[curr_row, :] == 0)[0].size / col_inds.size) for curr_row in range(row_inds.size)])
    	good_rows = np.where(row_sparsity < row_sparse_thresh)[0]
    	filtered_plot_mat = filtered_plot_mat[good_rows, :]
    	row_inds = row_inds[good_rows]

    if symmetric:
        ax_ret = ax.imshow(filtered_plot_mat, aspect=(filtered_plot_mat.shape[1]/filtered_plot_mat.shape[0]), cmap=cmap_name, clim=[-clim, clim])
    else:
        ax_ret = ax.imshow(filtered_plot_mat, aspect=(filtered_plot_mat.shape[1]/filtered_plot_mat.shape[0]), cmap=cmap_name, clim=[0, clim])

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
    ax.set_title(title)

    ax.tick_params(width=0.1)
    ax.set_ylim([filtered_plot_mat.shape[0]-0.5, -0.5])

    fig.savefig(save_pn, bbox_inches='tight', dpi=DPI)
    plt.close()

def calc_corr_matrix(mat):

	corr_mat = pd.DataFrame(mat).corr().to_numpy()
	corr_mat[np.isnan(corr_mat)] = 0

	return corr_mat

def cluster_by_kmeans(mat, num_clusters):

	km = KMeans(n_clusters=num_clusters, random_state=3)
	km.fit(mat)

	ret_order = np.array([], dtype=np.int)
	for clust_i in np.arange(num_clusters):
		curr_clust_inds = np.where(clust_i == km.labels_)[0]
		ret_order = np.append(ret_order, curr_clust_inds.astype(np.int))

	assert(ret_order.size == km.labels_.size)

	return ret_order, km.labels_[ret_order]

def cluster_mat_fastcluster(mat, type='complete'):

	row_linkage = fastcluster.linkage(mat, method=type)
	col_linkage = fastcluster.linkage(mat.T, method=type)
	row_dend = dendrogram(row_linkage, no_plot=True)
	col_dend = dendrogram(col_linkage, no_plot=True)
	row_inds = np.array(row_dend['leaves'])
	col_inds = np.array(col_dend['leaves'])

	return row_inds, col_inds


def calculate_inertia(mat, max_clusters=20):

	clusters_to_test = np.arange(max_clusters)+1

	all_distances = np.zeros(clusters_to_test.size)

	for cluster_i, num_clusters in enumerate(clusters_to_test):
		km = KMeans(n_clusters=num_clusters)
		km.fit(mat)
		all_distances[cluster_i] = km.inertia_

	return clusters_to_test, all_distances

def remove_sparse_rows_and_columns(mat, row_names, col_names, approx_zero = 0.02, row_sparse_thresh = 0.7, col_sparse_thresh = 0.7):

	row_sparsity = np.array([ (np.where( np.abs(mat[curr_row, :]) <= approx_zero)[0].size / col_names.size) for curr_row in range(row_names.size)])
	good_rows = np.where(row_sparsity < row_sparse_thresh)[0]
	filtered_rows = row_names[good_rows]

	col_sparsity = np.array([ (np.where( np.abs(mat[:, curr_col]) <= approx_zero)[0].size / row_names.size) for curr_col in range(col_names.size)])
	good_cols = np.where(col_sparsity < col_sparse_thresh)[0]
	filtered_cols = col_names[good_cols]

	filtered_mat = mat[good_rows, :][:, good_cols]

	return filtered_mat, filtered_rows, filtered_cols

def main():

    full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
    fig_pn = os.path.join(full_exp_pn, 'figures/manuscript_figures/Fig2')
    CROP_pn = os.path.join(full_exp_pn, 'CROP')
    DE_pn = os.path.join(full_exp_pn, 'DE/DE_csvs')

    adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))

    high_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'High')[0]
    low_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Low')[0]
    control_inds = np.where(adata_RNA_CITE.obs['Condition'] == 'Control')[0]
    cond_arr = adata_RNA_CITE.obs['Condition'].to_numpy().astype(str)
    features = adata_RNA_CITE.var.index.to_numpy().astype(str)

	# Isolate to cells with enriched sgRNA
    num_cells_enriched = 17
    cells_per_target = pd.read_csv(os.path.join(CROP_pn, 'cells_per_target.csv'), index_col=0)
    sgRNA_names = cells_per_target.index.to_numpy().astype(str)
    highE_guides = sgRNA_names[np.where(cells_per_target['HLA High'].to_numpy()>num_cells_enriched)[0]]
    lowE_guides = sgRNA_names[np.where(cells_per_target['HLA Low'].to_numpy()>num_cells_enriched)[0]]

    cell_guides = np.array([curr_sgRNA.split('_')[0] for curr_sgRNA in adata_RNA_CITE.obs['sgRNA'].to_numpy().astype(str)])
    highE_cells = np.intersect1d(np.where(np.isin(cell_guides, highE_guides))[0], np.where(adata_RNA_CITE.obs['Condition'].to_numpy().astype(str)=='High'))
    lowE_cells = np.intersect1d(np.where(np.isin(cell_guides, lowE_guides))[0], np.where(adata_RNA_CITE.obs['Condition'].to_numpy().astype(str)=='Low'))

    expr_mat = adata_RNA_CITE.X
    # cell_inds = np.random.choice(expr_mat.shape[0], size=10000, replace=False)
    cell_inds = np.union1d(highE_cells, lowE_cells)
    enriched_cell_guides = cell_guides[cell_inds]

    # HLA program
    FDR_val = 0.25
    LFC_val = 0.2
    highE_lowE_df = pd.read_csv(os.path.join(DE_pn, 'RNA_DE_results_MAST_highE_lowE.csv'))
    highE_lowE_features = highE_lowE_df['Unnamed: 0'].loc[ (highE_lowE_df['p_val_adj'] < FDR_val) & (np.abs(highE_lowE_df['avg_logFC']) > LFC_val)].to_numpy().astype(np.str)
    HLA_gene_inds = np.where(np.isin(features, highE_lowE_features))[0]

    downsample_expr = expr_mat[cell_inds, :][:, HLA_gene_inds]
    downsample_expr = sp_stats.zscore(downsample_expr, axis=0)
    downsample_conds = cond_arr[cell_inds]
    downsample_features = features[HLA_gene_inds]

    cell_clusters = 8
    features_clusters = 6

    cells_order, cells_cluster_labels = cluster_by_kmeans(downsample_expr, cell_clusters)

    unique_clusters = np.unique(cells_cluster_labels)
    low_frequency = np.zeros(unique_clusters.size)
    for cluster_i, unique_cluster in enumerate(unique_clusters):
        curr_cell_inds = cells_order[np.where(cells_cluster_labels==unique_cluster)[0]]

        curr_high_count = np.where(downsample_conds[curr_cell_inds] == 'High')[0].size
        curr_low_count = np.where(downsample_conds[curr_cell_inds] == 'Low')[0].size
        assert((curr_high_count+curr_low_count)==curr_cell_inds.size)

        low_frequency[cluster_i] = curr_low_count / curr_cell_inds.size

    sorted_clusters_by_low = unique_clusters[np.argsort(low_frequency)]

    refined_cells_order = np.array([], dtype=int)
    cells_cluster_labels_refined = np.array([], dtype=int)
    for unique_cluster in sorted_clusters_by_low:
        curr_cell_inds = cells_order[np.where(cells_cluster_labels==unique_cluster)[0]]
        curr_sorted_inds = curr_cell_inds[np.argsort(downsample_conds[curr_cell_inds])]

        curr_high_inds = np.where(downsample_conds[curr_sorted_inds] == 'High')[0]
        curr_low_inds = np.where(downsample_conds[curr_sorted_inds] == 'Low')[0]

        curr_high_order = np.argsort(enriched_cell_guides[curr_sorted_inds[curr_high_inds]])
        curr_low_order = np.argsort(enriched_cell_guides[curr_sorted_inds[curr_low_inds]])

        curr_inds_to_add = np.append(curr_sorted_inds[curr_high_inds[curr_high_order]], curr_sorted_inds[curr_low_inds[curr_low_order]])
        assert(np.all(np.isin(curr_cell_inds, curr_inds_to_add)))
        assert(curr_cell_inds.size == curr_inds_to_add.size)

        refined_cells_order = np.append(refined_cells_order, curr_inds_to_add)
        cells_cluster_labels_refined = np.append(cells_cluster_labels_refined, np.full(curr_cell_inds.size, unique_cluster))

    assert(cells_order.size==refined_cells_order.size)
    assert(cells_cluster_labels.size==cells_cluster_labels_refined.size)
    assert(np.all(np.isin(refined_cells_order, cells_order)))
    
    pd.DataFrame({'Cluster' : cells_cluster_labels_refined, 'Condition' : downsample_conds[refined_cells_order], 'Guide' : enriched_cell_guides[refined_cells_order]}).to_csv(os.path.join(fig_pn, '2C_expression_space_kmeans_%d_cells_LFC_%.2f_refined_order.csv' % (cell_clusters, LFC_val)), index=False)

    features_order, features_cluster_labels = cluster_by_kmeans(downsample_expr.T, features_clusters)
    pd.DataFrame({'Cluster' : features_cluster_labels, 'Covariate' : downsample_features[features_order]}).to_csv(os.path.join(fig_pn, '2C_expression_space_kmeans_%d_features_LFC_%.2f.csv' % (features_clusters, LFC_val)), index=False)

    plot_clustermap(downsample_expr, refined_cells_order, features_order, downsample_conds, downsample_features,
                    os.path.join(fig_pn, '2C_expression_space_kmeans_%d_cells_%d_features_LFC_%.2f_refined_order.pdf' % (cell_clusters, features_clusters, LFC_val)),
                    clim=3, xlabel='Features', ylabel='Cells', cmap_name='seismic', DPI=800)

    # Colorbar for condition
    cond_label_arr = np.full(refined_cells_order.size, 0)
    cond_label_arr[np.where(np.isin(downsample_conds[refined_cells_order], 'Control'))[0]] = 1
    cond_label_arr[np.where(np.isin(downsample_conds[refined_cells_order], 'High'))[0]] = 2
    cond_label_arr[np.where(np.isin(downsample_conds[refined_cells_order], 'Low'))[0]] = 3
    fig, ax = plt.subplots()
    ax_ret = ax.imshow(np.reshape(cond_label_arr, (-1, 1)), aspect=(1/100), cmap='Paired', clim=[0, 3])
    ax.set_xticks([]), ax.set_yticks([])
    fig.savefig(os.path.join(fig_pn, '2C_cond_label_%d_cells_%d_features_%.2f_refined_order.pdf' % (cell_clusters, features_clusters, LFC_val)), bbox_inches='tight', dpi=800)
    plt.close()

if __name__ == "__main__":
    main()
