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
    fig_pn = os.path.join(full_exp_pn, 'figures')
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
    highE_lowE_df = pd.read_csv(os.path.join(DE_pn, 'RNA_DE_results_MAST_highE_lowE.csv'))
    if False:
        highE_lowE_features = highE_lowE_df['Unnamed: 0'].loc[highE_lowE_df['p_val_adj'] < FDR_val].to_numpy().astype(np.str)
    else:
    	LFC_val = 0.2
    	highE_lowE_features = highE_lowE_df['Unnamed: 0'].loc[ (highE_lowE_df['p_val_adj'] < FDR_val) & (np.abs(highE_lowE_df['avg_logFC']) > LFC_val)].to_numpy().astype(np.str)

    HLA_gene_inds = np.where(np.isin(features, highE_lowE_features))[0]

    if False:
        sc.pp.highly_variable_genes(adata_RNA_CITE, n_top_genes=2000)
        highly_variable_gene_inds = np.where(adata_RNA_CITE.var['highly_variable'])[0]
        good_features = np.union1d(highly_variable_gene_inds, HLA_gene_inds)
    else:
        good_features = HLA_gene_inds

    downsample_expr = expr_mat[cell_inds, :][:, good_features]
    downsample_conds = cond_arr[cell_inds]
    downsample_features = features[good_features]

    if True:
        downsample_expr = sp_stats.zscore(downsample_expr, axis=0)

    fastcluster_rows, fastcluster_cols = cluster_mat_fastcluster(downsample_expr)

    cond_label_arr = np.full(fastcluster_rows.size, 0)
    cond_label_arr[np.where(np.isin(downsample_conds[fastcluster_rows], 'Control'))[0]] = 1
    cond_label_arr[np.where(np.isin(downsample_conds[fastcluster_rows], 'High'))[0]] = 2
    cond_label_arr[np.where(np.isin(downsample_conds[fastcluster_rows], 'Low'))[0]] = 3
    fig, ax = plt.subplots()
    ax_ret = ax.imshow(np.reshape(cond_label_arr, (-1, 1)), aspect=(1/100), cmap='Paired', clim=[0, 3])
    ax.set_xticks([]), ax.set_yticks([])
    fig.savefig(os.path.join(fig_pn, 'Clustering/expression/Cond_label_fastcluster.png'), bbox_inches='tight', dpi=800)
    plt.close()

    all_cell_clusters = np.arange(8, 12, 1)
    # cell_clusters = 3

    for cell_clusters in all_cell_clusters:
        cells_order, cells_cluster_labels = cluster_by_kmeans(downsample_expr, cell_clusters)
        pd.DataFrame({'Cluster' : cells_cluster_labels, 'Condition' : downsample_conds[cells_order], 'Guide' : enriched_cell_guides[cells_order]}).to_csv(os.path.join(fig_pn, 'Clustering/expression/expression_space_kmeans_%d_cells_LFC_%.2f.csv' % (cell_clusters, LFC_val)), index=False)

        all_features_clusters = np.arange(4, 8, 1)
        # features_clusters = 5

        for features_clusters in all_features_clusters:

            features_order, features_cluster_labels = cluster_by_kmeans(downsample_expr.T, features_clusters)
            pd.DataFrame({'Cluster' : features_cluster_labels, 'Covariate' : downsample_features[features_order]}).to_csv(os.path.join(fig_pn, 'Clustering/expression/expression_space_kmeans_%d_features_LFC_%.2f.csv' % (features_clusters, LFC_val)), index=False)

            plot_clustermap(downsample_expr, cells_order, features_order, downsample_conds, downsample_features,
                            os.path.join(fig_pn, 'Clustering/expression/expression_space_kmeans_%d_cells_%d_features_LFC_%.2f.png' % (cell_clusters, features_clusters, LFC_val)),
                            clim=3, xlabel='Features', ylabel='Cells', cmap_name='seismic', DPI=300)

            if cell_clusters==1:
                plot_clustermap(downsample_expr, fastcluster_rows, features_order, downsample_conds, downsample_features,
                				os.path.join(fig_pn, 'Clustering/expression/expression_space_fastcluster_%d_features.png' % features_clusters),
                				clim=3, xlabel='Features', ylabel='Cells', cmap_name='seismic')


            # Colorbar for HLA program
            # HLA_label_arr = np.full(features_order.size, 0)
            # HLA_label_arr[np.where(np.isin(downsample_features[features_order], highE_lowE_features))[0]] = 1
            # fig, ax = plt.subplots()
            # ax_ret = ax.imshow(np.reshape(HLA_label_arr, (-1, 1)), aspect=(1/100), cmap='Oranges', clim=[0, 1.4])
            # ax.set_xticks([]), ax.set_yticks([])
            # fig.savefig(os.path.join(fig_pn, 'Clustering/expression/HLA_label_%d_cells_%d_features_LFC_%.2f.png' % (cell_clusters, features_clusters, LFC_val)), bbox_inches='tight', dpi=800)
            # plt.close()

            # Colorbar for condition
            cond_label_arr = np.full(cells_order.size, 0)
            cond_label_arr[np.where(np.isin(downsample_conds[cells_order], 'Control'))[0]] = 1
            cond_label_arr[np.where(np.isin(downsample_conds[cells_order], 'High'))[0]] = 2
            cond_label_arr[np.where(np.isin(downsample_conds[cells_order], 'Low'))[0]] = 3
            fig, ax = plt.subplots()
            ax_ret = ax.imshow(np.reshape(cond_label_arr, (-1, 1)), aspect=(1/100), cmap='Paired', clim=[0, 3])
            ax.set_xticks([]), ax.set_yticks([])
            fig.savefig(os.path.join(fig_pn, 'Clustering/expression/Cond_label_%d_cells_%d_features_%.2f.png' % (cell_clusters, features_clusters, LFC_val)), bbox_inches='tight', dpi=800)
            plt.close()

if __name__ == "__main__":
    main()
