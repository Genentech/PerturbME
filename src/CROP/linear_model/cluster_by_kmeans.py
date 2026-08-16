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
					label_ticks=False, fig_size_mult=1, title=''):

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

    cells_per_target_thresh = 17
    plot = True
    save = True

    full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
    fig_pn = os.path.join(full_exp_pn, 'figures/CROP/linear_model')
    adata_arr_pn = os.path.join(full_exp_pn, 'adata_arrs')
    sequencing_info_pn = os.path.join(full_exp_pn, 'sequencing_info')
    CROP_pn = os.path.join(full_exp_pn, 'CROP')
    design_mat_pn = os.path.join(CROP_pn, 'design_mats')
    load_pn = os.path.join(CROP_pn, 'linear_model')
    DE_pn = os.path.join(full_exp_pn, 'DE/DE_csvs')

    target_names = np.load(os.path.join(design_mat_pn, 'target_names.npy'))
    cells_per_target_high = np.load(os.path.join(CROP_pn, 'cells_per_target_high.npy'))
    cells_per_target_low = np.load(os.path.join(CROP_pn, 'cells_per_target_low.npy'))

    cov_names = np.load(os.path.join(load_pn, '%s_cells_per_target/cov_names.npy' % cells_per_target_thresh))
    feature_names = np.load(os.path.join(load_pn, '%s_cells_per_target/feature_names.npy' % cells_per_target_thresh))
    B_mat = np.load(os.path.join(load_pn, '%s_cells_per_target/EN_B_EM.npy' % cells_per_target_thresh))

    assert(B_mat.shape[0] == feature_names.size)
    assert(B_mat.shape[1] == cov_names.size)

    filtered_B, filtered_features, filtered_covs = remove_sparse_rows_and_columns(B_mat, feature_names, cov_names, approx_zero = 0, col_sparse_thresh = 0.95, row_sparse_thresh = 0.9)

    B_correlation_covs = calc_corr_matrix(filtered_B)
    B_correlation_features = calc_corr_matrix(filtered_B.T)

    # Calculate inertia to determine number of clusters
    if False:
        K_covs, inertia_covs = calculate_inertia(B_correlation_covs, max_clusters=20)
        K_features, inertia_features = calculate_inertia(B_correlation_features, max_clusters=20)

        fig, ax = plt.subplots()
        ax.plot(K_covs, inertia_covs, marker='.')
        ax.set_xlabel('# Clusters'), ax.set_ylabel('Inertia'), ax.set_title('Covariates')
        ax.set_xticks(np.arange(0, 22, 2))
        fig.savefig(os.path.join(fig_pn, '%s_cells_per_target_find_k_covs.png' % cells_per_target_thresh), bbox_inches='tight', dpi=800)
        plt.close()
        fig, ax = plt.subplots()
        ax.plot(K_features, inertia_features, marker='.')
        ax.set_xlabel('# Clusters'), ax.set_ylabel('Inertia'), ax.set_title('Features')
        ax.set_xticks(np.arange(0, 22, 2))
        fig.savefig(os.path.join(fig_pn, '%s_cells_per_target_find_k_features.png' % cells_per_target_thresh), bbox_inches='tight', dpi=800)
        plt.close()



    num_clusters_covs = 7
    num_clusters_features = 5

    covs_order, covs_cluster_labels = cluster_by_kmeans(B_correlation_covs, num_clusters_covs)
    features_order, features_cluster_labels = cluster_by_kmeans(B_correlation_features, num_clusters_features)
    bp()
    pd.DataFrame({'Cluster' : covs_cluster_labels, 'Covariate' : filtered_covs[covs_order]}).to_csv(os.path.join(load_pn, 'beta_corr_covs_%d_clusters_list.csv' % (num_clusters_covs)), index=False)
    pd.DataFrame({'Cluster' : features_cluster_labels, 'Feature' : filtered_features[features_order]}).to_csv(os.path.join(load_pn, 'beta_corr_features_%d_clusters_list.csv' % (num_clusters_features)), index=False)
    if plot:
        plot_clustermap(filtered_B, features_order, covs_order, filtered_features, filtered_covs,
                        os.path.join(fig_pn, '%s_cells_per_target/beta_%d_features_%d_covs.png' % (cells_per_target_thresh, num_clusters_features, num_clusters_covs)),
                        clim=0.3, xlabel='Covariates', ylabel='Features', label_ticks=False, cmap_name='bwr')

    save_pn = os.path.join(load_pn, '%s_cells_per_target/GO_analysis/' % cells_per_target_thresh)

    if True:
        background_list_save_fp = open(os.path.join(save_pn, 'background_list.txt'), 'w+')
        _=[background_list_save_fp.write('%s\n' % element) for element in filtered_features]
        background_list_save_fp.close()

        for feature_cluster_num in np.arange(num_clusters_features):

            features_in_this_cluster = filtered_features[features_order[np.where(features_cluster_labels==feature_cluster_num)[0]]]

            save_fp = open(os.path.join(save_pn, '%s_feature_cluster.txt' % feature_cluster_num), 'w+')
            _=[save_fp.write('%s\n' % element) for element in features_in_this_cluster]
            save_fp.close()


    bp()
    # Colorbar for HLA program
    FDR_val = 0.25
    highE_lowE_df = pd.read_csv(os.path.join(DE_pn, 'RNA_DE_results_MAST_highE_lowE.csv'))
    highE_lowE_features = highE_lowE_df['Unnamed: 0'].loc[highE_lowE_df['p_val_adj'] < FDR_val].to_numpy().astype(np.str)

    HLA_label_arr = np.full(features_order.size, 0)
    HLA_label_arr[np.where(np.isin(filtered_features[features_order], highE_lowE_features))[0]] = 1

    fig, ax = plt.subplots()
    ax_ret = ax.imshow(np.reshape(HLA_label_arr, (-1, 1)), aspect=(1/100), cmap='Oranges', clim=[0, 1.4])
    ax.set_xticks([]), ax.set_yticks([])
    fig.savefig(os.path.join(fig_pn, '%s_cells_per_target/HLA_label_%d_features.png' % (cells_per_target_thresh, num_clusters_features)), bbox_inches='tight', dpi=800)
    plt.close()


    # Figure out how many clusters for covariates and features

    cov_clusters_to_test = np.arange(2, 11, 1)

    for num_clusters_covs in cov_clusters_to_test:

        # num_clusters_covs = 6
        covs_order, covs_cluster_labels = cluster_by_kmeans(B_correlation_covs, num_clusters_covs)
        if plot:
            plot_clustermap(B_correlation_covs, covs_order, covs_order, filtered_covs, filtered_covs,
            				os.path.join(fig_pn, '%s_cells_per_target/%s_cells_per_target_beta_corr_covs_%d_clusters.png' % (cells_per_target_thresh, cells_per_target_thresh, num_clusters_covs)),
            				clim=0.2, xlabel='Covariate', ylabel='Covariate', label_ticks=True, cmap_name='PRGn')

        if plot:
            fig, ax = plt.subplots()
            ax_ret = ax.imshow(np.reshape(covs_cluster_labels, (-1, 1)), aspect=(1/100), cmap='Paired', clim=[0, num_clusters_covs])
            ax.set_xticks([]), ax.set_yticks([])
            fig.savefig(os.path.join(fig_pn, '%s_cells_per_target/%s_cells_per_target_beta_corr_covs_%d_clusters_label.png' % (cells_per_target_thresh, cells_per_target_thresh, num_clusters_covs)), bbox_inches='tight', dpi=800)
            plt.close()

        # Add colorbar to show number of cells in high / low channels per target
        cells_high_label_arr = np.full(filtered_covs.size, 0)
        cells_low_label_arr = np.full(filtered_covs.size, 0)
        for cov_i, curr_cov in enumerate(filtered_covs[covs_order]):
            if curr_cov in target_names:
                curr_target_ind = np.where(curr_cov == target_names)[0]
                assert(curr_target_ind.size==1)
                cells_high_label_arr[cov_i] = cells_per_target_high[curr_target_ind]
                cells_low_label_arr[cov_i] = cells_per_target_low[curr_target_ind]

        if plot:
            fig, ax = plt.subplots()
            ax_ret = ax.imshow(np.reshape(cells_high_label_arr, (-1, 1)), aspect=(1/100), cmap='Reds', clim=[0, 30])
            ax.set_xticks([]), ax.set_yticks([])
            fig.savefig(os.path.join(fig_pn, '%s_cells_per_target/%s_cells_per_target_beta_corr_covs_%d_clusters_high_label_arr.png' % (cells_per_target_thresh, cells_per_target_thresh, num_clusters_covs)), bbox_inches='tight', dpi=800)
            plt.close()

            fig, ax = plt.subplots()
            ax_ret = ax.imshow(np.reshape(cells_low_label_arr, (-1, 1)), aspect=(1/100), cmap='Blues', clim=[0, 30])
            ax.set_xticks([]), ax.set_yticks([])
            fig.savefig(os.path.join(fig_pn, '%s_cells_per_target/%s_cells_per_target_beta_corr_covs_%d_clusters_low_label_arr.png' % (cells_per_target_thresh, cells_per_target_thresh, num_clusters_covs)), bbox_inches='tight', dpi=800)
            plt.close()

    if plot:
        fig, ax = plt.subplots(figsize=(1, 12))
        fig.subplots_adjust(bottom=0.5)
        norm = matplotlib.colors.Normalize(vmin=0, vmax=30)
        cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=plt.cm.Reds, norm=norm, orientation='vertical')
        cb1.set_label('Cells / Target High', fontsize=16, fontweight='bold')
        fig.savefig(os.path.join(fig_pn, '%s_cells_per_target/covs_%d_clusters_high_label_arr_cmap.png' % (cells_per_target_thresh, num_clusters_covs)), bbox_inches='tight', dpi=800)
        plt.close()

        fig, ax = plt.subplots(figsize=(1, 12))
        fig.subplots_adjust(bottom=0.5)
        norm = matplotlib.colors.Normalize(vmin=0, vmax=30)
        cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=plt.cm.Blues, norm=norm, orientation='vertical')
        cb1.set_label('Cells / Target Low', fontsize=16, fontweight='bold')
        fig.savefig(os.path.join(fig_pn, '%s_cells_per_target/covs_%d_clusters_low_label_arr_cmap.png' % (cells_per_target_thresh, num_clusters_covs)), bbox_inches='tight', dpi=800)
        plt.close()

    feature_clusters_to_test = np.arange(2, 11, 1)
    for num_clusters_features in feature_clusters_to_test:

        # num_clusters_features = 6
        features_order, features_cluster_labels = cluster_by_kmeans(B_correlation_features, num_clusters_features)
        if plot:
            plot_clustermap(B_correlation_features, features_order, features_order, filtered_features, filtered_features,
            				os.path.join(fig_pn, '%s_cells_per_target/%s_cells_per_target_beta_corr_features_%d_clusters.png' % (cells_per_target_thresh, cells_per_target_thresh, num_clusters_features)),
            				clim=0.2, xlabel='Target', ylabel='Target', label_ticks=False, cmap_name='PRGn')

        if plot:
            fig, ax = plt.subplots()
            ax_ret = ax.imshow(np.reshape(features_cluster_labels, (-1, 1)), aspect=(1/100), cmap='Paired', clim=[0, num_clusters_features])
            ax.set_xticks([]), ax.set_yticks([])
            fig.savefig(os.path.join(fig_pn, '%s_cells_per_target/%s_cells_per_target_beta_corr_features_%d_clusters_label.png' % (cells_per_target_thresh, cells_per_target_thresh, num_clusters_features)), bbox_inches='tight', dpi=800)
            plt.close()

if __name__ == "__main__":
    main()
