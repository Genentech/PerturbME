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

plot = True
save = True
cells_per_target_thresh = 17

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

# Load MOI12 sgRNA list
sgRNA_list = np.load(os.path.join(adata_arr_pn, 'adata_parts/sgRNA_names.npy'))
MOI12_inds = np.where( (adata_RNA_CITE.obs['MOI'] == 1) | (adata_RNA_CITE.obs['MOI'] == 2) )[0]

assert(np.where(sgRNA_list[MOI12_inds] == '')[0].size==0)

cells_per_target_high = np.load(os.path.join(full_exp_pn, 'CROP/cells_per_target_high.npy'))
cells_per_target_low = np.load(os.path.join(full_exp_pn, 'CROP/cells_per_target_low.npy'))

assert(cells_per_target_high.size == target_names.size)
assert(cells_per_target_low.size == target_names.size)

target_inds_to_keep = np.where( (cells_per_target_high > cells_per_target_thresh) | (cells_per_target_low > cells_per_target_thresh) )[0]
target_names_to_keep = target_names[target_inds_to_keep]

target_names_to_keep = np.delete(target_names_to_keep, np.where(target_names_to_keep=='NO_SITE')[0])
target_names_to_keep = np.delete(target_names_to_keep, np.where(target_names_to_keep=='ONE_NON-GENE_SITE')[0])

MOI2_count = 0
enriched_combination_count = 0
enriched_singlet_count = 0
target_dict = {}
num_combos = 0
# cov_mod_df = pd.read_csv(os.path.join(save_pn, '17_cells_per_target/beta_corr_covs_7_clusters_list.csv'))
cov_mod_df = pd.read_csv(os.path.join(save_pn, 'grouped_threshold/5_cells_per_target/NMF_U_mat_neg_corr_targets_7_clusters_list.csv'))
module_cov_list = cov_mod_df['Covariate'].to_numpy().astype(np.str)
enriched_mods_count = 0
singlet_mod_count = 0
for curr_sgRNA in sgRNA_list:

    if '--' in curr_sgRNA:
        MOI2_count += 1

        target_1 = curr_sgRNA.split('--')[0].split('_')[0]
        target_2 = curr_sgRNA.split('--')[1].split('_')[0]

        mod_1 = cov_mod_df['Cluster'].loc[cov_mod_df['Covariate']==target_1].to_numpy()
        mod_2 = cov_mod_df['Cluster'].loc[cov_mod_df['Covariate']==target_2].to_numpy()

        if (target_1 in target_names_to_keep) and (target_2 in target_names_to_keep):
            enriched_combination_count += 1

            dict_key = np.sort(np.array([target_1, target_2], dtype=str))[0] + '--' + np.sort(np.array([target_1, target_2], dtype=str))[1]
            if dict_key not in target_dict:
                target_dict[dict_key] = 1

        elif (target_1 in target_names_to_keep) or (target_2 in target_names_to_keep):
            enriched_singlet_count += 1

        if (mod_1.size == 1) and (mod_2.size == 1):
            enriched_mods_count += 1
        elif (mod_1.size == 1) or (mod_2.size == 1):
            singlet_mod_count += 1

bp()
