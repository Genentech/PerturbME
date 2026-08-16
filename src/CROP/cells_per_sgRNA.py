import numpy as np
import pandas as pd
import seaborn as sns
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt

READS_THRESH = 1
READS_PERC_THRESH = 200
# PERC_THRESH = 0.2
PERC_THRESH = 2 # PERC THRESH OFF
UMI_COUNT_THRESH = 10e3
UMI_READS_THRESH = 10e3
FILTERED_UMI_THRESH = 1

exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/20200809_exp/'
figures_pn = os.path.join(exp_pn, 'figures')
sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
CROP_dir = os.path.join(exp_pn, 'CROP')
design_dir = os.path.join(CROP_dir, 'design_mats/')
filtered_RNA_pn = os.path.join(exp_pn, 'filtered_RNA')

# Load sgRNA information
gw_guides_pn = os.path.join(sequencing_info_pn, 'gw_lib_guides.csv')
target_sgRNA_df = pd.read_csv(gw_guides_pn)
target_sgRNA_barcodes = target_sgRNA_df['sgRNA_barcode'].to_numpy().astype(np.str)
target_sgRNA_names =  target_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
assert(target_sgRNA_barcodes.size == target_sgRNA_names.size)
control_sgRNA_df = pd.read_csv(os.path.join(sequencing_info_pn, 'human_TF_CTL_guides.txt'), sep='\t', header=None, names=['sgRNA_name', 'sgRNA_barcode'])
control_sgRNA_barcodes = control_sgRNA_df['sgRNA_barcode'].to_numpy().astype(np.str)
control_sgRNA_names =  control_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
assert(control_sgRNA_barcodes.size == control_sgRNA_names.size)

conds = ['High', 'Low', 'CTL']
all_design_mats = glob.glob(os.path.join(design_dir, '*.design_mat.%d_%d_%.2f_%d_%d_%d.npy' % (READS_THRESH,READS_PERC_THRESH,PERC_THRESH,UMI_COUNT_THRESH,UMI_READS_THRESH,FILTERED_UMI_THRESH)))

plot_df_perc = pd.DataFrame()
plot_df_abs = pd.DataFrame()
CTL_cpg = np.zeros(control_sgRNA_names.size)
HLA_high_cpg = np.zeros(target_sgRNA_names.size)
HLA_low_cpg = np.zeros(target_sgRNA_names.size)
num_CTL_cells, num_high_cells, num_low_cells = 0, 0, 0
for curr_mat in all_design_mats:

	cond = curr_mat.split('/')[-1].split('.')[0]
	print('Processing', cond)

	design_mat = np.load(curr_mat)

	MOI = design_mat.sum(1)
	MOI1_cells = np.where(MOI == 1)[0]

	cells_per_sgRNA = design_mat[MOI1_cells, :].sum(0)

	if 'CTL' in cond:
		assert(design_mat.shape[1] == control_sgRNA_names.size)
		CTL_cpg += cells_per_sgRNA
		num_CTL_cells += MOI1_cells.size
	elif 'High' in cond:
		assert(design_mat.shape[1] == target_sgRNA_names.size)
		HLA_high_cpg += cells_per_sgRNA
		num_high_cells += MOI1_cells.size
	elif 'Low' in cond:
		assert(design_mat.shape[1] == target_sgRNA_names.size)
		HLA_low_cpg += cells_per_sgRNA
		num_low_cells += MOI1_cells.size
	else:
		bp()

bp()

sgRNAs_to_plot = 100
xtick_fontsize = 4
save_fig_pn = os.path.join(figures_pn, 'CROP/Low_cells_per_sgRNA_sorted_%d_sgRNAs.png' % sgRNAs_to_plot)
fig, ax = plt.subplots()
sgRNA_plot_inds = np.flip(np.argsort(HLA_low_cpg))[:sgRNAs_to_plot]
ax.bar(np.arange(sgRNAs_to_plot), HLA_low_cpg[sgRNA_plot_inds])
ax.set_xticks(np.arange(sgRNAs_to_plot)), ax.set_xticklabels(target_sgRNA_names[sgRNA_plot_inds])
ax.set_ylabel('# Cells'), plt.xticks(rotation=90, fontsize=xtick_fontsize)
ax.set_title('HLA Low')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()
save_fig_pn = os.path.join(figures_pn, 'CROP/High_cells_per_sgRNA_sorted_%d_sgRNAs.png' % sgRNAs_to_plot)
fig, ax = plt.subplots()
sgRNA_plot_inds = np.flip(np.argsort(HLA_high_cpg))[:sgRNAs_to_plot]
ax.bar(np.arange(sgRNAs_to_plot), HLA_high_cpg[sgRNA_plot_inds])
ax.set_xticks(np.arange(sgRNAs_to_plot)), ax.set_xticklabels(target_sgRNA_names[sgRNA_plot_inds])
ax.set_ylabel('# Cells'), plt.xticks(rotation=90, fontsize=xtick_fontsize)
ax.set_title('HLA High')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()
save_fig_pn = os.path.join(figures_pn, 'CROP/CTL_cells_per_sgRNA_sorted_%d_sgRNAs.png' % sgRNAs_to_plot)
fig, ax = plt.subplots()
sgRNA_plot_inds = np.flip(np.argsort(CTL_cpg))[:sgRNAs_to_plot]
ax.bar(np.arange(sgRNAs_to_plot), CTL_cpg[sgRNA_plot_inds])
ax.set_xticks(np.arange(sgRNAs_to_plot)), ax.set_xticklabels(control_sgRNA_names[sgRNA_plot_inds])
ax.set_ylabel('# Cells'), plt.xticks(rotation=90, fontsize=xtick_fontsize)
ax.set_title('CTL')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()

# Targets with > 1 sgRNA with > 10 cells
np.sort(target_sgRNA_names[np.where(HLA_low_cpg > 10)[0]])
np.sort(target_sgRNA_names[np.where(HLA_high_cpg > 10)[0]])


save_fig_pn = os.path.join(figures_pn, 'CROP/cells_per_sgRNA_line_plot.png')
fig, ax = plt.subplots()
linewidth=1
ax.plot(Input_cells, label='Input', linewidth=linewidth, marker='None')
ax.plot(HLA_low_cells, label='HLA Low', linewidth=linewidth, marker='None')
ax.plot(HLA_high_cells, label='HLA High', linewidth=linewidth, marker='None')
ax.legend()
ax.set_xlabel('sgRNA'), ax.set_ylabel('# Cells')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()

save_fig_pn = os.path.join(figures_pn, 'CROP/perc_cells_per_sgRNA_line_plot.png')
fig, ax = plt.subplots()
linewidth=1
ax.plot(Input_percs, label='Input', linewidth=linewidth, marker='None')
ax.plot(HLA_high_percs, label='HLA Low', linewidth=linewidth, marker='None')
ax.plot(HLA_low_percs, label='HLA High', linewidth=linewidth, marker='None')
ax.legend()
ax.set_xlabel('sgRNA'), ax.set_ylabel('Percent Cells')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()


save_fig_pn = os.path.join(figures_pn, 'CROP/perc_cells_per_sgRNA_line_plot_sorted_individually.png')
fig, ax = plt.subplots()
sgRNA_order = np.flip(np.argsort(Input_percs))
ax.plot(Input_percs[sgRNA_order], label='Input', linewidth=linewidth, marker='None')
sgRNA_order = np.flip(np.argsort(HLA_high_percs))
ax.plot(HLA_high_percs[sgRNA_order], label='HLA Low', linewidth=linewidth, marker='None')
sgRNA_order = np.flip(np.argsort(HLA_low_percs))
ax.plot(HLA_low_percs[sgRNA_order], label='HLA High', linewidth=linewidth, marker='None')
ax.set_xlabel('sgRNA'), ax.set_ylabel('Percent Cells')
ax.legend()
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()

save_fig_pn = os.path.join(figures_pn, 'CROP/all_channels_cells_per_sgNA.png')
fig, ax = plt.subplots(nrows=3, ncols=1)
nbins=50
xlims = [0, 20]
ax[0].hist(Input_cells, align='left', bins=nbins)
ax[0].set_xlim(xlims), ax[0].set_xticks([])
ax[0].set_ylabel('Input')
ax[1].hist(HLA_high_cells, align='left', bins=nbins)
ax[1].set_xlim(xlims), ax[1].set_xticks([])
ax[1].set_ylabel('HLA High')
ax[2].hist(HLA_low_cells, align='left', bins=nbins)
ax[2].set_xlim(xlims)
ax[2].set_ylabel('HLA Low'), ax[2].set_xlabel('Cells / sgRNA')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()

save_fig_pn = os.path.join(figures_pn, 'CROP/all_channels_perc_cells_per_sgNA.png')
fig, ax = plt.subplots(nrows=3, ncols=1)
nbins=50
xlims = [0, 0.1]
ax[0].hist(Input_percs, align='left', bins=nbins)
ax[0].set_xlim(xlims), ax[0].set_xticks([])
ax[0].set_ylabel('Input')
ax[1].hist(HLA_high_percs, align='left', bins=nbins)
ax[1].set_xlim(xlims), ax[1].set_xticks([])
ax[1].set_ylabel('HLA High')
ax[2].hist(HLA_low_percs, align='left', bins=nbins)
ax[2].set_xlim(xlims)
ax[2].set_ylabel('HLA Low'), ax[2].set_xlabel('Percent Cells / sgRNA')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()


diff_high = np.abs(Input_percs - HLA_high_percs)
diff_low = np.abs(Input_percs - HLA_low_percs)
diff_high_low = np.abs(HLA_high_percs - HLA_low_percs)

sorted_diff_high = np.flip(np.argsort(diff_high))
sorted_diff_low = np.flip(np.argsort(diff_low))
sorted_diff_high_low = np.flip(np.argsort(diff_high_low))

diffs_to_plot = 20

save_fig_pn = os.path.join(figures_pn, 'CROP/enriched_sgRNAs.png')
plot_inds = sorted_diff_high_low[:diffs_to_plot]
abs_diffs = diff_high_low[sorted_diff_high_low[:diffs_to_plot]]
plot_sgRNAs = sgRNA_names[plot_inds]
curr_plot_df = plot_df_perc[plot_df_perc['sgRNA'].isin(plot_sgRNAs)].copy()
curr_plot_df['Difference'] = np.append(abs_diffs, np.append(abs_diffs,abs_diffs))
curr_plot_df.sort_values('Difference', ascending=False, inplace=True)
fig, ax = plt.subplots()
sns.barplot(data=curr_plot_df, x='sgRNA', y='Percent Cells', hue='Channel', ax=ax)
plt.xticks(rotation=-90, fontsize=6)
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()


save_fig_pn = os.path.join(figures_pn, 'CROP/Input_perc_cells_per_sgNA.png')
fig, ax = plt.subplots()
ax.hist(Input_percs, bins=np.unique(Input_percs).size)
ax.set_xlim([0, 0.16])
ax.set_title('Input'), ax.set_xlabel('Percent Cells / sgRNA'), ax.set_ylabel('Count')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()

Input_percs_mean = np.mean(Input_percs)
Input_percs_std = np.std(Input_percs)
HLA_high_z_scores = (HLA_high_percs-Input_percs_mean) / Input_percs_std
HLA_low_z_scores = (HLA_low_percs-Input_percs_mean) / Input_percs_std

plot_order = np.flip(np.argsort(HLA_high_z_scores))
plot_mat = np.column_stack((HLA_high_z_scores[plot_order], HLA_low_z_scores[plot_order]))

save_fig_pn = os.path.join(figures_pn, 'CROP/Zscores_perc_cells_per_sgNA.png')
fig, ax = plt.subplots(figsize=[1*6.4, 1*4.8])
ax_ret = ax.imshow(plot_mat, aspect=(plot_mat.shape[1]/plot_mat.shape[0]), cmap='seismic', clim=[-5, 5])
cbar = fig.colorbar(ax_ret, ax=ax, extend='both'), cbar.minorticks_on()
ax.set_title('Z Scores'), ax.set_ylabel('sgRNA'), ax.set_xlabel('Condition')
ax.set_xticks([0, 1]), ax.set_xticklabels(['HLA_high', 'HLA_low'])
ax.set_ylim([plot_mat.shape[0]-0.5, -0.5]), ax.tick_params(width=0.1)
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()

# ax.set_yticks(np.arange(sgRNA_names.size)), ax.set_yticklabels(sgRNA_names, fontsize=1)

nontarg_inds = np.where(np.char.find(sgRNA_names, 'NO_SITE') == 0)[0]
Input_mean_nontarg = np.mean(Input_percs[nontarg_inds])

save_fig_pn = os.path.join(figures_pn, 'CROP/Input_perc_cells_per_sgNA_nontarg.png')
fig, ax = plt.subplots()
ax.hist(Input_percs[nontarg_inds], bins=np.unique(Input_percs[nontarg_inds]).size)
ax.set_title('Input'), ax.set_xlabel('Percent Cells / sgRNA'), ax.set_ylabel('Count')
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()
ax.set_xlim([0, 0.16])

diff_high_signed = HLA_high_percs-Input_percs
diff_low_signed =  HLA_low_percs-Input_percs

plot_order = np.flip(np.argsort(diff_high_signed))
plot_mat = np.column_stack((diff_high_signed[plot_order], diff_low_signed[plot_order]))

save_fig_pn = os.path.join(figures_pn, 'CROP/Mean_diffs_perc_cells_per_sgNA.png')
fig, ax = plt.subplots(figsize=[1*6.4, 1*4.8])
ax_ret = ax.imshow(plot_mat, aspect=(plot_mat.shape[1]/plot_mat.shape[0]), cmap='seismic', clim=[-0.05, 0.05])
cbar = fig.colorbar(ax_ret, ax=ax, extend='both')
cbar.minorticks_on()
ax.set_ylabel('sgRNA'), ax.set_xlabel('Condition')
ax.set_xticks([0, 1]), ax.set_xticklabels(['HLA_high', 'HLA_low'])
ax.set_ylim([plot_mat.shape[0]-0.5, -0.5]), ax.tick_params(width=0.1)
fig.savefig(save_fig_pn, bbox_inches='tight', dpi=800)
plt.close()
