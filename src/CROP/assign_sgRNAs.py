import glob,os
import pandas as pd
import numpy as np
import pickle
from pdb import set_trace as bp
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp_sparse

READS_THRESH = 1
READS_PERC_THRESH = 200
# PERC_THRESH = 0.2
PERC_THRESH = 2 # PERC THRESH OFF
UMI_COUNT_THRESH = 10e3
UMI_READS_THRESH = 10e3
FILTERED_UMI_THRESH = 2

FIG_DPI = 800
CBC_LENGTH = 16

def main():

	exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	figures_pn = os.path.join(exp_pn, 'figures')
	sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
	CROP_pn = os.path.join(exp_pn, 'CROP')
	filtered_RNA_pn = os.path.join(exp_pn, 'filtered_RNA')
	dict_dir = os.path.join(exp_pn, 'CROP/CBC_UMI_dicts/')
	save_pn = os.path.join(CROP_pn, 'design_mats/')

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

	combined_barcodes = np.append(target_sgRNA_barcodes, control_sgRNA_barcodes) # Targeting channels have both control and targeting sgRNAs

	# cell_per_sgRNA = np.zeros(sgRNA_barcodes.size)
	ALL_MOI = np.array([])
	plot_df = pd.DataFrame()
	CBC_UMI_dicts = glob.glob(os.path.join(dict_dir, '*.pkl'))
	for CBC_UMI_dict in CBC_UMI_dicts:

		channel = CBC_UMI_dict.split('/')[-1].split('CBCUMI')[0][:-1]
		print('Processing channel %s' % channel)

		if 'CTL' in channel:
			sgRNA_barcodes = control_sgRNA_barcodes
		else:
			sgRNA_barcodes = combined_barcodes

		curr_dict = pickle.load(open(CBC_UMI_dict, 'rb'))	# CBC-sgRNA -> UMI reads dictionary for this channel
		filtered_barcs = np.load(os.path.join(filtered_RNA_pn, '%s.barcodes.npy' % channel)) # Filtered CBC barcodes for this channel

		design_mat = np.zeros((filtered_barcs.size, sgRNA_barcodes.size))

		UMI_counts = np.zeros(len(curr_dict))
		for key_i, key in enumerate(curr_dict.keys()):

			total_UMIs = len(curr_dict[key])
			total_UMI_reads = sum([curr_dict[key][UMI_i][1] for UMI_i in range(len(curr_dict[key]))])

			for UMI_i in range(len(curr_dict[key])):
				curr_reads = curr_dict[key][UMI_i][1] # Reads for this UMI
				curr_perc = curr_reads/total_UMI_reads # % of total reads this UMI accounts for

				if curr_reads > READS_THRESH: # Count UMI if more than some threshold
					UMI_counts[key_i] += 1
				elif curr_reads > READS_PERC_THRESH and curr_perc > PERC_THRESH: # Count UMI if large percent + many reads
					UMI_counts[key_i] += 1

			if UMI_counts[key_i] > FILTERED_UMI_THRESH or total_UMIs > UMI_COUNT_THRESH or total_UMI_reads > UMI_READS_THRESH: # Assign sgRNA to cell
				curr_CBC = key[:CBC_LENGTH]
				curr_feature_barc = key[CBC_LENGTH:]

				CBC_ind = np.where(np.char.find(filtered_barcs, curr_CBC) == 0)[0]
				assert(CBC_ind.size < 2)

				if CBC_ind.size == 1: # Not all cells exist (some were filtered)

					feature_ind = np.where(curr_feature_barc == sgRNA_barcodes)[0]

					# assert(feature_ind.size == 1)
					if feature_ind.size == 1: # This library has duplicated sgRNA barcodes, so only assign sgRNAs with unique barcode
						# cell_per_sgRNA[feature_ind] += 1		# Count cell for this sgRNA
						design_mat[CBC_ind, feature_ind] = 1 	# Assign sgRNA to this cell

			if (key_i+1) % 10e3 == 0:
				print('%.2f%% complete (%d entries of %d total entries)' % (100*((key_i+1)/len(curr_dict)),key_i+1,len(curr_dict)))

		# np.save(os.path.join(save_pn, channel+'.design_mat.%d_%d_%.2f_%d_%d_%d.npy' % (READS_THRESH,READS_PERC_THRESH,PERC_THRESH,UMI_COUNT_THRESH,UMI_READS_THRESH,FILTERED_UMI_THRESH)), design_mat) # Dense
		sp_sparse.save_npz(os.path.join(save_pn, channel+'.sgRNA_design_mat.%d_%d_%.2f_%d_%d_%d.npz' % (READS_THRESH,READS_PERC_THRESH,PERC_THRESH,UMI_COUNT_THRESH,UMI_READS_THRESH,FILTERED_UMI_THRESH)), sp_sparse.csr_matrix(design_mat)) # Sparse

		curr_MOI = design_mat.sum(1)

		plot_df = plot_df.append(pd.DataFrame({'Sample' : [channel], 'MOI' : ['0'], 'Percent Cells' : [100*(np.where(curr_MOI == 0)[0].size/curr_MOI.size)]}), ignore_index=True)
		plot_df = plot_df.append(pd.DataFrame({'Sample' : [channel], 'MOI' : ['1'], 'Percent Cells' : [100*(np.where(curr_MOI == 1)[0].size/curr_MOI.size)]}), ignore_index=True)
		plot_df = plot_df.append(pd.DataFrame({'Sample' : [channel], 'MOI' : ['2'], 'Percent Cells' : [100*(np.where(curr_MOI == 2)[0].size/curr_MOI.size)]}), ignore_index=True)

		ALL_MOI = np.append(ALL_MOI, curr_MOI)

		print('MOI = 0 : %d total cells (%.2f%% of all cells)' % (np.where(ALL_MOI == 0)[0].size, 100*(np.where(ALL_MOI == 0)[0].size/ALL_MOI.size)))
		print('MOI = 1 : %d total cells (%.2f%% of all cells)' % (np.where(ALL_MOI == 1)[0].size, 100*(np.where(ALL_MOI == 1)[0].size/ALL_MOI.size)))
		print('MOI = 2 : %d total cells (%.2f%% of all cells)' % (np.where(ALL_MOI == 2)[0].size, 100*(np.where(ALL_MOI == 2)[0].size/ALL_MOI.size)))

	bp()

	plot_df = plot_df.sort_values(by=['Sample'])
	fig_save_pn = os.path.join(figures_pn, 'MOI_by_channel.png')
	fig, ax = plt.subplots()
	sns.barplot(data=plot_df, x='Sample', y='Percent Cells', hue='MOI', ax=ax)
	plt.xticks(rotation=90), ax.set_title('MOI')
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
	plt.close()

	save_fig_pn = os.path.join(figures_pn, 'ALL_MOI_%d_%d_%.2f_%d_%d_%d.png' % (READS_THRESH,READS_PERC_THRESH,PERC_THRESH,UMI_COUNT_THRESH,UMI_READS_THRESH,FILTERED_UMI_THRESH))
	fig, ax = plt.subplots()
	ax.hist(ALL_MOI, bins=np.max(ALL_MOI).astype(int), align='left')
	ax.set_xlabel('MOI'), ax.set_ylabel('Count')
	fig.savefig(save_fig_pn, bbox_inches='tight', dpi=FIG_DPI)
	plt.close()

	# save_fig_pn = os.path.join(figures_pn, 'cells_per_sgRNA_%d_%d_%.2f_%d_%d_%d.png' % (READS_THRESH,READS_PERC_THRESH,PERC_THRESH,UMI_COUNT_THRESH,UMI_READS_THRESH,FILTERED_UMI_THRESH))
	# fig, ax = plt.subplots()
	# ax.bar(np.arange(cell_per_sgRNA.size),cell_per_sgRNA)
	# ax.set_xlabel('sgRNA'), ax.set_ylabel('# Cells')
	# fig.savefig(save_fig_pn, bbox_inches='tight', dpi=FIG_DPI)
	# plt.close()


if __name__ == "__main__":
    main()
