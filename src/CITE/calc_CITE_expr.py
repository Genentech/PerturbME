import os
from pdb import set_trace as bp
import numpy as np
import pandas as pd
import glob

# Function to make dictionary mapping antibody names to their IgG control antibody
def make_control_dict(CITE_control_gene_df):
	ret_dict = {}
	for antibody in CITE_control_gene_df['Antibody'].to_numpy():
		ctrl_antibody = CITE_control_gene_df.loc[CITE_control_gene_df['Antibody']==antibody]['Control'].to_numpy()[0]
		ret_dict[antibody] = ctrl_antibody
	return ret_dict

exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
fig_pn = os.path.join(exp_pn, 'figures')
barc_txt_path = os.path.join(exp_pn, 'terra/Barcodes')
CITE_path = os.path.join(exp_pn, 'CITE')
terra_CITE_path = os.path.join(CITE_path, 'terra_csvs')
save_pn = os.path.join(CITE_path, 'CITE_expr')

CITE_control_gene_df = pd.read_csv(os.path.join(exp_pn, 'sequencing_info/antibody_to_control.csv'))
control_map_dict = make_control_dict(CITE_control_gene_df)
all_txts = glob.glob(os.path.join(barc_txt_path, '*.txt'))

antibody_names_pn = os.path.join(save_pn, 'antibody_names.npy')

if os.path.exists(antibody_names_pn):
	antibody_names = np.load(antibody_names_pn)
else:
	antibody_names = CITE_control_gene_df['Antibody'].to_numpy().astype(np.str)
	np.save(antibody_names_pn, antibody_names.astype(np.str))

ctrl_antibody_names = np.unique(CITE_control_gene_df['Control'].to_numpy()).astype(np.str)
num_antibodies = antibody_names.size
num_ctrl_antibodies = ctrl_antibody_names.size
for txt in all_txts:

	channel = txt.split('/')[-1].split('.')[0]
	print('Processing %s' % channel)

	curr_save_path = os.path.join(save_pn, '%s.CITE_norm_expr.npy' % channel)

	if not os.path.exists(curr_save_path):
		
		barcodes_rna = np.loadtxt(txt, dtype=str)
		np.save(os.path.join(save_pn, '%s.CITE_barcodes.npy' % channel), barcodes_rna)
		CITE_csv = pd.read_csv(os.path.join(terra_CITE_path, channel+'.CITE.csv'))
		expr_labels = CITE_csv['Antibody'].to_numpy()

		assert(CITE_csv.shape[0] == (num_antibodies + num_ctrl_antibodies))

		all_CITE_expr = np.zeros((num_antibodies, barcodes_rna.size))
		for barcode_i, rna_barcode in enumerate(barcodes_rna):

			if rna_barcode in CITE_csv.columns:

				# Antibody expression output by cumulus workflow
				curr_expr = CITE_csv[rna_barcode].to_numpy()

				assert(curr_expr.size == expr_labels.size)

				# Normalize to control - formula max(0, log( (expr+1) / (ctrl_expr+1) ) )
				norm_expr = np.zeros(num_antibodies)
				for antibody_i, antibody in enumerate(antibody_names):
					targ_ind = np.where(expr_labels == antibody)[0]
					assert(targ_ind.size == 1)
					ctrl_ind = np.where(expr_labels == control_map_dict[antibody])[0]
					assert(ctrl_ind.size == 1)

					log_expr = np.log( (curr_expr[targ_ind]+1) / (curr_expr[ctrl_ind]+1) )
					norm_expr[antibody_i] = np.max([0, log_expr])

				all_CITE_expr[:, barcode_i] = norm_expr

			else:
				print('Failed to find rna barcode in CITE csv for channel %s' % channel)

		np.save(curr_save_path, all_CITE_expr)
