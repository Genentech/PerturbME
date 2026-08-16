import os
from pdb import set_trace as bp
import numpy as np
import pandas as pd
import glob

exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/20200809_exp/'
figures_pn = os.path.join(exp_pn, 'figures')
sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
gDNA_pn = os.path.join(exp_pn, 'guide_seq')

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

sgRNA_names = np.append(target_sgRNA_names, control_sgRNA_names)
sgRNA_barcodes = np.append(target_sgRNA_barcodes, control_sgRNA_barcodes)
assert(sgRNA_names.size == sgRNA_barcodes.size)

save_pn = os.path.join(gDNA_pn, 'mageck_ALL_sgRNAs.csv')

fp = open(save_pn,'w')
for sgRNA_i,sgRNA in enumerate(sgRNA_names):

	if 'NO_SITE' in sgRNA:
		curr_target = 'NO_SITE'
	elif 'ONE_NON-GENE_SITE' in sgRNA:
		curr_target = 'ONE_NON-GENE_SITE'
	else:
		curr_target = sgRNA.split('_')[0]

	_=fp.write('%s,%s,%s\n' % (sgRNA,sgRNA_barcodes[sgRNA_i],curr_target))

fp.close()

save_pn = os.path.join(gDNA_pn, 'control_sgRNAs.txt')
fp = open(save_pn,'w')
for sgRNA in sgRNA_names:
	if 'NO_SITE' in sgRNA or 'ONE_NON-GENE_SITE' in sgRNA:
		_=fp.write('%s\n' % sgRNA)

fp.close()
