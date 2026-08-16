import os
from pdb import set_trace as bp
import numpy as np

gb_address = 'gs://fc-38170c52-8a3c-4612-9433-416a1cfe147e'

rna_out_pn = os.path.join(gb_address, '202008_terra_output/terra_output_RNA/')
save_path = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/terra/terra_h5s'

all_channels = np.array(['CTL_A1a', 'CTL_A1b', 'High_A1a', 'High_A1b', 'High_A2a', 'High_A2b', 'High_A3a',
						 'High_B1a', 'High_B2a', 'High_B2b', 'High_B3a', 'Low_A1a', 'Low_A21', 'Low_A2a', 'Low_A2b',
						 'Low_B1a', 'Low_B21', 'Low_B2a', 'Low_B2b', 'Low_B3a'])

for channel in all_channels:
	source_str = os.path.join(rna_out_pn, os.path.join(channel,'filtered_feature_bc_matrix.h5'))
	dest_str = os.path.join(save_path, channel+'.h5')

	if not os.path.exists(dest_str):
		os.system('gsutil -m cp -r %s %s' % (source_str, dest_str))
