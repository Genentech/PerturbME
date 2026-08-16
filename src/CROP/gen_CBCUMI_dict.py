import gzip
import glob,os
import pandas as pd
import numpy as np
import pickle
from pdb import set_trace as bp
import sys
# Conda environment pm1

def editDistance(str1, str2):
    mismatch = 0
    for char_i in range(len(str1)):
        if str1[char_i] != str2[char_i]:
            mismatch += 1
    return mismatch

def find_feature_barc(R2_seq, sgRNA_barcodes, sgRNA_dict, mismatch=3):

    if mismatch == 0:
        return R2_seq in sgRNA_dict

    found = False
    for sgRNA_barcode in sgRNA_barcodes:
        if editDistance(R2_seq, sgRNA_barcode) <= mismatch:
            if found:
                return False # sgRNA collision - can't confidenty assign sgRNA
            found = True

    return found

def make_dict_from_arr(np_arr):
    ret_dict = {}
    for entry in np_arr:
        ret_dict[entry] = 1
    return ret_dict

CBC_LENGTH = 16
UMI_LENGTH = 12
COORD_WIDTH = 5 # Longest possible string defining coordinate in fastq header
MISMATCH = 0

def main():

	exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	fig_pn = os.path.join(exp_pn, 'figures')
	terra_pn = os.path.join(exp_pn, 'terra')
	sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
	fastq_dir1 = os.path.join(exp_pn, 'CROP_DO_fastqs_run1')
	fastq_dir2 = os.path.join(exp_pn, 'CROP_DO_fastqs_run2')
	barcode_txt_pn = os.path.join(terra_pn, 'Barcodes')

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

	target_sgRNA_dict = make_dict_from_arr(combined_barcodes)
	control_sgRNA_dict = make_dict_from_arr(control_sgRNA_barcodes)

	assert(len(control_sgRNA_barcodes[0]) == len(target_sgRNA_barcodes[0]))
	sgRNA_barcode_length = len(target_sgRNA_barcodes[0])

	R2_fps_1 = glob.glob(os.path.join(fastq_dir1, '*_R2_*'))
	R2_fps_2 = glob.glob(os.path.join(fastq_dir2, '*_R2_*'))
	R2_fps = np.append(R2_fps_1, R2_fps_2)
	R2_fps = [fp for fp in R2_fps if 'Undetermined' not in fp]

	for R2_fp in R2_fps:

		sample_prefix = R2_fp.split('/')[-1].split('_')[0]
		sample_suffix = R2_fp.split('/')[-1].split('S')[0].split('_')[1]
		channel = sample_prefix+'_'+sample_suffix
		print('Processing channel %s' % channel)

		if sample_prefix == 'CTL':
			sgRNA_dict = control_sgRNA_dict
			sgRNA_barcodes = control_sgRNA_barcodes
		else:
			sgRNA_dict = target_sgRNA_dict
			sgRNA_barcodes = combined_barcodes

		CBC_UMI_dict_savepath = os.path.join(exp_pn, 'CROP/CBC_UMI_dicts/%s_CBCUMI_map_mismatch%d.pkl' % (channel, MISMATCH))
		if os.path.exists(CBC_UMI_dict_savepath): # Dict already exists, load for appending and then delete
			CBC_UMI_dict = pickle.load(open(CBC_UMI_dict_savepath, 'rb'))
			os.system('rm %s' % CBC_UMI_dict_savepath)
			print('Loaded and deleted %s' % CBC_UMI_dict_savepath)
		else:
			CBC_UMI_dict = {} # Initialize new dictionary

		RNA_barcodes = np.loadtxt(os.path.join(barcode_txt_pn, channel+'.txt'), dtype=str)
		RNA_dict = make_dict_from_arr(RNA_barcodes)

		open_R1_fp = gzip.open(R2_fp.replace('R2', 'R1'))
		open_R2_fp = gzip.open(R2_fp)

		loop_count = 0
		while True:

		    R1_head = open_R1_fp.readline().decode().strip() # Header
		    R1_seq = open_R1_fp.readline().decode().strip() # Sequence
		    R1_line3 = open_R1_fp.readline().decode().strip()
		    R1_quality = open_R1_fp.readline().decode().strip()

		    R2_head = open_R2_fp.readline().decode().strip() # Header
		    R2_seq = open_R2_fp.readline().decode().strip() # Sequence
		    R2_line3 = open_R2_fp.readline().decode().strip()
		    R2_quality = open_R2_fp.readline().decode().strip()

		    if R1_head == '': # End of file
		        assert(R2_head == '')
		        break

		    # Ensure R1 and R2 align
		    R1_coord = int(R1_head.split(':')[4].rjust(COORD_WIDTH, '0') + R1_head.split(':')[5].rjust(COORD_WIDTH, '0') + R1_head.split(':')[6].split(' ')[0].rjust(COORD_WIDTH, '0')) # Parse header
		    R2_coord = int(R2_head.split(':')[4].rjust(COORD_WIDTH, '0') + R2_head.split(':')[5].rjust(COORD_WIDTH, '0') + R2_head.split(':')[6].split(' ')[0].rjust(COORD_WIDTH, '0')) # Parse header
		    assert(R1_coord==R2_coord)

		    # Extract info from read sequences
		    tenx_CBC = R1_seq[:CBC_LENGTH]
		    tenx_UMI = R1_seq[CBC_LENGTH:]
		    assert(len(tenx_UMI) == UMI_LENGTH)

		    R2_feature_barc = R2_seq[:sgRNA_barcode_length]

		    # Add relevant elements to dictionary
		    if tenx_CBC in RNA_dict: # Exact match for CBC barcode

		        if find_feature_barc(R2_feature_barc, sgRNA_barcodes, sgRNA_dict, mismatch=MISMATCH):

		            CBC_UMI_dict_key = tenx_CBC + R2_feature_barc # Key for dictionary

		            if CBC_UMI_dict_key in CBC_UMI_dict: # CBC-sgRNA pair already exists
		                found_UMI = False
		                for UMI_i in range(len(CBC_UMI_dict[CBC_UMI_dict_key])):
		                    if tenx_UMI in CBC_UMI_dict[CBC_UMI_dict_key][UMI_i]: # UMI already found
		                        CBC_UMI_dict[CBC_UMI_dict_key][UMI_i][1] += 1 # So increment UMI count
		                        found_UMI = True
		                if not found_UMI: # New UMI
		                    CBC_UMI_dict[CBC_UMI_dict_key].append([tenx_UMI, 1]) # Add new UMI to list
		            else: # New CBC-sgRNA pair, create new entry in dictionary
		                CBC_UMI_dict[CBC_UMI_dict_key] = [[tenx_UMI, 1]]

		    if loop_count % 10e6 == 0:
		        print('Processed %d lines' % loop_count)

		    loop_count += 1

		CBC_UMI_fp = open(CBC_UMI_dict_savepath, 'wb') # Save dictionary (override if already exists)
		pickle.dump(CBC_UMI_dict, CBC_UMI_fp)
		print('Saved dictionary %s with size %d bytes' % (CBC_UMI_dict_savepath,sys.getsizeof(CBC_UMI_dict)))

		# Close all file pointers
		open_R1_fp.close()
		open_R2_fp.close()
		CBC_UMI_fp.close()

if __name__ == "__main__":
    main()
