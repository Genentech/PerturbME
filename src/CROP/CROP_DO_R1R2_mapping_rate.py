import gzip
import glob,os
import pandas as pd
import numpy as np
import seaborn as sns
from pdb import set_trace as bp
import matplotlib.pyplot as plt
from tqdm import tqdm

# Number sgRNA names (sgRNAs with same target have same name in .txt file)
def number_sgRNA_names(guide_names):
	guide_dict = {}

	sgRNA_names = np.array([])
	for guide in tqdm(guide_names):

		if guide in guide_dict:
			list_len = len(guide_dict[guide])
			sgRNA_names = np.append(sgRNA_names, guide + '_' + str(int(list_len+1)))
			guide_dict[guide].append(guide)
		else:
			sgRNA_names = np.append(sgRNA_names, guide + '_1')
			guide_dict[guide] = [guide]

	assert(sgRNA_names.size == guide_names.size)
	return sgRNA_names

def make_dict_from_arr(np_arr):
    ret_dict = {}
    for entry in np_arr:
        ret_dict[entry] = 1
    return ret_dict


CBC_LENGTH = 16
UMI_LENGTH = 12

exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
fig_pn = os.path.join(exp_pn, 'figures')
terra_pn = os.path.join(exp_pn, 'terra')
sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
fastq_dir1 = os.path.join(exp_pn, 'CROP_DO_fastqs_run1')
fastq_dir2 = os.path.join(exp_pn, 'CROP_DO_fastqs_run2')
barcode_txt_pn = os.path.join(terra_pn, 'Barcodes')

R2_fps_1 = glob.glob(os.path.join(fastq_dir1, '*_R2_*'))
R2_fps_2 = glob.glob(os.path.join(fastq_dir2, '*_R2_*'))
R2_fps = np.append(R2_fps_1, R2_fps_2)
R2_fps = [fp for fp in R2_fps if 'Undetermined' not in fp]

gw_guides_pn = os.path.join(sequencing_info_pn, 'gw_lib_guides.csv')
if not os.path.exists(gw_guides_pn):
	target_sgRNA_df = pd.read_csv(os.path.join(sequencing_info_pn, 'gw_lib_guides.txt'), sep='\t', header=None, names=['sgRNA_name', 'sgRNA_barcode'])
	target_sgRNA_barcodes = target_sgRNA_df['sgRNA_barcode'].to_numpy().astype(np.str)
	target_sgRNA_dict = make_dict_from_arr(target_sgRNA_barcodes)
	target_sgRNA_names = number_sgRNA_names(target_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str))
	pd.DataFrame({'sgRNA_name' : target_sgRNA_names, 'sgRNA_barcode' : target_sgRNA_barcodes}).to_csv(gw_guides_pn, index=False)
else:
	target_sgRNA_df = pd.read_csv(gw_guides_pn)
	target_sgRNA_barcodes = target_sgRNA_df['sgRNA_barcode'].to_numpy().astype(np.str)
	target_sgRNA_names =  target_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
	assert(target_sgRNA_barcodes.size == target_sgRNA_names.size)

control_sgRNA_df = pd.read_csv(os.path.join(sequencing_info_pn, 'human_TF_CTL_guides.txt'), sep='\t', header=None, names=['sgRNA_name', 'sgRNA_barcode'])
control_sgRNA_barcodes = control_sgRNA_df['sgRNA_barcode'].to_numpy().astype(np.str)
control_sgRNA_names =  control_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
assert(control_sgRNA_barcodes.size == control_sgRNA_names.size)

target_sgRNA_dict = make_dict_from_arr(target_sgRNA_barcodes)
control_sgRNA_dict = make_dict_from_arr(control_sgRNA_barcodes)

assert(len(control_sgRNA_barcodes[0]) == len(target_sgRNA_barcodes[0]))
sgRNA_barcode_length = len(target_sgRNA_barcodes[0])

plot_df = pd.DataFrame([])
for R2_fp in R2_fps:

	sample_prefix = R2_fp.split('/')[-1].split('_')[0]
	sample_suffix = R2_fp.split('/')[-1].split('S')[0].split('_')[1]
	channel = sample_prefix+'_'+sample_suffix
	print('Processing channel %s' % channel)

	if sample_prefix == 'CTL':
		sgRNA_dict = control_sgRNA_dict
	else:
		sgRNA_dict = target_sgRNA_dict

	tenx_barcodes = np.loadtxt(os.path.join(barcode_txt_pn, channel+'.txt'), dtype=str)
	tenx_dict = make_dict_from_arr(tenx_barcodes)

	open_R1_fp = gzip.open(R2_fp.replace('R2', 'R1'))
	open_R2_fp = gzip.open(R2_fp)
	R1_hits, R2_hits, loop_count = 0, 0, 0
	while True:

		R1_head = open_R1_fp.readline().decode().strip() # Header
		R1_seq = open_R1_fp.readline().decode().strip() # Sequence
		R1_line3 = open_R1_fp.readline().decode().strip()
		R1_quality = open_R1_fp.readline().decode().strip()

		R2_head = open_R2_fp.readline().decode().strip() # Header
		R2_seq = open_R2_fp.readline().decode().strip() # Sequence
		R2_line3 = open_R2_fp.readline().decode().strip()
		R2_quality = open_R2_fp.readline().decode().strip()

		tenx_CBC = R1_seq[:CBC_LENGTH]
		tenx_UMI = R1_seq[CBC_LENGTH:]

		R2_feature_barc = R2_seq[:sgRNA_barcode_length]

		if tenx_CBC in tenx_dict:
		    R1_hits += 1

		if R2_feature_barc in sgRNA_dict:
		    R2_hits += 1

		loop_count += 1

		if R1_head == '':
		    assert(R2_head == '')
		    break

	open_R1_fp.close()
	open_R2_fp.close()

	R1_map_rate = 100*(R1_hits / loop_count)
	R2_map_rate = 100*(R2_hits / loop_count)
	print('R1:', R1_map_rate, 'R2', R2_map_rate)

	plot_df = plot_df.append(pd.DataFrame({'Sample' : channel, 'Map Rate' : [R1_map_rate], 'Read' : 'R1'}), ignore_index=True)
	plot_df = plot_df.append(pd.DataFrame({'Sample' : channel, 'Map Rate' : [R2_map_rate], 'Read' : 'R2'}), ignore_index=True)

plot_df = plot_df.sort_values(by=['Sample'])

bp()

fig_save_pn = os.path.join(fig_pn, 'R1R2_mapping_rate.png')
fig, ax = plt.subplots()
sns.barplot(data=plot_df, x='Sample', y='Map Rate', hue='Read', ax=ax)
ax.set_ylabel('Map Rate (%)')
ax.set_title('CROP Dial-out Mapping Rate')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2])
plt.xticks(rotation=90)
fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
plt.close()
