import numpy as np
import pandas as pd
import seaborn as sns
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from adjustText import adjust_text

def main():

	exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	figures_pn = os.path.join(exp_pn, 'figures')
	sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
	CROP_dir = os.path.join(exp_pn, 'CROP')
	design_dir = os.path.join(CROP_dir, 'design_mats/')

	cells_per_target_high = np.load(os.path.join(CROP_dir, 'cells_per_target_high.npy'))
	cells_per_target_low = np.load(os.path.join(CROP_dir, 'cells_per_target_low.npy'))
	target_names = np.load(os.path.join(design_dir, 'target_names.npy'))

	# pd.DataFrame(data=np.column_stack((cells_per_target_high, cells_per_target_low)), index=target_names, columns=['HLA High', 'HLA Low']).to_csv(os.path.join(CROP_dir, 'cells_per_target.csv'))

	assert(cells_per_target_high.size == target_names.size)
	assert(cells_per_target_low.size == target_names.size)

	# Vary threshold of targets to include + see number of targets
	thresholds = np.arange(1, 50, 1)
	included_targets_high = np.zeros(thresholds.size)
	included_targets_low = np.zeros(thresholds.size)
	for threshold_i, threshold in enumerate(thresholds):
		included_targets_high[threshold_i] = np.where(cells_per_target_high > threshold)[0].size
		included_targets_low[threshold_i] = np.where(cells_per_target_low > threshold)[0].size

	fig, ax = plt.subplots()
	ax.plot(thresholds, included_targets_high, label='HLA High')
	ax.plot(thresholds, included_targets_low, label='HLA Low')
	ax.set_xlabel('Cells / Target Threshold'), ax.set_ylabel('# Targets')
	ax.legend()
	ax.set_ylim([0, 2000])
	fig.savefig(os.path.join(figures_pn, 'targets_to_include.png'), dpi=800, bbox_inches='tight')
	plt.close()

	# Scatter of cells / target high vs low
	fig, ax = plt.subplots()
	ax.scatter(cells_per_target_high, cells_per_target_low, s=1)
	ax.set_xlabel('High'), ax.set_ylabel('Low')
	ax.set_title('Cells / Target')
	label_inds = np.where( (cells_per_target_low > 40) | (cells_per_target_high > 19) )[0]
	text_objs = [ax.text(cells_per_target_high[ind], cells_per_target_low[ind], target_names[ind], fontsize=4, fontweight='bold') for ind in label_inds]
	adjust_text(text_objs, arrowprops=dict(arrowstyle='-', color='red', lw=0.5))
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_scatter.png'), dpi=800, bbox_inches='tight')
	plt.close()

	# Histrogram of cells / target distributions
	xlims = [0, 40]
	fig, ax = plt.subplots()
	ax.hist(cells_per_target_high, bins=40)
	ax.set_xlabel('Cells / Target'), ax.set_ylabel('Count')
	ax.set_title('HLA High')
	ax.set_xlim(xlims)
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_high_hist.png'), dpi=800, bbox_inches='tight')
	plt.close()
	fig, ax = plt.subplots()
	ax.hist(cells_per_target_low, bins=175)
	ax.set_xlabel('Cells / Target'), ax.set_ylabel('Count')
	ax.set_title('HLA Low')
	ax.set_xlim(xlims)
	fig.savefig(os.path.join(figures_pn, 'cells_per_target_low_hist.png'), dpi=800, bbox_inches='tight')
	plt.close()


if __name__ == "__main__":
    main()
