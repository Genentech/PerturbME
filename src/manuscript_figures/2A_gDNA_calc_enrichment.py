import os
from pdb import set_trace as bp
import numpy as np
import pandas as pd
import seaborn as sns
import glob
import matplotlib.pyplot as plt
from scipy.stats import norm
from adjustText import adjust_text

full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
figures_pn = os.path.join(full_exp_pn, 'figures/manuscript_figures/Fig2')

# Calculate number of cells for given coverage
target_coverage = 50
percent_impactful_sgRNAs = np.arange(1, 8.1, 0.01)
enrichment_of_impactful_sgRNAs = np.arange(1, 9.1, 0.01)
total_number_sgRNAs = 60e3

total_num_cells = np.zeros((percent_impactful_sgRNAs.size, enrichment_of_impactful_sgRNAs.size))
for perc_i, curr_percent in enumerate(percent_impactful_sgRNAs):
	for enrich_i, curr_enrichment in enumerate(enrichment_of_impactful_sgRNAs):
		num_impactful_sgRNAs = (curr_percent/100) * total_number_sgRNAs
		num_boring_sgRNAs = total_number_sgRNAs - num_impactful_sgRNAs
		cells_for_impactful_sgRNAs = num_impactful_sgRNAs*target_coverage
		num_cells_per_boring_sgRNA = target_coverage / curr_enrichment
		cells_for_boring_sgRNAs = num_boring_sgRNAs * num_cells_per_boring_sgRNA

		total_num_cells[perc_i, enrich_i] = cells_for_impactful_sgRNAs+cells_for_boring_sgRNAs


fig_save_pn = os.path.join(figures_pn, '2A_num_cells_%d_coverage.pdf' % target_coverage)
fig, ax = plt.subplots()
ax_ret = ax.imshow(total_num_cells, aspect=(total_num_cells.shape[1]/total_num_cells.shape[0]), cmap='tab20')
cbar = fig.colorbar(ax_ret, ax=ax, extend='both')
cbar.minorticks_on()
ax.set_ylim([total_num_cells.shape[0]-0.5, -0.5])
ax.set_ylabel('Percent Impactful sgRNAs')
ax.set_xlabel('Enrichment of Impactful sgRNAs')
ax.set_title('# Cells for %dx Coverage' % target_coverage)
ylabel_inds = np.where(np.around(percent_impactful_sgRNAs, decimals=2)==np.round(percent_impactful_sgRNAs))[0]
ax.set_yticks(ylabel_inds), ax.set_yticklabels(np.round(percent_impactful_sgRNAs[ylabel_inds]))
xlabel_inds = np.where(np.around(enrichment_of_impactful_sgRNAs, decimals=2)==np.round(enrichment_of_impactful_sgRNAs))[0]
ax.set_xticks(xlabel_inds), ax.set_xticklabels(np.round(enrichment_of_impactful_sgRNAs[xlabel_inds]))
fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
plt.close()