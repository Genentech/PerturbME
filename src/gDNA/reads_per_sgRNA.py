import os
from pdb import set_trace as bp
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import scipy.stats as sp_stats
from sklearn.decomposition import PCA

exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/20200809_exp/'
figures_pn = os.path.join(exp_pn, 'figures')
sequencing_info_pn = os.path.join(exp_pn, 'sequencing_info')
gDNA_pn = os.path.join(exp_pn, 'guide_seq')
mageck_out_pn = os.path.join(gDNA_pn, 'mageck_out')

# all_samples_count = pd.read_csv(os.path.join(mageck_out_pn, 'all_samples.count.txt'), sep='\t')
all_samples_count = pd.read_csv(os.path.join(mageck_out_pn, 'Input_samples.count.txt'), sep='\t')
reads_per_sgRNA = all_samples_count.iloc[:, 2:].to_numpy().sum(1)

fig, ax = plt.subplots()
ax.hist(reads_per_sgRNA, bins=300)
ax.set_xlabel('Reads / sgRNA'), ax.set_ylabel('Count')
ax.set_xlim([0, 30e3])
fig.savefig(os.path.join(figures_pn, 'mageck_reads_per_sgRNA_Input'), bbox_inches='tight', dpi=800)
plt.close()
