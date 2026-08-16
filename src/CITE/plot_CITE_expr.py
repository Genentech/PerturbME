import os,glob
from pdb import set_trace as bp
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
barc_txt_path = os.path.join(exp_pn, 'terra/Barcodes')
CITE_path = os.path.join(exp_pn, 'CITE')
fig_pn = os.path.join(exp_pn, 'figures/CITE')
load_pn = os.path.join(CITE_path, 'CITE_expr')

antibody_names = np.load(os.path.join(load_pn, 'antibody_names.npy'))
all_txts = glob.glob(os.path.join(barc_txt_path, '*.txt'))

plot_df = pd.DataFrame()
for txt in all_txts:

	channel = txt.split('/')[-1].split('.')[0]
	print('Processing %s' % channel)
	curr_CITE_expr = np.load(os.path.join(load_pn, '%s.CITE_norm_expr.npy' % channel))

	if 'High' in channel:
		cond = 'High'
	elif 'Low' in channel:
		cond = 'Low'
	elif 'CTL' in channel:
		cond = 'Control'
	else:
		print('INVALID CONDITION')
		bp()

	curr_df = pd.DataFrame()
	for antibody_i, antibody in enumerate(antibody_names):
		curr_df = curr_df.append(pd.DataFrame({'Antibody' : antibody, 'Condition' : cond, 'Expression' : curr_CITE_expr[antibody_i, :]}), ignore_index=True)

	plot_df = plot_df.append(curr_df, ignore_index=True)

fig_save_pn = os.path.join(fig_pn, 'CITE_expr_3conditions.png')
fig, ax = plt.subplots(figsize = [3*6.4, 3*4.8])
sns.violinplot(data=plot_df, x='Expression', y='Antibody', hue='Condition', ax=ax)
ax.set_xlabel('Log(Normalized Expression)')
fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
plt.close()

# Make separate plot for each antibody
for antibody in antibody_names:
	fig_save_pn = os.path.join(fig_pn, '%s.png' % antibody)

	if not os.path.exists(fig_save_pn):
		curr_plot_df = plot_df.loc[plot_df['Antibody'] == antibody]

		fig, ax = plt.subplots()
		sns.violinplot(data=curr_plot_df, x='Expression', y='Antibody', hue='Condition', ax=ax)
		ax.set_xlabel('Log(Normalized Expression)')
		ax.set_ylabel('')
		ax.set_xlim([-1, 12])
		fig.savefig(fig_save_pn, bbox_inches='tight', dpi=800)
		plt.close()

bp()
