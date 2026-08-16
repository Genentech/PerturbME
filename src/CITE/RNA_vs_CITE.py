import numpy as np
import pandas as pd
import glob,os
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
import scanpy as sc

def main():

	full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
	fig_pn = os.path.join(full_exp_pn, 'figures/CITE/RNA_vs_CITE')

	adata_RNA_CITE = anndata.read(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))
	feature_names = adata_RNA_CITE.var.index.to_numpy().astype(np.str)
	cond_labels = adata_RNA_CITE.obs['Condition'].to_numpy().astype(np.str)

	antibody_to_gene_df = pd.read_csv(os.path.join(full_exp_pn, 'sequencing_info/antibody_to_control.csv'))
	antibodies = antibody_to_gene_df['Antibody'].to_numpy().astype(np.str)
	genes_per_antibody = antibody_to_gene_df['Gene name'].to_numpy().astype(np.str)

	assert(antibodies.size == genes_per_antibody.size)

	color_dict = {'Low' : '#e1812c', 'High' : '#3a923a', 'Control' : '#3274a1'}
	for antibody_i, antibody in enumerate(antibodies):

		CITE_oi = antibody
		RNA_oi = genes_per_antibody[antibody_i]

		if len(RNA_oi.split(',')) > 1: # For now just look at first transcript for antibodies with multiple genes (e.g. HLA-A,B,C and HLA-DR)
			RNA_oi = RNA_oi.split(',')[0]

		CITE_expr = adata_RNA_CITE.X[:, np.where(feature_names == 'CITE-'+CITE_oi)[0]].flatten()
		RNA_expr = adata_RNA_CITE.X[:, np.where(feature_names == RNA_oi)[0]].flatten()

		plot_df_CITE = pd.DataFrame()
		plot_df_CITE = plot_df_CITE.append(pd.DataFrame({'Assay' : 'CITE', 'Expression' : CITE_expr, 'Condition' : cond_labels}), ignore_index=True)
		plot_df_RNA = pd.DataFrame()
		plot_df_RNA = plot_df_RNA.append(pd.DataFrame({'Assay' : 'RNA', 'Expression' : RNA_expr, 'Condition' : cond_labels}), ignore_index=True)

		bw_arr = np.arange(0.1, 1.1, 0.1)

		for bw in bw_arr:
			fig_save_pn_RNA = os.path.join(fig_pn, '%s_RNA_vs_CITE_violin_bw%.2f_RNA.png' % (CITE_oi, bw))
			fig_save_pn_CITE = os.path.join(fig_pn, '%s_RNA_vs_CITE_violin_bw%.2f_CITE.png' % (CITE_oi, bw))
			if not os.path.exists(fig_save_pn_RNA):
				fig, ax = plt.subplots()
				sns.violinplot(data=plot_df_CITE, y='Expression', x='Condition', ax=ax, bw=bw, palette=color_dict)
				ax.set_xlabel(CITE_oi)
				ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
				fig.savefig(fig_save_pn_CITE, bbox_inches='tight', dpi=800)
				plt.close()

				fig, ax = plt.subplots()
				sns.violinplot(data=plot_df_RNA, y='Expression', x='Condition', ax=ax, bw=bw, palette=color_dict)
				ax.set_xlabel(CITE_oi)
				ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
				fig.savefig(fig_save_pn_RNA, bbox_inches='tight', dpi=800)
				plt.close()

if __name__ == "__main__":
    main()
