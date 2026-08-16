import numpy as np
import pandas as pd
import os,glob
import scipy.sparse as sp_sparse
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
# Conda environment cs_scvelo

full_exp_pn = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/'
filtered_rna_pn = os.path.join(full_exp_pn, 'filtered_RNA')
sequencing_info_pn = os.path.join(full_exp_pn, 'sequencing_info')
CITE_expr_pn = os.path.join(full_exp_pn, 'CITE/CITE_expr')
design_mat_pn = os.path.join(full_exp_pn, 'CROP/design_mats')
fig_pn = os.path.join(full_exp_pn, 'figures')

gene_inds = np.load(os.path.join(filtered_rna_pn, 'filtered_gene_inds_100_cells.npy'))
gene_names = np.load(os.path.join(filtered_rna_pn, 'filtered_gene_names_100_cells.npy'))
assert(gene_inds.size == gene_names.size)
antibody_names = np.load(os.path.join(CITE_expr_pn, 'antibody_names.npy'))

# Load sgRNA information
gw_guides_pn = os.path.join(sequencing_info_pn, 'gw_lib_guides.csv')
target_sgRNA_df = pd.read_csv(gw_guides_pn)
target_sgRNA_names =  target_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)
control_sgRNA_df = pd.read_csv(os.path.join(sequencing_info_pn, 'human_TF_CTL_guides.txt'), sep='\t', header=None, names=['sgRNA_name', 'sgRNA_barcode'])
control_sgRNA_names =  control_sgRNA_df['sgRNA_name'].to_numpy().astype(np.str)

combined_sgRNA_names = np.append(target_sgRNA_names, control_sgRNA_names) # Targeting channels have both control and targeting sgRNAs

all_expr_RNA = np.array([])
all_expr_CITE = np.array([])
cond_categories = np.array([])
sample_categories = np.array([])
sgRNA_categories = np.array([])
MOI_num = np.array([])
cell_barcodes = np.array([])

all_designs = glob.glob(os.path.join(design_mat_pn, '*.sgRNA_design_mat.*.npz'))
for design in all_designs:

	curr_channel = design.split('/')[-1].split('.')[0]
	print('Processing %s' % curr_channel)

	if 'High' in curr_channel:
		cond = 'High'
		sgRNA_names = combined_sgRNA_names
	elif 'Low' in curr_channel:
		cond = 'Low'
		sgRNA_names = combined_sgRNA_names
	elif 'CTL' in curr_channel:
		cond = 'Control'
		sgRNA_names = control_sgRNA_names
	else:
		print('INVALID CONDITION')
		bp()

	# Load RNA
	curr_RNA_expr = np.load(os.path.join(filtered_rna_pn, curr_channel+'.expr_mat.npy'))[:, gene_inds]
	curr_RNA_barcs = np.load(os.path.join(filtered_rna_pn, curr_channel+'.barcodes.npy'))
	assert(curr_RNA_barcs.size == curr_RNA_expr.shape[0])
	assert(gene_names.size == curr_RNA_expr.shape[1])
	print('RNA successfully loaded')

	# Load CITE
	curr_CITE_expr = np.load(os.path.join(CITE_expr_pn, curr_channel+'.CITE_norm_expr.npy'))
	curr_CITE_barcs = np.load(os.path.join(CITE_expr_pn, curr_channel+'.CITE_barcodes.npy'))
	assert(curr_CITE_barcs.size == curr_CITE_expr.shape[1])
	print('CITE successfully loaded')

	# Load CROP
	curr_design_mat = np.array(sp_sparse.load_npz(design).todense())
	assert(curr_design_mat.shape[0] == curr_RNA_barcs.size)
	assert(curr_design_mat.shape[1] == sgRNA_names.size)
	print('CROP successfully loaded')

	# Only look at cells that passed QC for CITE
	good_CITE_inds = np.where(np.isin(curr_CITE_barcs, curr_RNA_barcs))[0]
	assert(np.all(np.char.equal(curr_CITE_barcs[good_CITE_inds], curr_RNA_barcs)))

	# Find MOI = 1 cells
	curr_MOI = curr_design_mat.sum(1)
	MOI1_inds = np.where(curr_MOI == 1)[0]
	MOI2_inds = np.where(curr_MOI == 2)[0]

	assert(np.intersect1d(MOI1_inds, MOI2_inds).size == 0)

	# Find sgRNAs
	MOI12_inds = np.union1d(MOI1_inds, MOI2_inds)
	assert(MOI12_inds.size == (MOI1_inds.size + MOI2_inds.size))

	curr_sgRNAs = np.full(curr_RNA_barcs.size, '', dtype=object)

	for MOI12_ind in MOI12_inds:
		sgRNA_ind = np.where(curr_design_mat[MOI12_ind, :] != 0)[0]
		assert(sgRNA_ind.size == 1 or sgRNA_ind.size == 2)

		if sgRNA_ind.size == 1:
			curr_sgRNAs[MOI12_ind] = sgRNA_names[sgRNA_ind[0]]
		elif sgRNA_ind.size == 2:
			curr_sgRNAs[MOI12_ind] = sgRNA_names[sgRNA_ind[0]] + '--' + sgRNA_names[sgRNA_ind[1]]

	curr_sgRNAs = curr_sgRNAs.astype(np.str)

	curr_annotated_barcodes = [curr_channel + '.' + barc for barc in curr_RNA_barcs]

	if all_expr_RNA.size == 0:
		all_expr_RNA = curr_RNA_expr
		all_expr_CITE = curr_CITE_expr.T[good_CITE_inds, :]
		cond_categories = np.full(curr_RNA_barcs.size, cond)
		sample_categories = np.full(curr_RNA_barcs.size, curr_channel)
		sgRNA_categories = curr_sgRNAs
		MOI_num = curr_MOI
		cell_barcodes = np.array(curr_annotated_barcodes).astype(np.str)
	else:
		all_expr_RNA = np.vstack((all_expr_RNA, curr_RNA_expr))
		all_expr_CITE = np.vstack((all_expr_CITE, curr_CITE_expr.T[good_CITE_inds, :]))
		cond_categories = np.append(cond_categories, np.full(curr_RNA_barcs.size, cond))
		sample_categories = np.append(sample_categories, np.full(curr_RNA_barcs.size, curr_channel))
		sgRNA_categories = np.append(sgRNA_categories, curr_sgRNAs)
		MOI_num = np.append(MOI_num, curr_MOI)
		cell_barcodes = np.append(cell_barcodes, np.array(curr_annotated_barcodes).astype(np.str))

assert(all_expr_RNA.shape[0] == all_expr_CITE.shape[0])
assert(all_expr_RNA.shape[0] == cell_barcodes.size)
assert(all_expr_RNA.shape[0] == cond_categories.size)
assert(all_expr_RNA.shape[0] == sample_categories.size)
assert(all_expr_RNA.shape[0] == sgRNA_categories.size)
assert(all_expr_RNA.shape[0] == MOI_num.size)
assert(all_expr_RNA.shape[1] == gene_names.size)
assert(all_expr_CITE.shape[1] == antibody_names.size)

bp()

np.save(os.path.join(full_exp_pn, 'adata_arrs/adata_parts/sgRNA_names.npy'), sgRNA_categories)

# UMI_save_pn = os.path.join(full_exp_pn, 'adata_arrs/UMI_count.npy')
# RNA_CITE_hstack = np.hstack((all_expr_RNA, all_expr_CITE))
# np.save(UMI_save_pn, RNA_CITE_hstack.sum(1)) # Technically CITE has already been normalized, but should have negligible effect

# Setup anndata structure
adata = anndata.AnnData(X=all_expr_RNA)
adata.var = pd.DataFrame(index=gene_names)
adata.obs = pd.DataFrame(index=cell_barcodes)

# Preprocess rna abundance data
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

# Append CITE expression data to log normalized rna expression matrix
adata_RNA_CITE = anndata.AnnData(X=np.hstack((adata.X, all_expr_CITE)))
del adata
adata_RNA_CITE.var = pd.DataFrame(index=np.append(gene_names, np.char.add('CITE-', antibody_names)))
adata_RNA_CITE.obs = pd.DataFrame(index=cell_barcodes)
adata_RNA_CITE.obs['Condition'] = pd.Categorical(cond_categories)
adata_RNA_CITE.obs['Sample'] = pd.Categorical(sample_categories)
adata_RNA_CITE.obs['MOI'] = pd.Categorical(MOI_num)
adata_RNA_CITE.obs['sgRNA'] = pd.Categorical(sgRNA_categories)

# Calculate metrics / embeddings
sc.pp.pca(adata_RNA_CITE)
sc.pp.neighbors(adata_RNA_CITE)
# sc.tl.louvain(adata_RNA_CITE)
sc.tl.umap(adata_RNA_CITE)

# Save structure
adata_RNA_CITE.write(os.path.join(full_exp_pn, 'adata_RNA_CITE.h5ad'))
