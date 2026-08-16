"""
Stage 2 of the Resis / editorial "incomplete submission" package.

Writes every number behind every plotted panel to its own Excel file, named after the
panel it belongs to (Fig1c-1e_..., Fig2e_..., ExtFig1a-i_...), into 03_Source_Data/.
Resis does not accept Word tables, and the panel has to be identifiable from the
file name alone.

Panels that share one underlying table (e.g. Fig 1c/1d/1e are three views of the same
per-cell table) get one file whose name lists all of them, rather than three copies of
353,091 rows.

Panel inventory follows the captions in doc/manuscript.pdf: Fig 1 a-g, Fig 2 a-e,
Extended Data Fig 1 (CITE marker violins) and Extended Data Fig 2 (guide capture).
Fig 1a,1b and Fig 2a are schematics and have no numerical content.

Run:  conda run -n perturbme python src/manuscript_figures/build_source_data_xlsx.py
"""
import os
import shutil
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, pearsonr

warnings.filterwarnings('ignore')

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
BASE = '/home/wangh256/hanchen/Pert_PG/perturb-me/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp'
LM = os.path.join(BASE, 'CROP/linear_model/17_cells_per_target/all_features')
MGK = os.path.join(BASE, 'guide_seq/mageck_out')
CROP = os.path.join(BASE, 'CROP')
NF = os.path.join(REPO, 'figures/nature_figures')
OUT = os.path.join(REPO, 'submission/20260721_Resis_package/03_Source_Data')
LOGF = os.path.join(REPO, 'submission/20260721_Resis_package/.stage2.log')

COND_ORDER = ['Low', 'Control', 'High']
log, index = [], []

if os.path.isdir(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT)


def note(msg):
    log.append(msg)
    open(LOGF, 'w').write('\n'.join(log))


def emit(panel, name, df, contents, index_col=False, supp=''):
    """One file per panel. `name` is the file stem, `panel` the figure location."""
    # NB: never enable xlsxwriter's constant_memory here. pandas emits cells
    # column-by-column and constant_memory silently discards any write to a row
    # above the current one, dropping all but the last value of every column.
    path = os.path.join(OUT, name + '.xlsx')
    with pd.ExcelWriter(path, engine='xlsxwriter',
                        engine_kwargs={'options': {'nan_inf_to_errors': True}}) as xw:
        df.to_excel(xw, sheet_name=name[:31], index=index_col)
    index.append({'panel': panel, 'supplementary_item': supp, 'file': name + '.xlsx',
                  'rows': df.shape[0], 'columns': df.shape[1] + int(index_col),
                  'contents': contents})
    note(f'  {panel:16s} {name + ".xlsx":52s} {df.shape[0]:>8,} x {df.shape[1] + int(index_col):>5,}')


# ---------------------------------------------------------------- shared inputs
note('loading adata_RNA_CITE.h5ad ...')
adata = sc.read_h5ad(os.path.join(BASE, 'adata_RNA_CITE.h5ad'))
condition = adata.obs['Condition'].astype(str).to_numpy()
sample = adata.obs['Sample'].astype(str).to_numpy()
umap = adata.obsm['X_umap']
note(f'  adata {adata.shape[0]:,} cells x {adata.shape[1]:,} features\n')


def expr(feature):
    x = adata[:, feature].X
    return np.asarray(x.todense()).ravel() if hasattr(x, 'todense') else np.asarray(x).ravel()


def dist_stats(values, groups, group_order, value_name, extra=None):
    rows = []
    for g in group_order:
        v = values[groups == g]
        rows.append({
            **(extra or {}),
            'group': g, 'n_cells': len(v),
            'mean': np.mean(v), 'sd': np.std(v, ddof=1), 'min': np.min(v),
            'q1': np.percentile(v, 25), 'median': np.median(v),
            'q3': np.percentile(v, 75), 'max': np.max(v),
            'measure': value_name,
        })
    return pd.DataFrame(rows)


# ================================================================== FIGURE 1
hla_prot = expr('CITE-HLA_A')
hla_rna = expr('HLA-A')
emit('Fig 1c,1d,1e', 'Fig1c-1e_per_cell_UMAP_protein_mRNA', pd.DataFrame({
    'cell_barcode': adata.obs.index.to_numpy().astype(str),
    'Sample': sample,
    'Condition': condition,
    'UMAP_1': umap[:, 0].round(6),
    'UMAP_2': umap[:, 1].round(6),
    'CITE_HLA_ABC_protein_log': hla_prot.round(6),
    'HLA_A_mRNA_log': hla_rna.round(6),
}), 'One row per cell: UMAP coordinates (1c, 1d), HLA-A,B,C protein (1d, 1e) and HLA-A mRNA (1e)')

emit('Fig 1e', 'Fig1e_violin_statistics', pd.concat([
    dist_stats(hla_prot, condition, COND_ORDER, 'CITE-seq HLA-A,B,C protein (log)'),
    dist_stats(hla_rna, condition, COND_ORDER, 'HLA-A mRNA (log)'),
], ignore_index=True), 'n, mean, sd, min, Q1, median, Q3, max per condition for both violins')

rank_rows, gene_rank_rows = [], []
for gate, fn in [('HLA-Low vs Input', 'Low_vs_Input.sgrna_summary.txt'),
                 ('HLA-High vs Input', 'High_vs_Input.sgrna_summary.txt')]:
    d = pd.read_csv(os.path.join(MGK, fn), sep='\t')
    d = d[~d['Gene'].str.contains('NO_SITE|NON-GENE', na=False)].copy()
    med = d.groupby('Gene')['LFC'].median().sort_values()
    order = pd.concat([med.head(5), med.tail(15)])          # bottom 5 + top 15, as plotted
    gene_rank_rows.append(pd.DataFrame({'panel': gate, 'Gene': order.index,
                                        'median_LFC': order.to_numpy().round(6),
                                        'row_position_bottom_to_top': range(1, len(order) + 1)}))
    sub = d[d['Gene'].isin(order.index)][['Gene', 'sgrna', 'LFC', 'control_count', 'treatment_count',
                                          'control_mean', 'treat_mean', 'p.low', 'p.high', 'p.twosided', 'FDR']]
    sub.insert(0, 'panel', gate)
    rank_rows.append(sub.sort_values(['Gene', 'sgrna']))
emit('Fig 1f', 'Fig1f_plotted_sgRNA_log2FC', pd.concat(rank_rows, ignore_index=True),
     'Every sgRNA tick drawn in the two rank views, with its MAGeCK statistics')
emit('Fig 1f', 'Fig1f_gene_ranking', pd.concat(gene_rank_rows, ignore_index=True),
     'Median LFC per gene, used to select the 20 genes shown and their row order')

all_genes = np.load(os.path.join(MGK, 'all_genes.npy'), allow_pickle=True).astype(str)
control_mask = np.isin(all_genes, ['NO_SITE', 'ONE_NON-GENE_SITE'])
g1_rows, g1_stats = [], []
for gate, thr in [('HLA-Low', 40), ('HLA-High', 25)]:
    tag = gate.split('-')[1].lower()
    reads = np.load(os.path.join(MGK, f'reads_per_gene_{tag}.npy'))
    cells = np.load(os.path.join(MGK, f'cells_per_gene_{tag}.npy'))
    keep = ~control_mask
    g1_rows.append(pd.DataFrame({
        'panel': gate, 'gene_target': all_genes[keep],
        'bulk_gDNA_reads_per_target': reads[keep],
        'single_cells_recovered_per_target': cells[keep],
        'labelled_in_panel': cells[keep] > thr,
    }))
    rho, p_rho = spearmanr(reads[keep], cells[keep])
    r, p_r = pearsonr(reads[keep], cells[keep])
    g1_stats.append({'panel': gate, 'n_gene_targets': int(keep.sum()),
                     'spearman_rho': round(rho, 6), 'spearman_p': p_rho,
                     'pearson_r': round(r, 6), 'pearson_p': p_r,
                     'label_threshold_cells': thr,
                     'note': 'non-targeting controls (NO_SITE, ONE_NON-GENE_SITE) excluded'})
emit('Fig 1g', 'Fig1g_reads_vs_cells_per_target', pd.concat(g1_rows, ignore_index=True),
     'Every plotted point: bulk gDNA reads vs single cells recovered, per gene-target per gate')
emit('Fig 1g', 'Fig1g_correlation_statistics', pd.DataFrame(g1_stats),
     "Spearman rho annotated on each panel, with Pearson r and n for reference")

# ================================================================== FIGURE 2
cov_names = np.load(os.path.join(LM, 'cov_names.npy'), allow_pickle=True).astype(str)
feat_names = np.load(os.path.join(LM, 'feature_names.npy'), allow_pickle=True).astype(str)
B_mat = np.load(os.path.join(LM, 'EN_B_EM.npy'))


def remove_sparse(mat, rows, cols, approx_zero=0.05, row_thresh=0.25, col_thresh=0.8):
    row_sp = np.mean(np.abs(mat) <= approx_zero, axis=1)
    col_sp = np.mean(np.abs(mat) <= approx_zero, axis=0)
    return mat[row_sp < row_thresh][:, col_sp < col_thresh], rows[row_sp < row_thresh], cols[col_sp < col_thresh]


filt_B, filt_features, filt_covs = remove_sparse(B_mat, feat_names, cov_names)
note(f'\n  beta matrix {B_mat.shape} -> plotted submatrix {filt_B.shape}')

tm = pd.read_csv(os.path.join(NF, 'Fig2/Fig2B_target_modules.csv'))
gp = pd.read_csv(os.path.join(NF, 'Fig2/Fig2B_gene_programs.csv'))
col_ord = np.argsort([{t: i for i, t in enumerate(tm['target'])}[c] for c in filt_covs])
row_ord = np.argsort([{g: i for i, g in enumerate(gp['gene'])}[g] for g in filt_features])
plot_B = filt_B[row_ord, :][:, col_ord]
plot_features, plot_covs = filt_features[row_ord], filt_covs[col_ord]
corr_t = pd.DataFrame(filt_B[:, col_ord]).corr().fillna(0).to_numpy()
corr_f = pd.DataFrame(filt_B[row_ord, :].T).corr().fillna(0).to_numpy()

emit('Fig 2b', 'Fig2b_beta_matrix',
     pd.DataFrame(plot_B.round(6), index=pd.Index(plot_features, name='gene'), columns=plot_covs),
     'Heatmap values: ElasticNet+EM beta, genes x perturbation targets, in plotted order', index_col=True)
emit('Fig 2b', 'Fig2b_feature_correlation',
     pd.DataFrame(corr_f.round(4), index=pd.Index(plot_features, name='gene'), columns=plot_features),
     'Gene x gene Pearson r across beta profiles, in plotted order', index_col=True)
emit('Fig 2b', 'Fig2b_target_correlation',
     pd.DataFrame(corr_t.round(4), index=pd.Index(plot_covs, name='target'), columns=plot_covs),
     'Target x target Pearson r across beta profiles, in plotted order', index_col=True)

mg_genes = np.load(os.path.join(MGK, 'all_genes.npy'), allow_pickle=True).astype(str)
mg_hi = np.load(os.path.join(MGK, 'cells_per_gene_high.npy'))
mg_lo = np.load(os.path.join(MGK, 'cells_per_gene_low.npy'))
midx = {g: i for i, g in enumerate(mg_genes)}
emit('Fig 2b', 'Fig2b_cells_per_target_sidebar', pd.DataFrame({
    'plot_order': np.arange(1, len(plot_covs) + 1),
    'target': plot_covs,
    'cells_HLA_high': [mg_hi[midx[t]] if t in midx else np.nan for t in plot_covs],
    'cells_HLA_low': [mg_lo[midx[t]] if t in midx else np.nan for t in plot_covs],
    'target_module': [int(tm.loc[tm['target'] == t, 'target_module'].iloc[0]) for t in plot_covs],
}), 'High/Low enrichment strip: cells recovered per target, in plotted order')
emit('Fig 2b', 'Fig2b_target_modules', tm, '7 target modules; k-means on the target correlation matrix',
     supp='Supplementary Tables 5-9')
emit('Fig 2b', 'Fig2b_gene_programs', gp, '9 gene programs; k-means on the gene correlation matrix',
     supp='Supplementary Tables 5-9')

for panel, stem, fn, contents, supp in [
    ('Fig 2c', 'Fig2c_target_module_top8', 'Fig2C_target_module_top8.csv',
     'The 8 targets printed in each module box; ranked by minimum MAGeCK pos-rank in either gate vs Input',
     'Supplementary Tables 5-9'),
    ('Fig 2c', 'Fig2c_gene_program_top8', 'Fig2C_gene_program_top8.csv',
     'The 8 genes printed in each program box; ranked by L2 norm of the beta row',
     'Supplementary Tables 5-9'),
    ('Fig 2c', 'Fig2c_module_program_edges', 'Fig2C_module_program_edges.csv',
     'Edges: module x program mean beta; drawn when |mean beta| > 0.15', 'Supplementary Tables 5-9'),
    ('Fig 2d', 'Fig2d_schematic_gene_lists', 'Fig2d_gene_reconciliation.csv',
     'Curated gene list printed on each compartment of the MHC-I schematic', ''),
]:
    emit(panel, stem, pd.read_csv(os.path.join(NF, 'Fig2', fn)), contents, supp=supp)

SANKEY = [
    ('TC5 JAK/STAT Core', 'FC7 HLA Program', 0.45, 'blue: KO -> HLA low'),
    ('TC5 JAK/STAT Core', 'FC5 Proliferation', 0.12, 'blue: KO -> HLA low'),
    ('TC3 Biosynthesis', 'FC2 ISR Stress', 0.19, 'blue: KO -> HLA low'),
    ('TC3 Biosynthesis', 'FC7 HLA Program', 0.29, 'blue: KO -> HLA low'),
    ('TC1 Proteostasis/UPS', 'FC5 Proliferation', 0.21, 'blue: KO -> HLA low'),
    ('TC2 Epigenetic Hub', 'FC8 Chromatin', 0.28, 'red: KO -> HLA high'),
    ('TC2 Epigenetic Hub', 'FC7 HLA Program', 0.31, 'red: KO -> HLA high'),
    ('TC4 QC/Folding', 'FC1 ER Stress', 0.24, 'red: KO -> HLA high'),
    ('TC4 QC/Folding', 'FC7 HLA Program', 0.19, 'red: KO -> HLA high'),
    ('TC6 Golgi Machinery', 'FC4 Trafficking', 0.26, 'red: KO -> HLA high'),
    ('TC6 Golgi Machinery', 'FC7 HLA Program', 0.29, 'red: KO -> HLA high'),
    ('TC6 Golgi Machinery', 'FC1 ER Stress', 0.11, 'red: KO -> HLA high'),
    ('TC0 Retromer/Recycling', 'FC6 Cytoskeletal', 0.22, 'red: KO -> HLA high'),
    ('TC0 Retromer/Recycling', 'FC7 HLA Program', 0.18, 'red: KO -> HLA high'),
    ('FC7 HLA Program', 'Phenotype: surface HLA', 0.50, 'blue: KO -> HLA low'),
    ('FC2 ISR Stress', 'Phenotype: surface HLA', 0.22, 'blue: KO -> HLA low'),
    ('FC5 Proliferation', 'Phenotype: surface HLA', 0.15, 'blue: KO -> HLA low'),
    ('FC8 Chromatin', 'Phenotype: surface HLA', 0.30, 'red: KO -> HLA high'),
    ('FC1 ER Stress', 'Phenotype: surface HLA', 0.25, 'red: KO -> HLA high'),
    ('FC4 Trafficking', 'Phenotype: surface HLA', 0.35, 'red: KO -> HLA high'),
    ('FC6 Cytoskeletal', 'Phenotype: surface HLA', 0.20, 'red: KO -> HLA high'),
]
emit('Fig 2e', 'Fig2e_sankey_links',
     pd.DataFrame(SANKEY, columns=['source_node', 'target_node', 'link_width', 'axis_colour']),
     'Sankey node pairs, link widths and axis colour')

SANKEY_NODES = [
    ('target module', 'TC0 Retromer/Recycling', ''),
    ('target module', 'TC1 Proteostasis/UPS', ''),
    ('target module', 'TC2 Epigenetic Hub', ''),
    ('target module', 'TC3 Biosynthesis Hub', ''),
    ('target module', 'TC4 QC/Folding', ''),
    ('target module', 'TC5 JAK/STAT Core', ''),
    ('target module', 'TC6 Golgi Machinery', ''),
    ('gene program', 'FC1 ER Stress', 'MANF, HSPA5, CALU, CANX'),
    ('gene program', 'FC2 ISR Stress', 'ATF4, ASNS, DDIT3/4'),
    ('gene program', 'FC4 Trafficking', 'GOLGA2, RAB1A, KDELR2, SEC61A1'),
    ('gene program', 'FC5 Proliferation', 'MKI67, CDC5L, CCNB1, TOP2A'),
    ('gene program', 'FC6 Cytoskeletal', 'ACTB, TUBA1B, TMSB4X, CFL1'),
    ('gene program', 'FC7 HLA Program', 'BST2, GBP1, HLA-A/B/C'),
    ('gene program', 'FC8 Chromatin', 'HMGB2, BRD2, BRD4, DAXX'),
    ('phenotype', 'Phenotype: surface HLA', ''),
]
emit('Fig 2e', 'Fig2e_sankey_nodes',
     pd.DataFrame(SANKEY_NODES, columns=['column', 'node_label', 'genes_printed_on_node']),
     'Sankey node labels and the representative genes printed inside each node')

emit('Fig 2b (source)', 'Fig2b_beta_matrix_unfiltered',
     pd.DataFrame(B_mat.round(6), index=pd.Index(feat_names, name='gene'), columns=cov_names),
     'Full 18,683 x 231 beta matrix before the sparsity filter, so the filter can be reproduced',
     index_col=True)

# ========================================================= EXTENDED DATA
# ED Fig 1 is now only the CITE-seq marker violin row; the earlier ED Fig 1a i-iii panels
# (cells per channel, per-sample HLA ridgeline, IFNg signature) are no longer in the paper.
ED1_MARKERS = [('CITE-HLA_D', 'HLA-DR'), ('CITE-CD274', 'CD274_PD-L1'), ('CITE-CD47', 'CD47'),
               ('CITE-CD49f', 'CD49f'), ('CITE-CD58', 'CD58')]
marker_vals = {label: expr(feat) for feat, label in ED1_MARKERS}

ed_cell = pd.DataFrame({
    'cell_barcode': adata.obs.index.to_numpy().astype(str),
    'Sample': sample,
    'Condition': condition,
})
for label in marker_vals:
    ed_cell['CITE_' + label] = marker_vals[label].round(6)
emit('Ext Fig 1', 'ExtFig1_per_cell_CITE_markers', ed_cell,
     'Violin input, one row per cell: HLA-DR, CD274 (PD-L1), CD47, CD49f and CD58')

emit('Ext Fig 1', 'ExtFig1_violin_statistics',
     pd.concat([dist_stats(v, condition, COND_ORDER, label) for label, v in marker_vals.items()],
               ignore_index=True),
     'n, mean, sd, min, Q1, median, Q3, max per condition for each of the five violins')

gc = pd.read_csv(os.path.join(NF, 'FigS_guide_capture/guide_capture_per_target.csv'))
emit('Ext Fig 2', 'ExtFig2_guide_capture_per_target', gc,
     'All three panels: bulk gDNA reads, single cells, guide UMIs and guide reads per gene-target per gate')

# Same filter as the plotting script (FigS_guide_capture_vs_bulk.py): drop controls,
# keep targets with >= 5 recovered cells in the gate. rho is reported both over all such
# targets and over the 221 targets that enter the linear model.
MINC = 5
model_targets = set(pd.read_csv(os.path.join(REPO, 'results/figure2_moi1/target_clusters_moi1.csv'))['target']) \
    - {'HIGH_CTL', 'LOW_CTL', 'NO_SITE', 'ONE_NON-GENE_SITE', 'G1'}
gc_stats = []
for gate in ['high', 'low']:
    sub = gc[(~gc['is_control']) & (gc[f'cells_{gate}'] >= MINC)]
    mdl = sub[sub['gene'].isin(model_targets)]
    upc, bpc, cells = sub[f'umis_per_cell_{gate}'], sub[f'bulkreads_per_cell_{gate}'], sub[f'cells_{gate}']
    gc_stats.append({
        'gate': f'HLA-{gate.capitalize()}',
        'min_cells_per_target': MINC,
        'n_gene_targets': len(sub), 'n_model_targets': len(mdl),
        'median_guide_umis_per_cell': round(float(upc.median()), 3),
        'sd_guide_umis_per_cell': round(float(upc.std()), 3),
        'cv_guide_umis_per_cell': round(float(upc.std() / upc.mean()), 3),
        'cv_bulk_reads_per_cell': round(float(bpc.std() / bpc.mean()), 3),
        'spearman_rho_all_targets': round(float(spearmanr(cells, upc, nan_policy='omit')[0]), 3),
        'spearman_rho_model_targets': round(float(spearmanr(mdl[f'cells_{gate}'],
                                                           mdl[f'umis_per_cell_{gate}'],
                                                           nan_policy='omit')[0]), 3),
    })
emit('Ext Fig 2', 'ExtFig2_guide_capture_statistics', pd.DataFrame(gc_stats),
     'The numbers quoted in the caption: median and s.d. of guide UMIs per cell, its Spearman '
     'correlation with recovered cells, and its CV against bulk reads per cell')

FG = os.path.join(NF, 'FigS_FACS_gating')
if os.path.exists(os.path.join(FG, 'ExtFig_FACS_gating_statistics.csv')):
    emit('Ext Fig 3', 'ExtFig3_FACS_gating_statistics',
         pd.read_csv(os.path.join(FG, 'ExtFig_FACS_gating_statistics.csv')),
         'Every gate boundary and the events retained at each step, per sample and per sorter')
    emit('Ext Fig 3', 'ExtFig3_FACS_gating_events',
         pd.read_csv(os.path.join(FG, 'ExtFig_FACS_gating_events.csv')),
         'One row per recorded FACS event: FSC/SSC, mKate2, HLA-A,B,C and gate memberships')
else:
    note('  MISSING ExtFig3_FACS_gating: run src/manuscript_figures/ExtFig_FACS_gating.py first')

# Separate ED figure - cells per guide / per gene coverage distributions (Aviv's request).
# Precomputed on the MOI=1 singlet data by ExtFig_cells_per_guide_gene.py so the source
# numbers match the plotted histograms exactly (the ~4 singlet-cells/guide baseline is
# per-guide, MOI=1).
CPG = os.path.join(NF, 'FigS_cells_per_guide_gene')
if os.path.exists(os.path.join(CPG, 'cells_per_guide.csv')):
    emit('Ext Fig (coverage)', 'ExtFig_cells_per_guide',
         pd.read_csv(os.path.join(CPG, 'cells_per_guide.csv')),
         'Coverage-distribution histogram: sgRNA-singlet cells (MOI=1) per guide, one row per guide')
    emit('Ext Fig (coverage)', 'ExtFig_cells_per_gene',
         pd.read_csv(os.path.join(CPG, 'cells_per_gene.csv')),
         'Coverage-distribution histogram: sgRNA-singlet cells (MOI=1) per gene-target, one row per gene')
else:
    note('  MISSING ExtFig_cells_per_guide_gene: run src/manuscript_figures/ExtFig_cells_per_guide_gene.py first')

tn = np.load(os.path.join(CROP, 'design_mats/target_names.npy'), allow_pickle=True).astype(str)
cth = np.load(os.path.join(CROP, 'cells_per_target_high.npy')).astype(int)
ctl = np.load(os.path.join(CROP, 'cells_per_target_low.npy')).astype(int)
emit('Methods', 'Supporting_target_retention_per_gate', pd.DataFrame({
    'target': tn, 'cells_HLA_high': cth, 'cells_HLA_low': ctl, 'cells_total': cth + ctl,
    'passes_filter_gt17_in_either_gate': (cth > 17) | (ctl > 17),
}), 'Cells per gene-target per gate and the >17-in-either-gate retention filter')

# Supplementary Tables 1-9 are not re-emitted here: they ship verbatim in
# doc/Supplementary Files (copied into the package as 04_Supplementary). Tables 5-9 are
# byte-identical to the Fig 2b/2c tables above; Table 4 is byte-identical to the MAGeCK
# output behind Fig 1f. They are listed in the index as pointers.
for supp, where, contents in [
    ('Supplementary Table 1', 'TableS1_CITE_antibody_panel.csv/.xlsx', 'CITE-seq antibody panel'),
    ('Supplementary Table 2', 'Table S2 CRISPR Library.xlsx', 'Guide sequences for both libraries'),
    ('Supplementary Table 3', 'TableS3_bulk_vs_perturbme_concordance.csv', 'Bulk screen vs Perturb-ME concordance'),
    ('Supplementary Table 4', 'Supplementary Table 4 - MAGeCK bulk screen/', 'Full MAGeCK output; same files as Fig 1f'),
    ('Supplementary Table 5', 'Supplementary Table 5 - Fig2B target modules.csv', 'identical to Fig2b_target_modules.xlsx'),
    ('Supplementary Table 6', 'Supplementary Table 6 - Fig2B gene programs.csv', 'identical to Fig2b_gene_programs.xlsx'),
    ('Supplementary Table 7', 'Supplementary Table 7 - Fig2C target modules.csv', 'identical to Fig2c_target_module_top8.xlsx'),
    ('Supplementary Table 8', 'Supplementary Table 8 - Fig2C module-program edges.csv', 'identical to Fig2c_module_program_edges.xlsx'),
    ('Supplementary Table 9', 'Supplementary Table 9 - Fig2C gene programs.csv', 'identical to Fig2c_gene_program_top8.xlsx'),
]:
    index.append({'panel': 'Supp table', 'supplementary_item': supp,
                  'file': '04_Supplementary/' + where, 'rows': '', 'columns': '',
                  'contents': contents})

# ---------------------------------------------------------------- index
idx = pd.DataFrame(index)
idx.to_csv(os.path.join(OUT, '00_INDEX.csv'), index=False)
with pd.ExcelWriter(os.path.join(OUT, '00_INDEX.xlsx'), engine='xlsxwriter') as xw:
    idx.to_excel(xw, sheet_name='INDEX', index=False)

note('\ndone')
for f in sorted(os.listdir(OUT)):
    note(f'  {os.path.getsize(os.path.join(OUT, f))/1e6:8.1f} MB  {f}')
print('\n'.join(log))
