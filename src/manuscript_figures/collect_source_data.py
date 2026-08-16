"""
Collect the tables/CSVs behind every plotted panel of the manuscript figures into
one folder (figures/source_data/), organised by figure.

Two kinds of entry:
  COPY   - a table that already exists on disk and was read/written by the plotting
           notebook or script; copied verbatim.
  EXPORT - a panel whose plotted values live in .npy arrays; re-exported here as a
           labelled CSV (row/column names attached) so the numbers are readable.

Panels driven straight off adata_RNA_CITE.h5ad (per-cell UMAP / violin / histogram
panels) have no intermediate table; they are listed in the manifest as H5AD-only.

Run:  conda run -n perturbme python src/manuscript_figures/collect_source_data.py
"""
import os, shutil
import numpy as np, pandas as pd

REPO  = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
BASE  = '/home/wangh256/hanchen/Pert_PG/perturb-me/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp'
LM    = os.path.join(BASE, 'CROP/linear_model/17_cells_per_target/all_features')
MGK   = os.path.join(BASE, 'guide_seq/mageck_out')
CROP  = os.path.join(BASE, 'CROP')
PERM  = '/home/wangh256/hanchen/Pert_PG/perturb-me/src/CROP/linear_model/permutation_test/results'
OUT   = os.path.join(REPO, 'figures/source_data')

log = []


def copy(panel, src, dst_name, note):
    dst_dir = os.path.join(OUT, panel.split('|')[0])
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, dst_name)
    if not os.path.exists(src):
        log.append((panel, dst_name, 'MISSING', src, note)); return
    shutil.copy2(src, dst)
    log.append((panel, dst_name, 'COPY', src, note))


def export(panel, dst_name, df, src, note, index=True):
    dst_dir = os.path.join(OUT, panel.split('|')[0])
    os.makedirs(dst_dir, exist_ok=True)
    df.to_csv(os.path.join(dst_dir, dst_name), index=index)
    log.append((panel, dst_name, 'EXPORT', src, note))


# ----------------------------------------------------------------- Figure 1
copy('Fig1|1f', os.path.join(MGK, 'Low_vs_Input.sgrna_summary.txt'),
     'Fig1f_MAGeCK_Low_vs_Input.sgrna_summary.txt', 'per-sgRNA LFC, HLA-Low vs Input')
copy('Fig1|1f', os.path.join(MGK, 'High_vs_Input.sgrna_summary.txt'),
     'Fig1f_MAGeCK_High_vs_Input.sgrna_summary.txt', 'per-sgRNA LFC, HLA-High vs Input')
copy('Fig1|1g', os.path.join(REPO, 'figures/fig1g_recompute/reads_vs_cells_per_target.csv'),
     'Fig1g_reads_vs_cells_per_target.csv', 'plotted points: bulk gDNA reads vs cells recovered, per gene-target per gate')
copy('Fig1|1g', os.path.join(REPO, 'figures/fig1g_recompute/correlations.csv'),
     'Fig1g_correlations.csv', 'Spearman/Pearson values annotated on the panel')
copy('Fig1|1g', os.path.join(REPO, 'figures/fig1g_recompute/corr_vs_depth_threshold.csv'),
     'Fig1g_corr_vs_depth_decile.csv', 'decile inset: correlation vs bulk read-depth threshold')

# ----------------------------------------------------------------- Figure 2
tn  = np.load(os.path.join(CROP, 'design_mats/target_names.npy'), allow_pickle=True).astype(str)
cth = np.load(os.path.join(CROP, 'cells_per_target_high.npy'))
ctl = np.load(os.path.join(CROP, 'cells_per_target_low.npy'))
export('Fig2|2a', 'Fig2a_cells_per_target.csv',
       pd.DataFrame({'target': tn, 'cells_high': cth.astype(int), 'cells_low': ctl.astype(int),
                     'cells_total': (cth + ctl).astype(int)}),
       f'{CROP}/{{design_mats/target_names,cells_per_target_high,cells_per_target_low}}.npy',
       'cells per gene-target per gate (retention threshold >17 in >=1 gate)', index=False)

cov  = np.load(os.path.join(LM, 'cov_names.npy'), allow_pickle=True).astype(str)
feat = np.load(os.path.join(LM, 'feature_names.npy'), allow_pickle=True).astype(str)
B    = np.load(os.path.join(LM, 'EN_B_EM.npy'))
export('Fig2|2b', 'Fig2b_beta_matrix_EN_B_EM.csv',
       pd.DataFrame(B, index=pd.Index(feat, name='feature'), columns=cov),
       os.path.join(LM, 'EN_B_EM.npy'),
       'ElasticNet+EM beta matrix, 18,683 features x 231 covariates (221 gene-targets + controls/covariates)')

ss = np.load(os.path.join(PERM, 'signed_significance.npy'))
ep = np.load(os.path.join(PERM, 'empirical_P.npy'))
export('Fig2|2b', 'Fig2b_signed_significance.csv',
       pd.DataFrame(ss, index=pd.Index(feat, name='feature'), columns=cov),
       os.path.join(PERM, 'signed_significance.npy'), 'permutation-test signed significance, same shape as beta')
export('Fig2|2b', 'Fig2b_empirical_P.csv',
       pd.DataFrame(ep, index=pd.Index(feat, name='feature'), columns=cov),
       os.path.join(PERM, 'empirical_P.npy'), 'permutation-test empirical P, same shape as beta')

for f, note in [('Fig2B_target_modules.csv', '7 target modules: gene-target -> module assignment'),
                ('Fig2B_gene_programs.csv',  '9 gene programs: RNA gene -> program assignment'),
                ('Fig2C_target_module_top8.csv', 'Fig2c labels: top-8 representative gene-targets per module'),
                ('Fig2C_gene_program_top8.csv',  'Fig2c labels: top-8 representative RNA genes per program'),
                ('Fig2C_module_program_edges.csv', 'Fig2c edges: module x program mean |beta|'),
                ('Fig2d_gene_reconciliation.csv', 'Fig2d curated gene lists vs algorithmic Fig2c lists')]:
    copy('Fig2|2b-2d', os.path.join(REPO, 'figures/nature_figures/Fig2', f), f, note)

copy('Fig2|2c', os.path.join(BASE, 'DE/DE_csvs/RNA_DE_results_MAST_highE_lowE.csv'),
     'Fig2c_RNA_DE_MAST_highE_vs_lowE.csv', 'MAST DE used to pick the genes shown in the expression heatmap')
copy('Fig2|2a', os.path.join(CROP, 'cells_per_target.csv'),
     'Fig2a_cells_per_target_per_channel.csv', 'per-channel cells per gene-target (source of the .npy above)')

# ----------------------------------------------------------------- Figure S1 / Extended Data
S1 = os.path.join(REPO, 'figures/nature_figures/FigS1')
for f, note in [('feature_clusters_moi1.csv', 'MOI=1 rerun: RNA feature clusters'),
                ('target_clusters_moi1.csv', 'MOI=1 rerun: gene-target clusters (reference for module comparison)'),
                ('TableS1_CITE_antibody_panel.csv', 'CITE-seq antibody panel'),
                ('TableS3_bulk_vs_perturbme_concordance.csv', 'bulk screen vs Perturb-ME concordance table'),
                ('Reconciliation_Fig2_counts.csv', 'reconciliation of target/feature counts across Fig2 panels')]:
    copy('FigS1', os.path.join(S1, f), f, note)

for f, note in [('High_vs_Input.gene_summary.txt', 'MAGeCK gene-level, HLA-High vs Input'),
                ('Low_vs_Input.gene_summary.txt',  'MAGeCK gene-level, HLA-Low vs Input'),
                ('High_vs_CTL.gene_summary.txt',   'MAGeCK gene-level, HLA-High vs control guides'),
                ('Low_vs_CTL.gene_summary.txt',    'MAGeCK gene-level, HLA-Low vs control guides'),
                ('High_vs_Low.gene_summary.txt',   'MAGeCK gene-level, HLA-High vs HLA-Low'),
                ('CTL_vs_Input.gene_summary.txt',  'MAGeCK gene-level, control guides vs Input'),
                ('all_samples.count.txt',          'MAGeCK raw sgRNA counts, all samples'),
                ('all_samples.count_normalized.txt', 'MAGeCK normalised sgRNA counts, all samples'),
                ('all_samples.countsummary.txt',   'MAGeCK per-sample count QC summary')]:
    copy('FigS1|MAGeCK', os.path.join(S1, 'MAGeCK_bulk_screen', f), f, note)

# ----------------------------------------------------- Extended Data: guide capture
copy('FigS_guide_capture', os.path.join(REPO, 'figures/nature_figures/FigS_guide_capture/guide_capture_per_target.csv'),
     'guide_capture_per_target.csv',
     'all three panels: bulk reads, cells, guide UMIs and guide reads per gene-target per gate')

# ----------------------------------------------------------------- manifest
H5AD_ONLY = [
    ('Fig 1b', 'HLA-I protein histogram with sorted tails'),
    ('Fig 1c', 'UMAP coloured by condition'),
    ('Fig 1d', 'UMAP coloured by HLA-I protein'),
    ('Fig 1e', 'RNA vs CITE-seq protein violins'),
    ('Fig 2d', 'HLA program score by condition (adata_arrs/HLA_composite_score_50_bins.npy, 353,091 cells)'),
    ('Fig S1a i-iii', 'cells per sample / per-sample HLA ridgeline / IFNG signature'),
    ('Fig S1c', 'marker violin row'),
]

rows = ['# Source data for the Perturb-ME manuscript figures', '',
        'Every table/CSV behind a plotted panel, copied into one place.', '',
        '`COPY` = file already existed and was read/written by the plotting code (copied verbatim).',
        '`EXPORT` = panel plotted from `.npy` arrays; re-exported here as a labelled CSV.', '',
        '| Figure | Panel | File | Kind | Provenance | Contents |',
        '|---|---|---|---|---|---|']
for panel, name, kind, src, note in log:
    fig, _, pan = panel.partition('|')
    rows.append(f'| {fig} | {pan or "-"} | `{fig}/{name}` | {kind} | `{src}` | {note} |')

rows += ['', '## Panels with no intermediate table (plotted per-cell from `adata_RNA_CITE.h5ad`)', '',
         '| Panel | Contents |', '|---|---|']
rows += [f'| {p} | {d} |' for p, d in H5AD_ONLY]
rows += ['', f'`adata_RNA_CITE.h5ad` lives at `{BASE}/adata_RNA_CITE.h5ad` (353,091 cells).',
         'Regenerate this folder with `conda run -n perturbme python src/manuscript_figures/collect_source_data.py`.', '']

os.makedirs(OUT, exist_ok=True)
open(os.path.join(OUT, 'MANIFEST.md'), 'w').write('\n'.join(rows))

miss = [l for l in log if l[2] == 'MISSING']
print(f'{len(log)} entries written to {OUT}  ({len(miss)} missing)')
for m in miss:
    print('  MISSING:', m[0], m[3])
