"""
Extended Data panel (Aviv's request): distribution of # cells per guide and # cells per
gene-target in the phenotype-enriched Perturb-ME data (MOI=1, 277,152 cells).

Supports the Results sentence: a conventional genome-scale Perturb-seq at this depth would
give ~4 sgRNA-singlet cells per guide (277,152 / ~60,000-guide library, uniform sampling),
whereas phenotype-based enrichment recovers enriched guides in tens of singlet cells, the
strongest guides in ~100 cells (max 96; ~21x uniform) and the strongest gene-targets in
~175 cells (~12x uniform). All counts are sgRNA-singlet cells (MOI=1), matching the ~4
baseline.

Two panels: (a) cells per guide (sgRNA), (b) cells per gene-target. Targeting guides/genes
only (NO_SITE / ONE_NON-GENE_SITE controls excluded). Log-spaced bins, log x-axis.

Writes PDF only (per request); PNG preview -> scratchpad.
"""
import os
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator, ScalarFormatter, NullLocator

warnings.filterwarnings('ignore')
rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 11})

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
DRAFT = f'{REPO}/doc/Fig_stats_drafts'
NFDIR = f'{REPO}/figures/nature_figures/FigS_cells_per_guide_gene'   # canonical output
SCR = '/gpfs/scratchfs01/site/u/wangh256/tmp/claude-2090426/-gnet-is1-p01-shares-regevlab-hanchen-Pert-PG-perturb-me/8947fcd1-ac2d-42ac-9744-b3cd90d92285/scratchpad'
UNIFORM_GUIDE = 277152 / 60000    # ~4.6 cells/guide expected under uniform sampling
UNIFORM_GENE = 277152 / 19770     # ~14 cells/gene expected under uniform sampling
os.makedirs(NFDIR, exist_ok=True)

a = ad.read_h5ad(f'{REPO}/data/20251116_processed_1moi.h5ad', backed='r')
obs = a.obs[['sgRNA']].copy()
obs['sgRNA'] = obs['sgRNA'].astype(str)
obs['gene'] = obs['sgRNA'].str.rsplit('_', n=1).str[0]
tailnum = obs['sgRNA'].str.rsplit('_', n=1).str[-1]
ctrl = (~tailnum.str.isdigit()) | obs['sgRNA'].str.upper().str.contains('NO_SITE|NON-GENE')
tg = obs[~ctrl]

# distribution tables (kept with names, so the source-data builder can reuse them exactly)
guide_tbl = tg['sgRNA'].value_counts().rename_axis('sgRNA').reset_index(name='n_cells')
guide_tbl['gene'] = guide_tbl['sgRNA'].str.rsplit('_', n=1).str[0]
guide_tbl = guide_tbl[['sgRNA', 'gene', 'n_cells']].sort_values('n_cells', ascending=False)
gene_tbl = tg['gene'].value_counts().rename_axis('gene').reset_index(name='n_cells') \
    .sort_values('n_cells', ascending=False)
guide_tbl.to_csv(f'{NFDIR}/cells_per_guide.csv', index=False)
gene_tbl.to_csv(f'{NFDIR}/cells_per_gene.csv', index=False)

cpg = guide_tbl['n_cells'].values                # cells per guide
cpgene = gene_tbl['n_cells'].values              # cells per gene-target


def stat(x, name, uni):
    print(f'{name}: n={len(x)}  median={np.median(x):.0f}  mean={x.mean():.1f}  '
          f'Q1={np.percentile(x,25):.0f}  Q3={np.percentile(x,75):.0f}  max={x.max():.0f}  '
          f'uniform={uni:.1f}  fold(max/uniform)={x.max()/uni:.0f}x  >=50: {(x>=50).sum()}')
    return np.median(x), x.mean(), int(x.max())


med_g, mean_g, max_g = stat(cpg, 'cells/guide ', UNIFORM_GUIDE)
med_t, mean_t, max_t = stat(cpgene, 'cells/gene  ', UNIFORM_GENE)

panels = [
    ('a', cpg, 'Cells per guide (sgRNA)', 'Number of guides', med_g, max_g, UNIFORM_GUIDE, '#4a6fa5'),
    ('b', cpgene, 'Cells per gene-target', 'Number of gene-targets', med_t, max_t, UNIFORM_GENE, '#7b3f6f'),
]
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
for ax, (tag, x, title, ylab, med, mx, uni, col) in zip(axes, panels):
    bins = np.logspace(0, np.log10(x.max() + 1), 45)
    ax.hist(x, bins=bins, color=col, alpha=0.85, edgecolor='white', linewidth=0.3)
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(FixedLocator([1, 10, 100]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.ticklabel_format(axis='x', style='plain')
    ax.axvline(med, ls='--', c='gray', lw=1, label=f'median = {med:.0f}')
    ax.axvline(uni, ls=':', c='#999999', lw=1.2, label=f'uniform expectation ~ {uni:.0f}')
    ax.text(0.05, 0.96, f'n = {len(x):,}\nmax {mx:,}  (~{mx/uni:.0f}x uniform)',
            transform=ax.transAxes, va='top', ha='left', fontsize=9)
    ax.legend(fontsize=8, frameon=False, loc='upper right')
    ax.set_xlabel('Cells (sgRNA-singlets, MOI=1)')
    ax.set_ylabel(ylab)
    ax.set_title(title, fontweight='bold')
    ax.text(-0.16, 1.04, tag, transform=ax.transAxes, fontsize=15, fontweight='bold', va='bottom')

sns.despine(fig)
for ax in axes:
    ax.grid(False)
plt.tight_layout()
for d in (NFDIR, DRAFT):                          # canonical + draft copies
    plt.savefig(f'{d}/ExtFig_cells_per_guide_gene.pdf', bbox_inches='tight')
plt.savefig(f'{SCR}/ExtFig_cells_per_guide_gene.png', bbox_inches='tight', dpi=150)
plt.close()
print('\nsaved ExtFig_cells_per_guide_gene.pdf to', NFDIR, 'and', DRAFT)
print('saved cells_per_guide.csv / cells_per_gene.csv to', NFDIR)
