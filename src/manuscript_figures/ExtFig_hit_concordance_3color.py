"""
Extended Data figure (three-colour): per-gene bulk-gDNA reads vs single-cell counts,
with hits classified by which method(s) call them:
    bulk gDNA screen only  |  single-cell (gRNA) only  |  both

Single-cell hit set is recomputed genome-wide (src: scratchpad/sc_enrich_final.csv, from
sc_enrich3.py): per-gene binomial on Low-vs-High gate membership, genome-wide BH, requiring
>= 2 guides consistently skewed (recovers all 11 originally-reported single-cell low hits).
Bulk hit set: MAGeCK Low/High-vs-Input, FDR < 0.05 (40 low, 1 high).

Writes PDF only (per request); PNG preview goes to the scratchpad.
"""
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from scipy.stats import spearmanr
from adjustText import adjust_text

warnings.filterwarnings('ignore')
rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 11})

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
SCR = '/gpfs/scratchfs01/site/u/wangh256/tmp/claude-2090426/-gnet-is1-p01-shares-regevlab-hanchen-Pert-PG-perturb-me/8947fcd1-ac2d-42ac-9744-b3cd90d92285/scratchpad'
DRAFT = f'{REPO}/doc/Fig_stats_drafts'

GREY = '#d0d0d0'
C_BULK = '#c0392b'   # red   - bulk gDNA only
C_SC = '#2c6fbb'     # blue  - single-cell only
C_BOTH = '#6a3d9a'   # purple- both methods

coord = pd.read_csv(f'{REPO}/figures/fig1g_recompute/reads_vs_cells_per_target.csv')
sc = pd.read_csv(f'{SCR}/sc_enrich_final.csv', index_col=0)
conc = pd.read_csv(f'{REPO}/doc/Supplementary Files/TableS3_bulk_vs_perturbme_concordance.csv')
bulk = {'low': set(conc.loc[conc['bulk_low_FDR_lt_0.05'], 'target']),
        'high': set(conc.loc[conc['bulk_high_FDR_lt_0.05'], 'target'])}
scset = {'low': set(sc.index[sc['sc_low']]), 'high': set(sc.index[sc['sc_high']])}

nc = coord[~coord.is_control]

def category(gene, d):
    b, s = gene in bulk[d], gene in scset[d]
    return 'both' if (b and s) else 'bulk' if b else 'sc' if s else 'ns'

panels = [
    dict(d='low', rc='reads_low', cc='cells_low', xlim=(0, 1.5e5), ylim=(0, 100),
         xticks=[0, 5e4, 1e5, 1.5e5], yticks=[0, 25, 50, 75, 100]),
    dict(d='high', rc='reads_high', cc='cells_high', xlim=(0, 1e5), ylim=(0, 75),
         xticks=[0, 5e4, 1e5], yticks=[0, 25, 50, 75]),
]
CMAP = {'bulk': C_BULK, 'sc': C_SC, 'both': C_BOTH}

fig, axes = plt.subplots(1, 2, figsize=(15, 7))
for ax, P in zip(axes, panels):
    d = P['d']
    sub = nc.dropna(subset=[P['rc'], P['cc']]).copy()
    sub['cat'] = [category(g, d) for g in sub['gene']]
    ns = sub[sub['cat'] == 'ns']
    ax.scatter(ns[P['rc']], ns[P['cc']], s=7, c=GREY, alpha=0.5, rasterized=True, linewidths=0)
    hits = sub[sub['cat'] != 'ns']
    ax.scatter(hits[P['rc']], hits[P['cc']], s=42, c=[CMAP[c] for c in hits['cat']],
               linewidths=0, alpha=0.95, zorder=3)
    ax.set_xlim(P['xlim']); ax.set_ylim(P['ylim'])
    ax.set_xticks(P['xticks']); ax.set_yticks(P['yticks'])
    ax.set_xlabel('Reads per gene-target (bulk gDNA screen)')
    ax.set_ylabel('Cells per gene-target (single-cell)')
    n_b = int((hits['cat'] == 'bulk').sum())
    n_s = int((hits['cat'] == 'sc').sum())
    n_bo = int((hits['cat'] == 'both').sum())
    ax.set_title(f'HLA-{d}   (bulk-only {n_b}, single-cell-only {n_s}, both {n_bo})', fontsize=11)

    # rho annotations (reads vs cells) per category
    def rho(genes):
        s = sub[sub['gene'].isin(genes)]
        return spearmanr(s[P['rc']], s[P['cc']])[0] if len(s) > 2 else np.nan
    r_all = rho(set(sub['gene']))
    r_sc = rho(scset[d])
    if d == 'low':
        r_bulk = rho(bulk[d])
        txt = (f'all targets  $\\rho$ = {r_all:.2f}\nbulk hits  $\\rho$ = {r_bulk:.2f}\n'
               f'single-cell hits  $\\rho$ = {r_sc:.2f}')
    else:
        txt = f'all targets  $\\rho$ = {r_all:.2f}\nsingle-cell hits  $\\rho$ = {r_sc:.2f}'
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', fontsize=9)

    # label only hits within the (zoomed) axis range; a few high-value both-methods low hits
    # fall outside the frame and are listed in the caption (also labelled in Fig. 1g)
    lab = hits[(hits[P['rc']] <= P['xlim'][1]) & (hits[P['cc']] <= P['ylim'][1])]
    texts = [ax.text(r[P['rc']], r[P['cc']], r['gene'], fontsize=7, fontstyle='italic',
                     color=CMAP[r['cat']]) for _, r in lab.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.4))

for ax, loc in [(axes[0], 'lower right'), (axes[1], 'lower right')]:
    ax.scatter([], [], s=42, c=C_BULK, label='bulk gDNA screen only')
    ax.scatter([], [], s=42, c=C_SC, label='single-cell (gRNA) only')
    ax.scatter([], [], s=42, c=C_BOTH, label='both methods')
    ax.scatter([], [], s=7, c=GREY, label='all targets (n.s.)')
    ax.legend(loc=loc, frameon=False, fontsize=9)

sns.despine(fig)
for ax in axes:
    ax.grid(False)
plt.tight_layout()
plt.savefig(f'{DRAFT}/ExtFig_hit_concordance_3color.pdf', bbox_inches='tight')
plt.savefig(f'{SCR}/ExtFig_hit_concordance_3color.png', bbox_inches='tight', dpi=150)
plt.close()
print('saved ExtFig_hit_concordance_3color.pdf to', DRAFT)
