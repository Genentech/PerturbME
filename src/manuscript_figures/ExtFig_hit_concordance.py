"""
Extended Data figure: bulk-vs-single-cell concordance for the significant hits.

Promotes Aviv's requested Fig 1g "significant hits only" inset to its own Extended Data
figure, where there is room to label every hit and to show both sort directions.

Hit definition (single, rigorous set): the gene-targets significant in the bulk gDNA
screen (MAGeCK Low/High-vs-Input, FDR < 0.05) -> 40 HLA-low, 1 HLA-high (VPS35). Every
target the single-cell enrichment flags is already in this bulk set (11 of the 40 low
hits; 0 high), so the bulk set is also the union of the two methods among tested targets.
We do NOT colour-code by method here: the single-cell FDR was only ever computed for the
bulk hits, so we cannot claim a complete "single-cell-only" category. The 11 co-significant
low hits are stated in the caption instead.

Panel a (HLA-low): 40 hits; among them Spearman rho = 0.83 vs rho = 0.35 across all 19,561
targets -- the modest genome-wide correlation reflects dilution by non-hit targets in the
detection-noise floor, not a noisy single-cell readout.
Panel b (HLA-high): a single hit (VPS35), highlighted against the full distribution to make
the 40-vs-1 asymmetry explicit.

Hit membership: doc/Supplementary Files/TableS3_bulk_vs_perturbme_concordance.csv
Coordinates:    figures/fig1g_recompute/reads_vs_cells_per_target.csv

Run: conda run -n perturbme python src/manuscript_figures/ExtFig_hit_concordance.py
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
COORD = f'{REPO}/figures/fig1g_recompute/reads_vs_cells_per_target.csv'
CONC = f'{REPO}/doc/Supplementary Files/TableS3_bulk_vs_perturbme_concordance.csv'
NFDIR = f'{REPO}/figures/nature_figures/FigS_hit_concordance'
DRAFT = f'{REPO}/doc/Fig_stats_drafts'
os.makedirs(NFDIR, exist_ok=True)
os.makedirs(DRAFT, exist_ok=True)

GREY = '#c8c8c8'   # all targets (backdrop)
HIT = '#c0392b'    # significant bulk-screen hit (Fig 1F/G accent red)

d = pd.read_csv(COORD)
conc = pd.read_csv(CONC)
low_hits = set(conc.loc[conc['bulk_low_FDR_lt_0.05'], 'target'])
high_hits = set(conc.loc[conc['bulk_high_FDR_lt_0.05'], 'target'])
also_sc_low = sorted(conc.loc[conc['bulk_low_FDR_lt_0.05'] & conc['PertME_low_FDR_lt_0.05'], 'target'])

nc = d[~d.is_control]

panels = [
    dict(cond='low', rc='reads_low', cc='cells_low', xlim=(0, 3e5), ylim=(0, 200),
         xticks=[0, 1e5, 2e5, 3e5], yticks=[0, 50, 100, 150, 200], hits=low_hits),
    dict(cond='high', rc='reads_high', cc='cells_high', xlim=(0, 15e4), ylim=(0, 100),
         xticks=[0, 5e4, 1e5, 1.5e5], yticks=[0, 25, 50, 75, 100], hits=high_hits),
]

fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))
for ax, P in zip(axes, panels):
    sub = nc.dropna(subset=[P['rc'], P['cc']])
    ax.scatter(sub[P['rc']], sub[P['cc']], s=8, c=GREY, alpha=0.5, rasterized=True, linewidths=0)
    hh = sub[sub['gene'].isin(P['hits'])].copy()
    ax.scatter(hh[P['rc']], hh[P['cc']], s=45, c=HIT, linewidths=0, alpha=0.95, zorder=3)
    ax.set_xlim(P['xlim'])
    ax.set_ylim(P['ylim'])
    ax.set_xticks(P['xticks'])
    ax.set_yticks(P['yticks'])
    ax.set_xlabel('Reads per gene-target (bulk gDNA screen)')
    ax.set_ylabel('Cells per gene-target (single-cell)')
    ax.set_title(f'HLA-{P["cond"]}')

    rho_all = spearmanr(sub[P['rc']], sub[P['cc']])[0]
    if len(hh) >= 3:
        rho_h, p_h = spearmanr(hh[P['rc']], hh[P['cc']])
        ax.text(0.04, 0.96, f'hits: $\\rho$ = {rho_h:.2f}  (n = {len(hh)})\n'
                            f'all targets: $\\rho$ = {rho_all:.2f}',
                transform=ax.transAxes, va='top', fontsize=10.5)
        print(f'HLA-{P["cond"]}: n_hits={len(hh)}  rho_hits={rho_h:.3f}  p={p_h:.2e}  rho_all={rho_all:.3f}')
    else:
        ax.text(0.04, 0.96, f'{len(hh)} significant hit ({", ".join(sorted(P["hits"]))})\n'
                            f'all targets: $\\rho$ = {rho_all:.2f}',
                transform=ax.transAxes, va='top', fontsize=10.5)
        print(f'HLA-{P["cond"]}: n_hits={len(hh)} ({sorted(P["hits"])})  rho_all={rho_all:.3f}')

    texts = [ax.text(r[P['rc']], r[P['cc']], r['gene'], fontsize=8.5, fontstyle='italic')
             for _, r in hh.iterrows()]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#888888', lw=0.6))

axes[0].scatter([], [], s=45, c=HIT, label='significant hit (bulk gDNA, FDR < 0.05)')
axes[0].scatter([], [], s=8, c=GREY, label='all targets (n.s.)')
axes[0].legend(loc='lower right', frameon=False, fontsize=9)

sns.despine(fig)
for ax in axes:
    ax.grid(False)
plt.tight_layout()
for stem in (f'{NFDIR}/ExtFig_hit_concordance', f'{DRAFT}/ExtFig_hit_concordance'):
    plt.savefig(stem + '.pdf', bbox_inches='tight')
    plt.savefig(stem + '.png', bbox_inches='tight', dpi=200)
plt.close()
print(f'also significant in single-cell (n={len(also_sc_low)}):', ', '.join(also_sc_low))
print('saved ExtFig_hit_concordance (2 panels, single-colour) to', NFDIR, 'and', DRAFT)
