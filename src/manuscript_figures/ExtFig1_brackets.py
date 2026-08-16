"""
Extended Data Fig 1 (CITE surface markers) with all three pairwise significance brackets.

Violins of 5 CITE surface markers across HLA-low / control / high; seaborn inner='box' gives
the black IQR bar + white median dot. Each panel is bracketed for all three pairwise contrasts
(Low-vs-Control and Control-vs-High adjacent, High-vs-Low spanning above), tested on
PER-CHANNEL mean expression (two-sided Mann-Whitney, Benjamini-Hochberg across all 15
comparisons = 5 markers x 3 pairs) -- pseudobulk, to avoid the n-inflation of testing 353k
individual cells. High-vs-Low is the robust 9-vs-9 contrast; the control comparisons use only
n = 2 channels and are floor-limited (best possible two-sided p ~ 0.036 for 2-vs-9 complete
separation). Stars are derived from the data, not hard-coded.

Run: conda run -n perturbme python src/manuscript_figures/ExtFig1_brackets.py
"""
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 11})

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
SRC = f'{REPO}/submission/20260721_Resis_package/03_Source_Data/ExtFig1_per_cell_CITE_markers.xlsx'
OUT = f'{REPO}/doc/Fig_stats_drafts'
os.makedirs(OUT, exist_ok=True)

MARKERS = [('CITE_HLA-DR', 'HLA-DR'), ('CITE_CD274_PD-L1', 'CD274 (PD-L1)'),
           ('CITE_CD47', 'CD47'), ('CITE_CD49f', 'CD49f'), ('CITE_CD58', 'CD58')]
ORDER = ['Low', 'Control', 'High']
COLORS = {'Low': '#e1812c', 'Control': '#6a8fbf', 'High': '#3a923a'}
cols = [c for c, _ in MARKERS]

print('reading', SRC, '...', flush=True)
df = pd.read_excel(SRC)

# ---- per-channel (pseudobulk) pairwise tests; two-sided Mann-Whitney on channel means,
# BH across all markers x comparisons (5 x 3 = 15) ----
ch = df.groupby(['Sample', 'Condition'], observed=True)[cols].mean().reset_index()
COMPS = [('High', 'Low'), ('High', 'Control'), ('Low', 'Control')]
recs = []
for c in cols:
    for g1, g2 in COMPS:
        a = ch.loc[ch.Condition == g1, c].values
        b = ch.loc[ch.Condition == g2, c].values
        recs.append((c, f'{g1} vs {g2}', mannwhitneyu(a, b, alternative='two-sided')[1]))
allfdr = multipletests([r[2] for r in recs], method='fdr_bh')[1]
fdrmap = {}
for (c, comp, _), f in zip(recs, allfdr):
    fdrmap.setdefault(c, {})[comp] = f


def stars(f):
    return '****' if f < 1e-4 else '***' if f < 1e-3 else '**' if f < 1e-2 else '*' if f < 0.05 else 'ns'


print('per-channel FDR (two-sided Mann-Whitney, BH across all 15 comparisons) / stars:')
for c, lab in MARKERS:
    print(' ', lab)
    for comp in ['Low vs Control', 'High vs Control', 'High vs Low']:
        f = fdrmap[c][comp]
        print(f'    {comp:18s} FDR={f:.4f}  {stars(f)}')

# ---- figure: three pairwise brackets per panel (Low-Control, Control-High adjacent;
# High-Low spanning above), placed above the FULL data range (keep the natural y-scale;
# do not truncate - the violins run to the data max as in the original) ----
fig, axes = plt.subplots(1, 5, figsize=(16, 4.2))
for ax, (col, label) in zip(axes, MARKERS):
    sns.violinplot(data=df, x='Condition', y=col, order=ORDER, hue='Condition',
                   hue_order=ORDER, palette=COLORS, legend=False, inner='box',
                   linewidth=0.8, cut=0, ax=ax)
    ax.set_title(label)
    ax.set_xlabel('')
    ax.set_ylabel('CITE expression' if col == cols[0] else '')
    lo, hi = ax.get_ylim()
    span = hi - lo
    h = span * 0.02

    def bracket(x1, x2, y, comp):
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c='black')
        ax.text((x1 + x2) / 2, y + h, stars(fdrmap[col][comp]),
                ha='center', va='bottom', fontsize=10)
    y1 = hi + span * 0.03
    y2 = y1 + span * 0.15
    bracket(0, 1, y1, 'Low vs Control')
    bracket(1, 2, y1, 'High vs Control')
    bracket(0, 2, y2, 'High vs Low')
    ax.set_ylim(top=y2 + span * 0.15)

sns.despine(fig)
for ax in axes:
    ax.grid(False)
plt.tight_layout()
plt.savefig(f'{OUT}/EDFig1_brackets.pdf', bbox_inches='tight')
plt.savefig('/gpfs/scratchfs01/site/u/wangh256/tmp/claude-2090426/'
            '-gnet-is1-p01-shares-regevlab-hanchen-Pert-PG-perturb-me/'
            '8947fcd1-ac2d-42ac-9744-b3cd90d92285/scratchpad/EDFig1_brackets_new.png',
            bbox_inches='tight', dpi=130)
plt.close()
print('saved EDFig1_brackets.pdf to', OUT)
