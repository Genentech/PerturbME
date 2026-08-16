"""
Fig 1e (HLA-I surface protein + mRNA) with all three pairwise significance brackets.

Two violins -- HLA-A/B/C surface protein (CITE-seq) and HLA-A mRNA (scRNA-seq) -- across
HLA-low / control / high. Each panel is bracketed for all three pairwise contrasts
(Low-vs-Control and Control-vs-High adjacent, High-vs-Low spanning above), tested on
PER-CHANNEL mean expression (two-sided Mann-Whitney, BH across all 6 comparisons =
2 measures x 3 pairs) -- pseudobulk, avoiding the n-inflation of 353k individual cells.
High-vs-Low is the robust 9-vs-9 contrast; the control comparisons use only n = 2 channels
and are floor-limited (best possible two-sided p ~ 0.036 for 2-vs-9 complete separation).
Stars are derived from the data, not hard-coded.

Run: conda run -n perturbme python src/manuscript_figures/Fig1e_brackets.py
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
SRC = f'{REPO}/submission/20260721_Resis_package/03_Source_Data/Fig1c-1e_per_cell_UMAP_protein_mRNA.xlsx'
OUT = f'{REPO}/doc/Fig_stats_drafts'
SCR = ('/gpfs/scratchfs01/site/u/wangh256/tmp/claude-2090426/'
       '-gnet-is1-p01-shares-regevlab-hanchen-Pert-PG-perturb-me/'
       '8947fcd1-ac2d-42ac-9744-b3cd90d92285/scratchpad')
os.makedirs(OUT, exist_ok=True)

MEASURES = [('CITE_HLA_ABC_protein_log', 'HLA-A,B,C protein (CITE-seq)'),
            ('HLA_A_mRNA_log', 'HLA-A mRNA (scRNA-seq)')]
ORDER = ['Low', 'Control', 'High']
COLORS = {'Low': '#e1812c', 'Control': '#6a8fbf', 'High': '#3a923a'}
cols = [m for m, _ in MEASURES]

print('reading', SRC, '...', flush=True)
df = pd.read_excel(SRC)

# ---- per-channel (pseudobulk) tests; two-sided Mann-Whitney on channel means, BH across ALL
# comparisons run (2 measures x 3 pairs = 6). All three pairwise contrasts are displayed.
ch = df.groupby(['Sample', 'Condition'], observed=True)[cols].mean().reset_index()
COMPS = [('High', 'Low'), ('High', 'Control'), ('Low', 'Control')]
recs = []
for m in cols:
    for g1, g2 in COMPS:
        a = ch.loc[ch.Condition == g1, m].values
        b = ch.loc[ch.Condition == g2, m].values
        recs.append((m, f'{g1} vs {g2}', mannwhitneyu(a, b, alternative='two-sided')[1]))
allfdr = multipletests([r[2] for r in recs], method='fdr_bh')[1]
fdrmap = {}
for (m, comp, _), f in zip(recs, allfdr):
    fdrmap.setdefault(m, {})[comp] = f


def stars(f):
    return '****' if f < 1e-4 else '***' if f < 1e-3 else '**' if f < 1e-2 else '*' if f < 0.05 else 'ns'


print('per-channel FDR (two-sided Mann-Whitney, BH across all 6 comparisons) / stars:')
for m, lab in MEASURES:
    print(' ', lab)
    for comp in ['Low vs Control', 'High vs Control', 'High vs Low']:
        f = fdrmap[m][comp]
        print(f'    {comp:18s} FDR={f:.4f}  {stars(f)}')

# ---- figure: three pairwise brackets per panel (Low-Control, Control-High adjacent;
# High-Low spanning above) ----
fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.8))
for ax, (col, ylab) in zip(axes, MEASURES):
    sns.violinplot(data=df, x='Condition', y=col, order=ORDER, hue='Condition',
                   hue_order=ORDER, palette=COLORS, legend=False, inner='box',
                   linewidth=0.8, cut=0, ax=ax)
    ax.set_xlabel('')
    ax.set_ylabel(ylab)
    lo, hi = ax.get_ylim()
    span = hi - lo
    h = span * 0.02

    def bracket(x1, x2, y, comp):
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c='black')
        ax.text((x1 + x2) / 2, y + h, stars(fdrmap[col][comp]),
                ha='center', va='bottom', fontsize=11)
    y1 = hi + span * 0.03
    y2 = y1 + span * 0.14
    bracket(0, 1, y1, 'Low vs Control')
    bracket(1, 2, y1, 'High vs Control')
    bracket(0, 2, y2, 'High vs Low')
    ax.set_ylim(top=y2 + span * 0.14)

sns.despine(fig)
for ax in axes:
    ax.grid(False)
plt.tight_layout()
plt.savefig(f'{OUT}/Fig1e_brackets.pdf', bbox_inches='tight')
plt.savefig(f'{SCR}/Fig1e_brackets_new.png', bbox_inches='tight', dpi=130)
plt.close()
print('saved Fig1e_brackets.pdf to', OUT)
