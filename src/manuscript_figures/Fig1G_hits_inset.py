"""
Fig 1g + "significant hits only" inset (Aviv's round-2 request).

Aviv, p4: "I would add an inset showing only significant 'hits' (ie targets enriched
in at least one of the methods gDNA or gRNA based)."  The implicit question is whether
the low genome-wide correlation (Spearman rho ~ 0.35 / 0.31) is "just noise of the
unenriched". Restricting to true hits answers it: among the 40 HLA-low hits the
reads-vs-cells concordance is strong (rho = 0.83), so the single-cell readout faithfully
tracks bulk enrichment for real targets; the low genome-wide rho is dominated by the
~19,500 non-hit targets sitting in the detection-noise floor.

Hit definition (union of the two methods, exactly as Aviv phrased it):
  bulk gDNA   = MAGeCK Low/High_vs_Input, FDR < 0.05  (40 low, 1 high)
  single-cell = PerturbME per-gene enrichment, FDR < 0.05  (11 low, 0 high; all low
                already in the bulk set)
  -> HLA-low  : 40 hits (29 bulk-only + 11 significant in BOTH methods)
  -> HLA-high : 1 hit  (VPS35, bulk-only) -> too few for a correlation inset; annotated.

Source of truth for membership: doc/Supplementary Files/TableS3_bulk_vs_perturbme_concordance.csv
Coordinates: figures/fig1g_recompute/reads_vs_cells_per_target.csv

Run: conda run -n perturbme python src/manuscript_figures/Fig1G_hits_inset.py
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
rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 12})

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
COORD = f'{REPO}/figures/fig1g_recompute/reads_vs_cells_per_target.csv'
CONC = f'{REPO}/doc/Supplementary Files/TableS3_bulk_vs_perturbme_concordance.csv'
OUT = f'{REPO}/doc/Fig_stats_drafts'
os.makedirs(OUT, exist_ok=True)

GREY = '#4a4a4a'
BULK_ONLY = '#c0392b'   # red  - significant in bulk gDNA only
BOTH = '#2c6fbb'        # blue - significant in BOTH bulk and single-cell

d = pd.read_csv(COORD)
conc = pd.read_csv(CONC)
conc['hit_low'] = conc['bulk_low_FDR_lt_0.05'] | conc['PertME_low_FDR_lt_0.05']
conc['hit_high'] = conc['bulk_high_FDR_lt_0.05'] | conc['PertME_high_FDR_lt_0.05']
low_hits = set(conc.loc[conc['hit_low'], 'target'])
low_both = set(conc.loc[conc['bulk_low_FDR_lt_0.05'] & conc['PertME_low_FDR_lt_0.05'], 'target'])
high_hits = set(conc.loc[conc['hit_high'], 'target'])

nc_all = d[~d.is_control].copy()

# ---- main two-panel figure, matching the current published layout -------------------
configs = [('Low', 'reads_low', 'cells_low', (0, 3e5), (0, 200), 40),
           ('High', 'reads_high', 'cells_high', (0, 150e3), (0, 100), 25)]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, (cond, rc, cc, xlim, ylim, thr) in zip(axes, configs):
    nc = nc_all.dropna(subset=[rc, cc])
    r, c = nc[rc].values, nc[cc].values
    ax.scatter(r, c, s=12, c=GREY, alpha=0.5, rasterized=True, linewidths=0)
    rho = spearmanr(r, c)[0]
    ax.text(0.05, 0.95, f'$\\rho$ = {rho:.2f}', transform=ax.transAxes, va='top')
    ax.set_xlabel('Reads per gene-target')
    ax.set_ylabel('Cells per gene-target')
    ax.set_title(f'HLA-{cond}')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    lab = nc[nc[cc] > thr]
    texts = [ax.text(row[rc], row[cc], row['gene'], fontsize=11, fontstyle='italic')
             for _, row in lab.iterrows()]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color=BULK_ONLY, lw=0.8))

    if cond == 'Low':
        # ---- hits-only inset ----
        hh = nc[nc['gene'].isin(low_hits)].copy()
        hh['col'] = np.where(hh['gene'].isin(low_both), BOTH, BULK_ONLY)
        rho_h = spearmanr(hh[rc], hh[cc])[0]
        _, p_h = spearmanr(hh[rc], hh[cc])
        # far bottom-right, where the main scatter is empty (high reads / low cells)
        axi = ax.inset_axes([0.60, 0.07, 0.37, 0.37])
        axi.scatter(hh[rc], hh[cc], s=18, c=hh['col'], linewidths=0, alpha=0.9)
        axi.set_xlim(xlim)   # same frame as the main HLA-low panel
        axi.set_ylim(ylim)
        axi.set_xlabel('reads / gene-target', fontsize=7.5)
        axi.set_ylabel('cells / gene-target', fontsize=7.5)
        axi.tick_params(labelsize=6.5)
        axi.text(0.04, 0.97, f'significant hits only (n={len(hh)})\n$\\rho$ = {rho_h:.2f}',
                 transform=axi.transAxes, va='top', fontsize=8)
        n_both, n_only = len(low_both), len(low_hits) - len(low_both)
        axi.scatter([], [], s=18, c=BULK_ONLY, label=f'bulk gDNA only (n={n_only})')
        axi.scatter([], [], s=18, c=BOTH, label=f'bulk + single-cell (n={n_both})')
        axi.legend(fontsize=6, loc='lower right', frameon=False, handletextpad=0.2)
        sns.despine(ax=axi)
        axi.grid(False)
        print(f'LOW hits inset: n={len(hh)}  rho={rho_h:.3f}  p={p_h:.2e}  '
              f'(bulk-only {n_only}, both {n_both})')
    else:
        vps = nc[nc['gene'] == 'VPS35']
        ax.text(0.05, 0.87,
                'only 1 significant hit (VPS35)\nin either method — no hits-only inset',
                transform=ax.transAxes, va='top', fontsize=8.5, color=GREY)
        print(f'HIGH hits: {sorted(high_hits)}  (n={len(high_hits)} - inset not defined)')

sns.despine(fig)
for ax in axes:
    ax.grid(False)
plt.tight_layout()
plt.savefig(f'{OUT}/Fig1g_hits_inset_v3.pdf', bbox_inches='tight')
plt.savefig(f'{OUT}/Fig1g_hits_inset_v3.png', bbox_inches='tight', dpi=200)
plt.close()

# ---- standalone inset (drop straight into the Illustrator figure) -------------------
figi, axi = plt.subplots(figsize=(3.4, 3.0))
hh = nc_all.dropna(subset=['reads_low', 'cells_low'])
hh = hh[hh['gene'].isin(low_hits)].copy()
hh['col'] = np.where(hh['gene'].isin(low_both), BOTH, BULK_ONLY)
rho_h = spearmanr(hh['reads_low'], hh['cells_low'])[0]
axi.scatter(hh['reads_low'], hh['cells_low'], s=30, c=hh['col'], linewidths=0, alpha=0.9)
axi.set_xlim(0, 3e5)   # same frame as the main HLA-low panel
axi.set_ylim(0, 200)
axi.set_title(f'HLA-low, significant hits only (n={len(hh)})', fontsize=9)
axi.set_xlabel('Reads per gene-target', fontsize=9)
axi.set_ylabel('Cells per gene-target', fontsize=9)
axi.text(0.05, 0.96, f'$\\rho$ = {rho_h:.2f}', transform=axi.transAxes, va='top', fontsize=10)
n_both, n_only = len(low_both), len(low_hits) - len(low_both)
axi.scatter([], [], s=30, c=BULK_ONLY, label=f'bulk gDNA only (n={n_only})')
axi.scatter([], [], s=30, c=BOTH, label=f'bulk + single-cell (n={n_both})')
axi.legend(fontsize=7, loc='lower right', frameon=False)
sns.despine(ax=axi)
axi.grid(False)
plt.tight_layout()
plt.savefig(f'{OUT}/Fig1g_low_inset_standalone.pdf', bbox_inches='tight')
plt.savefig(f'{OUT}/Fig1g_low_inset_standalone.png', bbox_inches='tight', dpi=200)
plt.close()

print('saved Fig1g_hits_inset_v3 + Fig1g_low_inset_standalone to', OUT)
