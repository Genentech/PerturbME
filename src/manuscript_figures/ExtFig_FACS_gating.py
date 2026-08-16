"""
Extended Data figure: FACS gating strategy for the HLA-I enrichment sort.

Aviv's comment on the promax draft: "You need to show the FACS plots/gating in a supp
figure." The listmode data is archived at
  <transfer>/GREDXFER-5990/PerturbCITE_ICR/DTA_Google_drive/20200709_PERTURBciteSORT/
one folder per Sony MA900 ($CYT = LE-MA900FP, $DATE = 09-Jul-2020), matching the Methods
("sorted on two Sony MA900s concurrently").

The instrument gate definitions were not archived with the .fcs files, so the gates here are
re-derived from the events themselves and every boundary is reported in the source data:
  cells      FSC-A and SSC-A within the 2nd-98th percentile of the recorded events
  singlets   FSC-H/FSC-A ratio within median +/- 2.5 MAD
  mKate+     mCherry-A above the 99th percentile of the no-library NTC control
  HLA gates  5th and 95th percentile of FITC-A among mKate+ singlets

Run (note LD_LIBRARY_PATH - the conda pandas needs the env's libstdc++):
  LD_LIBRARY_PATH=/home/wangh256/miniforge3/envs/perturbme/lib \
  /home/wangh256/miniforge3/envs/perturbme/bin/python \
  src/manuscript_figures/ExtFig_FACS_gating.py
"""
import os

import fcsparser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
from matplotlib.lines import Line2D

rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['pdf.use14corefonts'] = True

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
FACS = ('/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/GREDXFER-5990/PerturbCITE_ICR/'
        'DTA_Google_drive/20200709_PERTURBciteSORT')
S1 = os.path.join(FACS, 'MZ_070920_ICR47_GW_10x_waffles_mike/Sample Group - 1')   # sorter 1
S2 = os.path.join(FACS, 'MZ_070920_ICR47_GW_10x_nala_maryann/Sample Group - 1')   # sorter 2
OUT = os.path.join(REPO, 'figures/nature_figures/FigS_FACS_gating')
os.makedirs(OUT, exist_ok=True)

LOW, CTL, HIGH = '#e1812c', '#3274a1', '#3a923a'
GREY = '#9a9a9a'

SAMPLES = {
    'GW+Ctl sorter 1': f'{S1}/GW+Ctl IFN2 CITE_[5 mL Tubes] Data Source - 1.fcs',
    'GW+Ctl sorter 2': f'{S2}/GW_IFN2_CITE_[5 mL Tubes] Data Source - 1.fcs',
    'Ctl library only': f'{S1}/Ctl IFN2 CITE_Data Source - 1.fcs',
    'NTC (no library)': f'{S1}/NTC IFN2_Data Source - 1.fcs',
    'NT IFN0': f'{S2}/NT_IFN0_HLA-FITC_Data Source - 1.fcs',
    'NT IFN2': f'{S2}/NT_IFN2_HLA-FITC _Data Source - 1.fcs',
}


def load(path):
    _, d = fcsparser.parse(path, reformat_meta=True)
    return d


def cell_gate(d):
    """Main scatter population: 2nd-98th percentile in both FSC-A and SSC-A."""
    lo_f, hi_f = np.percentile(d['FSC-A'], [2, 98])
    lo_s, hi_s = np.percentile(d['SSC-A'], [2, 98])
    m = d['FSC-A'].between(lo_f, hi_f) & d['SSC-A'].between(lo_s, hi_s)
    return m.to_numpy(), dict(fsc_lo=lo_f, fsc_hi=hi_f, ssc_lo=lo_s, ssc_hi=hi_s)


def singlet_gate(d, cells):
    """FSC-H/FSC-A ratio within median +/- 2.5 MAD of the cell-gated events."""
    r = (d['FSC-H'] / d['FSC-A']).to_numpy()
    med = np.median(r[cells])
    mad = np.median(np.abs(r[cells] - med))
    lo, hi = med - 2.5 * mad, med + 2.5 * mad
    return cells & (r >= lo) & (r <= hi), dict(ratio_lo=lo, ratio_hi=hi)


data = {k: load(v) for k, v in SAMPLES.items()}

# mKate+ threshold from the no-library control
ntc = data['NTC (no library)']
ntc_cells, _ = cell_gate(ntc)
ntc_singlets, _ = singlet_gate(ntc, ntc_cells)
MKATE_THR = float(np.percentile(ntc['mCherry-A'].to_numpy()[ntc_singlets], 99))

rows, gated = [], {}
for name, d in data.items():
    cells, cg = cell_gate(d)
    singlets, sg = singlet_gate(d, cells)
    mkate = singlets & (d['mCherry-A'].to_numpy() > MKATE_THR)
    # HLA statistics are quoted on the population the gate would actually act on: mKate2+
    # singlets for library-containing samples, all singlets for the no-library controls
    # (where mKate2+ is ~1% and its HLA distribution is meaningless).
    has_library = mkate.sum() / max(singlets.sum(), 1) > 0.5
    pop = mkate if has_library else singlets
    hla = d['FITC-A'].to_numpy()[pop]
    p5, p95 = (float(np.percentile(hla, 5)), float(np.percentile(hla, 95))) if hla.size else (np.nan, np.nan)
    gated[name] = dict(cells=cells, singlets=singlets, mkate=mkate, p5=p5, p95=p95, **cg, **sg)
    rows.append({
        'sample': name, 'events_recorded': len(d),
        'events_in_cell_gate': int(cells.sum()), 'pct_cells': round(100 * cells.mean(), 2),
        'events_in_singlet_gate': int(singlets.sum()), 'pct_singlets': round(100 * singlets.mean(), 2),
        'events_mKate_pos': int(mkate.sum()),
        'pct_mKate_pos_of_singlets': round(100 * mkate.sum() / max(singlets.sum(), 1), 2),
        'mKate_threshold_mCherry_A': round(MKATE_THR, 1),
        'HLA_stats_population': 'mKate2+ singlets' if has_library else 'all singlets (no library)',
        'HLA_FITC_A_5th_pct': round(p5, 1), 'HLA_FITC_A_median': round(float(np.median(hla)), 1) if hla.size else np.nan,
        'HLA_FITC_A_95th_pct': round(p95, 1),
        'FSC_A_gate_low': round(cg['fsc_lo'], 1), 'FSC_A_gate_high': round(cg['fsc_hi'], 1),
        'SSC_A_gate_low': round(cg['ssc_lo'], 1), 'SSC_A_gate_high': round(cg['ssc_hi'], 1),
        'FSC_H_over_FSC_A_low': round(sg['ratio_lo'], 4), 'FSC_H_over_FSC_A_high': round(sg['ratio_hi'], 4),
    })
summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT, 'ExtFig_FACS_gating_statistics.csv'), index=False)

# per-event table for the plotted samples (source data)
ev = []
for name, d in data.items():
    g = gated[name]
    ev.append(pd.DataFrame({
        'sample': name,
        'FSC_A': d['FSC-A'].to_numpy().round(1), 'FSC_H': d['FSC-H'].to_numpy().round(1),
        'SSC_A': d['SSC-A'].to_numpy().round(1),
        'mCherry_A_mKate2': d['mCherry-A'].to_numpy().round(1),
        'FITC_A_HLA_ABC': d['FITC-A'].to_numpy().round(1),
        'in_cell_gate': g['cells'], 'in_singlet_gate': g['singlets'], 'mKate_positive': g['mkate'],
    }))
pd.concat(ev, ignore_index=True).to_csv(
    os.path.join(OUT, 'ExtFig_FACS_gating_events.csv'), index=False)

# ------------------------------------------------------------------ figure
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
ref = 'GW+Ctl sorter 1'
d, g = data[ref], gated[ref]


def logbins(v, n=120):
    v = v[v > 0]
    return np.logspace(np.log10(np.percentile(v, 0.1)), np.log10(np.percentile(v, 99.9)), n)


# a) FSC-A vs SSC-A
ax = axes[0, 0]
ax.hexbin(d['FSC-A'], d['SSC-A'], gridsize=90, bins='log', cmap='Greys', mincnt=1,
          extent=(0, np.percentile(d['FSC-A'], 99.5), 0, np.percentile(d['SSC-A'], 99.5)))
ax.add_patch(plt.Rectangle((g['fsc_lo'], g['ssc_lo']), g['fsc_hi'] - g['fsc_lo'],
                           g['ssc_hi'] - g['ssc_lo'], fill=False, ec='#c0392b', lw=1.4))
ax.set_xlim(0, np.percentile(d['FSC-A'], 99.5)); ax.set_ylim(0, np.percentile(d['SSC-A'], 99.5))
ax.set_xlabel('FSC-A'); ax.set_ylabel('SSC-A')
ax.set_title(f"a  Cells\n{100 * g['cells'].mean():.1f}% of events", loc='left', fontsize=11)

# b) FSC-A vs FSC-H
ax = axes[0, 1]
c = g['cells']
ax.hexbin(d['FSC-A'][c], d['FSC-H'][c], gridsize=90, bins='log', cmap='Greys', mincnt=1)
xs = np.array([0, np.percentile(d['FSC-A'][c], 99.5)])
for k, ls in [(g['ratio_lo'], '--'), (g['ratio_hi'], '--')]:
    ax.plot(xs, k * xs, ls, color='#c0392b', lw=1.2)
ax.set_xlim(*xs); ax.set_ylim(0, np.percentile(d['FSC-H'][c], 99.5))
ax.set_xlabel('FSC-A'); ax.set_ylabel('FSC-H')
ax.set_title(f"b  Singlets\n{100 * g['singlets'].sum() / max(c.sum(), 1):.1f}% of cells",
             loc='left', fontsize=11)

# c) mKate2
ax = axes[0, 2]
bins = logbins(np.concatenate([data['NTC (no library)']['mCherry-A'].to_numpy(),
                               d['mCherry-A'].to_numpy()]))
ax.hist(data['NTC (no library)']['mCherry-A'][gated['NTC (no library)']['singlets']],
        bins=bins, color=GREY, alpha=.75, label='NTC (no library)')
ax.hist(d['mCherry-A'][g['singlets']], bins=bins, histtype='step', color=CTL, lw=1.6,
        label='genome-wide + control')
ax.axvline(MKATE_THR, ls='--', color='#c0392b', lw=1.2)
ax.set_xscale('log'); ax.set_xlabel('mKate2 (mCherry-A)'); ax.set_ylabel('Events')
ax.legend(frameon=False, fontsize=8, loc='upper left')
ax.set_title(f"c  mKate2$^+$ (library$^+$)\n{100 * g['mkate'].sum() / max(g['singlets'].sum(), 1):.1f}%"
             f" of singlets", loc='left', fontsize=11)

# d) HLA gates
ax = axes[0, 3]
hla = d['FITC-A'].to_numpy()[g['mkate']]
bins = logbins(hla)
n, e, patches = ax.hist(hla, bins=bins, color='#3a3a3a')
centres = (e[:-1] + e[1:]) / 2
for patch, x in zip(patches, centres):
    patch.set_facecolor(LOW if x <= g['p5'] else HIGH if x >= g['p95'] else '#3a3a3a')
for v in (g['p5'], g['p95']):
    ax.axvline(v, ls='--', color='gray', lw=1)
ax.set_xscale('log'); ax.set_xlabel('HLA-A,B,C (FITC-A)'); ax.set_ylabel('Events')
ax.set_title('d  Sort gates on mKate2$^+$ singlets\nbottom 5% / top 5%', loc='left', fontsize=11)
ax.legend(handles=[Line2D([], [], color=LOW, lw=6, label='HLA-low (bottom 5%)'),
                   Line2D([], [], color=HIGH, lw=6, label='HLA-high (top 5%)')],
          frameon=False, fontsize=8, loc='upper left')

# e) sorter concordance
ax = axes[1, 0]
allv = np.concatenate([data[s]['FITC-A'].to_numpy()[gated[s]['mkate']]
                       for s in ('GW+Ctl sorter 1', 'GW+Ctl sorter 2')])
bins = logbins(allv)
for s, col in [('GW+Ctl sorter 1', CTL), ('GW+Ctl sorter 2', '#7a5090')]:
    ax.hist(data[s]['FITC-A'][gated[s]['mkate']], bins=bins, histtype='step', lw=1.6,
            density=True, color=col, label=f"{s} (median {np.median(data[s]['FITC-A'][gated[s]['mkate']]):,.0f})")
ax.set_xscale('log'); ax.set_xlabel('HLA-A,B,C (FITC-A)'); ax.set_ylabel('Density')
ax.legend(frameon=False, fontsize=8)
ax.set_title('e  Two MA900 sorters, run concurrently', loc='left', fontsize=11)

# f) IFNg induction
ax = axes[1, 1]
allv = np.concatenate([data[s]['FITC-A'].to_numpy() for s in ('NT IFN0', 'NT IFN2')])
bins = logbins(allv)
for s, col, lab in [('NT IFN0', GREY, 'untreated'), ('NT IFN2', '#c0392b', 'IFN$\\gamma$ 2 ng/mL, 24 h')]:
    m = gated[s]['singlets']
    ax.hist(data[s]['FITC-A'][m], bins=bins, histtype='stepfilled', alpha=.55, density=True,
            color=col, label=f"{lab} (median {np.median(data[s]['FITC-A'][m]):,.0f})")
ax.set_xscale('log'); ax.set_xlabel('HLA-A,B,C (FITC-A)'); ax.set_ylabel('Density')
ax.legend(frameon=False, fontsize=8)
ax.set_title('f  IFN$\\gamma$ induction of surface HLA-I', loc='left', fontsize=11)

# g) control library only
ax = axes[1, 2]
s = 'Ctl library only'
hla_c = data[s]['FITC-A'].to_numpy()[gated[s]['mkate']]
bins = logbins(np.concatenate([hla_c, hla]))
ax.hist(hla_c, bins=bins, color=CTL, alpha=.8, density=True, label='control library (all sorted)')
ax.hist(hla, bins=bins, histtype='step', color='#3a3a3a', lw=1.4, density=True,
        label='genome-wide + control')
ax.set_xscale('log'); ax.set_xlabel('HLA-A,B,C (FITC-A)'); ax.set_ylabel('Density')
ax.legend(frameon=False, fontsize=8)
ax.set_title('g  Control-library arm\nentire distribution collected', loc='left', fontsize=11)

# h) events retained per gate
ax = axes[1, 3]
steps = ['Recorded', 'Cells', 'Singlets', 'mKate2$^+$']
for s, col in [('GW+Ctl sorter 1', CTL), ('GW+Ctl sorter 2', '#7a5090')]:
    gg = gated[s]
    tot = len(data[s])
    vals = [100, 100 * gg['cells'].mean(), 100 * gg['singlets'].mean(), 100 * gg['mkate'].sum() / tot]
    ax.plot(steps, vals, 'o-', color=col, lw=1.6, ms=5, label=s)
ax.set_ylabel('% of recorded events'); ax.set_ylim(0, 105)
ax.legend(frameon=False, fontsize=8)
ax.set_title('h  Gating hierarchy yield', loc='left', fontsize=11)

sns.despine(fig=fig)
for a in axes.flat:
    a.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'ExtFig_FACS_gating.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(OUT, 'ExtFig_FACS_gating.png'), bbox_inches='tight', dpi=300)

print(summary.to_string(index=False))
print(f'\nwrote {OUT}')
