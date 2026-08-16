"""
Fig 1g (updated): per-target single-cell representation vs bulk gDNA reads, with a
decile-binned inset showing the monotonic concordance.

Fixes the text<->figure mismatch Aviv flagged: the per-target Pearson (0.79/0.82)
reported in the text was inflated by the two pooled control pseudo-targets
(NO_SITE / ONE_NON-GENE_SITE), which carry ~241 guides each and sit 6-16x beyond any
real target. The robust statistics are per-target Spearman rho = 0.31/0.35 (controls
excluded) and a strong decile-binned trend (Spearman = 1.0, Pearson = 0.94/0.99).

Run: conda run -n perturbme python src/manuscript_figures/Fig1G_reads_vs_cells_deciles.py
"""
import os, warnings, numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
from matplotlib import rcParams
from scipy.stats import spearmanr, pearsonr
from adjustText import adjust_text
warnings.filterwarnings('ignore')
rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'pdf.use14corefonts':True,'font.size':12})

CSV='figures/fig1g_recompute/reads_vs_cells_per_target.csv'
OUT='figures/nature_figures/Fig1'; os.makedirs(OUT,exist_ok=True)
ACCENT='#c0392b'   # red, as in Fig1F/G label arrows
d=pd.read_csv(CSV)

configs=[('Low','reads_low','cells_low',(0,3e5),(0,200),40),
         ('High','reads_high','cells_high',(0,150e3),(0,100),25)]

fig,axes=plt.subplots(1,2,figsize=(14,6))
for ax,(cond,rc,cc,xlim,ylim,thr) in zip(axes,configs):
    full=d.dropna(subset=[rc,cc]); nc=full[~full.is_control]
    r,c=nc[rc].values,nc[cc].values
    ax.scatter(r,c,s=12,c='#4a4a4a',alpha=0.5,rasterized=True,linewidths=0)
    rho=spearmanr(r,c)[0]
    ax.text(0.05,0.95,f'$\\rho$ = {rho:.2f}  (per target)',transform=ax.transAxes,va='top')
    ax.set_xlabel('Reads per target'); ax.set_ylabel('Cells per target')
    ax.set_title(f'HLA-{cond}'); ax.set_xlim(xlim); ax.set_ylim(ylim)
    # gene labels for well-detected targets (as in published 1g)
    lab=nc[nc[cc]>thr]
    texts=[ax.text(row[rc],row[cc],row['gene'],fontsize=11,fontstyle='italic') for _,row in lab.iterrows()]
    if texts: adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',color=ACCENT,lw=0.8))

    # ---- decile-binned inset ----
    nc=nc.copy(); nc['dec']=pd.qcut(nc[rc].rank(method='first'),10,labels=False)
    g=nc.groupby('dec').agg(mr=(rc,'mean'),mc=(cc,'mean'))
    pr=pearsonr(g['mr'],g['mc'])[0]; sr=spearmanr(g['mr'],g['mc'])[0]
    ax_in=ax.inset_axes([0.56,0.12,0.40,0.42])
    ax_in.plot(g['mr'],g['mc'],'-o',color=ACCENT,ms=5,lw=1.5)
    ax_in.set_title('binned (deciles)',fontsize=10)
    ax_in.set_xlabel('mean reads',fontsize=9); ax_in.set_ylabel('mean cells',fontsize=9)
    ax_in.tick_params(labelsize=8)
    ax_in.text(0.05,0.95,f'$r$ = {pr:.2f}\n$\\rho$ = {sr:.1f}',transform=ax_in.transAxes,va='top',fontsize=9)
    sns.despine(ax=ax_in); ax_in.grid(False)

sns.despine(fig)
for ax in axes: ax.grid(False)
plt.tight_layout()
for fmt in ('pdf','png'):
    plt.savefig(os.path.join(OUT,f'Fig1G_reads_vs_cells_with_deciles.{fmt}'),bbox_inches='tight',dpi=300)
print('saved Fig1G_reads_vs_cells_with_deciles to',OUT)
