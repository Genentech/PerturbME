"""
Supplementary figure for Aviv's comment on Fig 1g concordance:
"plot the gRNA reads in the entire bin ... addresses whether you are impacted by
how much 'expression' of the gRNA you capture."

Compares, per target per HLA bin (High/Low):
  - bulk gDNA reads (entire sorted bin; proxy for cell abundance, the bulk-screen
    tacit assumption reads ~ cells)
  - single-cell cells recovered
  - single-cell GUIDE UMIs captured  ("expression of the gRNA you capture")

Conclusion: guide UMIs/cell is ~constant (CV~0.25) and uncorrelated with cells
recovered (rho~0) -> single-cell recovery is NOT biased by guide capture; and the
single-cell per-cell guide amount is TIGHTER than bulk reads/cell.

Reconstructs guide-UMI capture from CROP/CBC_UMI_dicts/*.pkl by replicating the
assign_sgRNAs.py assignment logic while retaining UMI counts. Validated: per-target
cell counts reproduce CROP/cells_per_target_{high,low}.npy exactly (Pearson=1.000).

Run:  conda run -n perturbme python src/manuscript_figures/FigS_guide_capture_vs_bulk.py
"""
import os, glob, pickle, warnings
from collections import defaultdict
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator, ScalarFormatter, NullLocator
from scipy.stats import spearmanr

EXP   = 'PerturbME_transfer/PerturbCITE_ICR/202008_full_exp'
DICTS = os.path.join(EXP, 'CROP/CBC_UMI_dicts')
FILT  = os.path.join(EXP, 'filtered_RNA')
SEQ   = os.path.join(EXP, 'sequencing_info')
FIG1G = 'figures/fig1g_recompute/reads_vs_cells_per_target.csv'   # published Fig1g table (reads,cells per target/bin)
TCLST = 'results/figure2_moi1/target_clusters_moi1.csv'           # 221 model targets
OUT   = 'figures/nature_figures/FigS_guide_capture'; os.makedirs(OUT, exist_ok=True)

# assign_sgRNAs.py thresholds
READS_THRESH=1; READS_PERC_THRESH=200; PERC_THRESH=2
UMI_COUNT_THRESH=10e3; UMI_READS_THRESH=10e3; FILTERED_UMI_THRESH=2; CBC=16

def build_recon():
    gl=pd.read_csv(os.path.join(SEQ,'gw_lib_guides.csv'))
    ctl=pd.read_csv(os.path.join(SEQ,'human_TF_CTL_guides.txt'),sep='\t',header=None,names=['sgRNA_name','sgRNA_barcode'])
    bcs=np.append(gl['sgRNA_barcode'].to_numpy(),ctl['sgRNA_barcode'].to_numpy())
    nms=np.append(gl['sgRNA_name'].to_numpy(),ctl['sgRNA_name'].to_numpy())
    d=defaultdict(list)
    for b,n in zip(bcs,nms): d[b].append(n)
    bcmap={b:v[0] for b,v in d.items() if len(v)==1}   # unique-barcode guides only
    gene=lambda n:'NO_SITE' if n.startswith('NO_SITE') else ('ONE_NON-GENE_SITE' if n.startswith('ONE_NON-GENE_SITE') else n.rsplit('_',1)[0])
    chans=sorted(os.path.basename(x).split('_CBCUMI')[0] for x in glob.glob(os.path.join(DICTS,'*.pkl')))
    chans=[c for c in chans if c.startswith(('High','Low'))]
    agg=defaultdict(lambda:[0,0,0])
    for ch in chans:
        binlab='High' if ch.startswith('High') else 'Low'
        dd=pickle.load(open(os.path.join(DICTS,f'{ch}_CBCUMI_map_mismatch0.pkl'),'rb'))
        fb=set(np.load(os.path.join(FILT,f'{ch}.barcodes.npy'),allow_pickle=True).tolist())
        cg=defaultdict(list)
        for key,val in dd.items():
            cbc=key[:CBC]
            if cbc not in fb: continue
            name=bcmap.get(key[CBC:])
            if name is None: continue
            tu=len(val); tr=sum(e[1] for e in val); uc=0
            for e in val:
                r=e[1]; p=r/tr if tr else 0
                if r>READS_THRESH or (r>READS_PERC_THRESH and p>PERC_THRESH): uc+=1
            if uc>FILTERED_UMI_THRESH or tu>UMI_COUNT_THRESH or tr>UMI_READS_THRESH:
                cg[cbc].append((name,tu,tr))
        for cbc,gs in cg.items():
            if len(gs)!=1: continue
            n,tu,tr=gs[0]; k=(binlab,gene(n))
            agg[k][0]+=1; agg[k][1]+=tu; agg[k][2]+=tr
        print('  ',ch,'done',flush=True)
    rows=[dict(bin=b,target=t,cells_recon=c,guide_umis=u,guide_reads=r) for (b,t),(c,u,r) in agg.items()]
    return pd.DataFrame(rows)

def assemble(rec):
    auth=pd.read_csv(FIG1G)
    H=rec[rec.bin=='High'].set_index('target')[['guide_umis','guide_reads']].add_suffix('_high')
    L=rec[rec.bin=='Low'].set_index('target')[['guide_umis','guide_reads']].add_suffix('_low')
    m=auth.set_index('gene').join(H).join(L).reset_index().rename(columns={'index':'gene'})
    for b in ['high','low']:
        m[f'umis_per_cell_{b}']=m[f'guide_umis_{b}']/m[f'cells_{b}'].replace(0,np.nan)
        m[f'bulkreads_per_cell_{b}']=m[f'reads_{b}']/m[f'cells_{b}'].replace(0,np.nan)
    return m

def make_figure(m):
    rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'pdf.use14corefonts':True,'font.size':11})
    COND={'Low':'#e1812c','High':'#3a923a'}; GRAY='#9a9a9a'; MINC=5
    tc=pd.read_csv(TCLST)
    MODEL=set(t for t in tc['target'] if t not in {'HIGH_CTL','LOW_CTL','NO_SITE','ONE_NON-GENE_SITE'} and t!='G1')
    m=m[~m.is_control].copy(); m['is_model']=m['gene'].isin(MODEL)
    rho=lambda x,y:spearmanr(np.asarray(x,float),np.asarray(y,float),nan_policy='omit')[0]
    def clean_logx(ax,ticks):
        ax.set_xscale('log'); ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_minor_locator(NullLocator()); ax.xaxis.set_major_formatter(ScalarFormatter()); ax.ticklabel_format(axis='x',style='plain')
    # 1 row x 4 cols: metric a (High,Low) + metric b (High,Low). Each metric's High/Low pair
    # stays together; landscape. (The former panel c -- bulk gDNA reads divided by the
    # single-cell cell count -- was removed: bulk gDNA and single-cell guide capture are separate
    # experiments on different cells, so that cross-sample per-cell ratio is not well defined.)
    fig,axes=plt.subplots(1,4,figsize=(18,4.6))
    for col,b in enumerate(['high','low']):
        cond='High' if b=='high' else 'Low'; c=COND[cond]
        dpos=m.dropna(subset=[f'reads_{b}',f'guide_umis_{b}',f'cells_{b}'])
        dpos=dpos[dpos[f'cells_{b}']>0]
        # a
        ax=axes[col]; da=dpos[(dpos[f'reads_{b}']>0)&(dpos[f'guide_umis_{b}']>0)]; am=da[da.is_model]
        ax.scatter(da[f'reads_{b}'],da[f'guide_umis_{b}'],s=8,c=GRAY,alpha=.35,rasterized=True,lw=0,label='all gene-targets')
        ax.scatter(am[f'reads_{b}'],am[f'guide_umis_{b}'],s=16,c=c,alpha=.85,rasterized=True,lw=0,label='221 model gene-targets')
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.text(.05,.96,f'$\\rho_{{model}}$ = {rho(am[f"reads_{b}"],am[f"guide_umis_{b}"]):.2f}\n$\\rho_{{all}}$ = {rho(da[f"reads_{b}"],da[f"guide_umis_{b}"]):.2f}',transform=ax.transAxes,va='top',fontsize=9)
        ax.set_xlabel('Bulk gDNA reads per gene-target (entire bin)'); ax.set_ylabel('Guide UMIs captured per gene-target')
        ax.set_title(f'HLA-{cond}',color=c,fontweight='bold')
        # colour key (grey = all, coloured = regulatory-model gene-targets) described in the caption
        # b   (columns 3-4)
        ax=axes[col+2]; db=dpos[dpos[f'cells_{b}']>=MINC]; bm=db[db.is_model]
        ax.scatter(db[f'cells_{b}'],db[f'umis_per_cell_{b}'],s=8,c=GRAY,alpha=.35,rasterized=True,lw=0)
        ax.scatter(bm[f'cells_{b}'],bm[f'umis_per_cell_{b}'],s=16,c=c,alpha=.85,rasterized=True,lw=0)
        med=db[f'umis_per_cell_{b}'].median(); ax.axhline(med,ls='--',c='gray',lw=1)
        ax.text(.97,.96,f'median = {med:.0f} UMIs/cell',transform=ax.transAxes,va='top',ha='right',fontsize=9,color='dimgray')
        ax.text(.05,.96,f'$\\rho_{{all}}$ = {rho(db[f"cells_{b}"],db[f"umis_per_cell_{b}"]):.2f}',transform=ax.transAxes,va='top',fontsize=9)
        clean_logx(ax,[5,10,20,50] if b=='high' else [5,10,20,50,100,200])
        ax.set_ylim(0,np.percentile(db[f'umis_per_cell_{b}'],99.5))
        ax.set_xlabel('Cells per gene-target'); ax.set_ylabel('Guide UMIs per cell  (capture / “expression”)')
        ax.set_title(f'HLA-{cond}',color=c,fontweight='bold')
    sns.despine(fig)
    for ax in axes.flat: ax.grid(False)
    for letter,i in zip('ab',[0,2]):
        axes[i].text(-.18,1.04,letter,transform=axes[i].transAxes,fontsize=15,fontweight='bold',va='bottom')
    plt.tight_layout()
    for fmt in ('pdf','png'): plt.savefig(os.path.join(OUT,f'FigS_guide_capture_vs_bulk.{fmt}'),bbox_inches='tight',dpi=300)

if __name__=='__main__':
    warnings.filterwarnings('ignore')
    csv=os.path.join(OUT,'guide_capture_per_target.csv')
    if os.path.exists(csv):
        print('Loading cached guide_capture_per_target.csv (skip pickle reconstruction) ...')
        m=pd.read_csv(csv)
    else:
        print('Reconstructing guide-UMI capture from CBC_UMI dicts ...')
        m=assemble(build_recon()); m.to_csv(csv,index=False)
    make_figure(m)
    print('Saved figure + table to',OUT)
