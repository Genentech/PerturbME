"""
Stage 1 of the Resis / editorial "incomplete submission" package.

Renders every manuscript figure to a high-resolution raster image (600 dpi TIFF-LZW
plus a 300 dpi PNG preview) and copies the vector originals alongside.

Resis will not accept figures embedded in a .docx/.pptx, so each figure has to ship
as its own image file.

Run:  conda run -n perturbme python src/manuscript_figures/build_submission_figures.py
"""
import os
import shutil
import subprocess

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
PKG = os.path.join(REPO, 'submission/20260721_Resis_package')
FIGDIR = os.path.join(PKG, '02_Figures')
VECDIR = os.path.join(FIGDIR, 'vector_originals')
NF = os.path.join(REPO, 'figures/nature_figures')

DPI_TIFF = 600
DPI_PNG = 300

# (source pdf, output basename)
# doc/Figure{1,2}.pdf are the current layouts: Fig 1 = a-g, Fig 2 = a-e.
# figures/final/Figure{1,2}.pdf are an older layout (Fig 2 still carries the
# expression heatmap that has since been archived) - not shipped.
MAIN = [
    (os.path.join(REPO, 'doc/Figure1.pdf'), 'Figure_1'),
    (os.path.join(REPO, 'doc/Figure2.pdf'), 'Figure_2'),
]

# Extended Data figures, per the captions in doc/manuscript.pdf:
#   ED Fig 1 = the CITE-seq marker violin row (formerly FigS1 panel c)
#   ED Fig 2 = per-target bulk reads vs cells vs guide UMIs
# The earlier FigS1 a i-iii panels (cells per channel, per-sample HLA ridgeline, IFNg
# signature) are no longer in the paper and are not shipped.
#   ED Fig 3 = FACS gating strategy (new, from src/manuscript_figures/ExtFig_FACS_gating.py)
# cells_per_guide_gene is a SEPARATE ED figure (coverage-depth distributions) - deliberately
# kept distinct from the guide-capture concordance figure above: different point, different
# structure, and cited in a different Results paragraph (the coverage sentence). ED number TBD.
PANELS = [
    (os.path.join(NF, 'FigS1/S1c_violin_row.pdf'), 'ExtendedData_Fig1_CITE_marker_violins'),
    (os.path.join(NF, 'FigS_guide_capture/FigS_guide_capture_vs_bulk.pdf'), 'ExtendedData_Fig2_guide_capture_vs_bulk'),
    (os.path.join(NF, 'FigS_FACS_gating/ExtFig_FACS_gating.pdf'), 'ExtendedData_Fig3_FACS_gating'),
    (os.path.join(NF, 'FigS_cells_per_guide_gene/ExtFig_cells_per_guide_gene.pdf'), 'ExtendedData_cells_per_guide_gene'),
]

log = []


def render(src, base, outdir):
    if not os.path.exists(src):
        log.append(('MISSING', src, ''))
        return
    tif = os.path.join(outdir, base)
    png = os.path.join(outdir, base)
    subprocess.run(['pdftoppm', '-tiff', '-tiffcompression', 'lzw', '-r', str(DPI_TIFF),
                    '-singlefile', src, tif], check=True)
    subprocess.run(['pdftoppm', '-png', '-r', str(DPI_PNG),
                    '-singlefile', src, png], check=True)
    shutil.copy2(src, os.path.join(VECDIR, base + '.pdf'))
    for ext in ('.tif', '.png'):
        f = tif + ext
        log.append(('OK', os.path.relpath(f, PKG), f'{os.path.getsize(f)/1e6:.1f} MB'))


if os.path.isdir(FIGDIR):
    shutil.rmtree(FIGDIR)
os.makedirs(VECDIR, exist_ok=True)
for src, base in MAIN:
    render(src, base, FIGDIR)

PANELDIR = os.path.join(FIGDIR, 'ExtendedData_panels')
os.makedirs(PANELDIR, exist_ok=True)
for src, base in PANELS:
    render(src, base, PANELDIR)

report = ['stage 1: figure rendering', '']
report += [f'  {s:8s} {p}  {sz}' for s, p, sz in log]
print('\n'.join(report))
open(os.path.join(PKG, '.stage1.log'), 'w').write('\n'.join(report))
