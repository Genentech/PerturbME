"""
Render Supplementary Code 1 (a Jupyter notebook) to PDF.

.ipynb is not one of the formats RESIS accepts; PDF is. There is no LaTeX or
LibreOffice on this box, so the route is nbconvert -> HTML -> headless Chrome.

Run:  conda run -n perturbme python src/manuscript_figures/nb_to_pdf.py
"""
import os
import subprocess
import sys

import nbformat
from nbconvert import HTMLExporter

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
NB = os.path.join(REPO, 'doc/Supplementary Files',
                  'Supplementary Code 1 - ElasticNet beta module analysis.ipynb')
BUILD = os.path.join(REPO, 'submission/build')
HTML = os.path.join(BUILD, 'SupplementaryCode1.html')
PDF = os.path.join(BUILD, 'Supplementary Code 1 - ElasticNet beta module analysis.pdf')

# Chrome's print path clips long source lines at the page edge instead of wrapping,
# so force wrapping and keep output images inside the page box.
PRINT_CSS = """
<style>
@page { size: A4; margin: 12mm; }
body { font-size: 11px; }
pre, code,
.highlight pre,
.jp-RenderedText pre,
.jp-OutputArea-output pre,
.jp-InputArea-editor {
  white-space: pre-wrap !important;
  overflow-wrap: anywhere !important;
  word-break: break-word !important;
  overflow: visible !important;
}
.jp-Cell, .jp-OutputArea-child { break-inside: avoid; }
img, svg { max-width: 100% !important; height: auto !important; }
</style>
</head>"""

os.makedirs(BUILD, exist_ok=True)
body, _ = HTMLExporter(template_name='lab').from_notebook_node(nbformat.read(NB, as_version=4))
open(HTML, 'w').write(body.replace('</head>', PRINT_CSS, 1))

subprocess.run(['google-chrome', '--headless=new', '--disable-gpu', '--no-sandbox',
                '--no-pdf-header-footer', '--run-all-compositor-stages-before-draw',
                '--virtual-time-budget=30000', f'--print-to-pdf={PDF}', f'file://{HTML}'],
               check=True, capture_output=True, timeout=600)
os.remove(HTML)

if not os.path.exists(PDF):
    sys.exit(f'failed to write {PDF}')
print(f'{PDF}  ({os.path.getsize(PDF) / 1e6:.1f} MB)')
