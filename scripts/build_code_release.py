"""
Assemble the public code release for Perturb-ME (GitHub -> Zenodo archive -> DOI).

The working tree cannot be published as-is:
  - it is not a git repository in any useful sense (.gitignore is a bare `*`, zero commits)
  - src/ carries ~198 MB of permutation-test .npy outputs that are results, not code
  - three notebooks exceed or approach GitHub's 100 MB per-file hard limit purely because of
    embedded outputs (the largest is 162 MB and would be rejected outright)

So this script selects code only, strips notebook outputs, and generates the scaffolding a
public repository needs. Actual code is ~0.8 MB across 116 files.

Run:
    LD_LIBRARY_PATH=/home/wangh256/miniforge3/envs/perturbme/lib \
    /home/wangh256/miniforge3/envs/perturbme/bin/python scripts/build_code_release.py
"""
import json
import os
import re
import shutil
import sys

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
OUT = os.path.join(REPO, 'submission/code_release')

CODE_EXT = {'.py', '.R', '.sh'}
SKIP_DIRS = {'__pycache__', '.ipynb_checkpoints', '.git', 'archive'}
# build_supplementary_files.py only renumbers the submission's supplementary items; it is
# submission logistics rather than analysis code, so it stays out of the public release.
SKIP_NAMES = {'.DS_Store', 'build_supplementary_files.py'}

# --- what counts as "part of the paper" -------------------------------------------------
# The working tree carries several exploratory branches that the final manuscript does not
# use: grepping the Methods for harmony / ComBat / totalVI / scVI / permutation / NMF /
# Seurat returns nothing. Those are excluded, with the reason recorded so the decision can
# be revisited.
INCLUDE_TREES = [
    'src/RNA',                  # cell/gene filtering, Terra transfer, QC
    'src/Hash',                 # hashtag demultiplexing (Cumulus)
    'src/CITE',                 # antibody-derived tag normalisation
    'src/gDNA',                 # bulk sgRNA enrichment / MAGeCK inputs
    'src/manuscript_figures',   # one module per published panel + submission tooling
    'scripts',
]
INCLUDE_FILES = [
    'src/output_adata.py',
    'src/output_sgRNA_combinations.py',
    'src/clustering/cluster_conditions.py',        # Fig 1c
    'src/clustering/cluster_raw_expression.py',    # expression heatmap
    'src/clustering/cluster_cell_cycle.py',        # Tirosh phase -> model covariate
    'src/CROP/linear_model/lm_target_EM.py',                        # the published model
    'src/CROP/linear_model/cluster_by_kmeans_genome_wide_features.py',  # modules + programs
    'src/CROP/linear_model/cluster_by_kmeans.py',
]
# src/CROP/*.py at the top level (guide calling and coverage) is included wholesale;
# src/CROP/linear_model/ is not, except for the three files listed above.
INCLUDE_CROP_TOPLEVEL = True

EXCLUDED_WITH_REASON = {
    'src/batch_correction': 'Harmony / scVI / ComBat - no batch correction in the final Methods',
    'src/DE': 'Seurat and scVI differential expression - not used; the model uses all 18,683 features',
    'src/Programs': 'HLA program scoring - feeds the archived Fig 2d composite-score panel',
    'src/CROP/linear_model (variants)': 'totalVI / scVI / NMF / interactions / preEM / '
                                        'all_features / grouped_threshold / factorize_recover / '
                                        'permutation_test - alternative fits not reported',
    'src/output_adata_scVI.py, src/run_totalVI.py': 'scVI / totalVI branches not used',
}

# Notebooks that write into figures/nature_figures (i.e. produced published panels),
# plus the one shipped with the paper as Supplementary Code 1.
INCLUDE_NOTEBOOKS = [
    '20260417_Figure1_Nature.ipynb',
    '20260417_Figure2_Nature.ipynb',
    '20260519_Figure2B_Rerun_BetaCorr.ipynb',   # canonical Fig 2b/2c partition
    '20260513_FigS1_Supplementary.ipynb',
    'Fig_2g_regulatory_map.ipynb',
]
# Deliberately NOT released: 20260220_ElasticNet_Beta_Module_Analysis.ipynb. It was shipped
# with the submission as "Supplementary Code 1", but it is the superseded February run --
# its feature filter yields 2,579 features (the paper reports 1,998) and its target partition
# agrees with the published one for only 193/221 targets. Supplementary Code 1 has been
# withdrawn in favour of this full codebase; including the notebook here would reintroduce
# exactly the contradiction that removing it resolves.

# Absolute paths that only exist on the institutional cluster. Not secrets, but they make the
# code unrunnable elsewhere and leak internal layout, so the README must call them out.
PATH_PATTERNS = ['/gstore', '/gnet', '/ahg/regevdata', '/home/wangh256']


def strip_notebook(src, dest):
    """Write the notebook with all outputs cleared. Keeps source, drops embedded images."""
    with open(src, encoding='utf-8') as f:
        nb = json.load(f)
    for cell in nb.get('cells', []):
        if cell.get('cell_type') == 'code':
            cell['outputs'] = []
            cell['execution_count'] = None
        cell.get('metadata', {}).pop('execution', None)
    nb.get('metadata', {}).pop('widgets', None)
    with open(dest, 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
        f.write('\n')


def copy_code_tree(rel_root):
    """Copy only source files under rel_root, preserving layout."""
    n, total = 0, 0
    for root, dirs, files in os.walk(os.path.join(REPO, rel_root)):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for name in sorted(files):
            if name in SKIP_NAMES or os.path.splitext(name)[1] not in CODE_EXT:
                continue
            src = os.path.join(root, name)
            dest = os.path.join(OUT, os.path.relpath(src, REPO))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copy2(src, dest)
            n += 1
            total += os.path.getsize(dest)
    print(f'  {rel_root}: {n} files, {total / 2**20:.2f} MB', flush=True)
    return n


def copy_single(rel):
    src = os.path.join(REPO, rel)
    if not os.path.exists(src):
        print(f'  MISSING: {rel}', file=sys.stderr)
        return 0
    dest = os.path.join(OUT, rel)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    shutil.copy2(src, dest)
    return 1


def copy_crop_toplevel():
    src_dir = os.path.join(REPO, 'src/CROP')
    n = 0
    for name in sorted(os.listdir(src_dir)):
        if os.path.splitext(name)[1] in CODE_EXT and name not in SKIP_NAMES:
            n += copy_single(f'src/CROP/{name}')
    print(f'  src/CROP (top level): {n} files', flush=True)
    return n


def copy_notebooks():
    src_dir = os.path.join(REPO, 'notebook')
    dest_dir = os.path.join(OUT, 'notebook')
    os.makedirs(dest_dir, exist_ok=True)
    n, before, after = 0, 0, 0
    for name in sorted(INCLUDE_NOTEBOOKS):
        if not os.path.exists(os.path.join(src_dir, name)):
            print(f'  MISSING notebook: {name}', file=sys.stderr)
            continue
        src = os.path.join(src_dir, name)
        dest = os.path.join(dest_dir, name)
        before += os.path.getsize(src)
        try:
            strip_notebook(src, dest)
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            print(f'  WARNING: could not strip {name} ({e}); skipped', file=sys.stderr)
            continue
        after += os.path.getsize(dest)
        n += 1
    print(f'  notebook: {n} notebooks, {before / 2**20:.0f} MB -> {after / 2**20:.1f} MB '
          f'(outputs stripped)', flush=True)
    return n


def audit_paths():
    """Report which released files still carry cluster-absolute paths."""
    hits = {p: [] for p in PATH_PATTERNS}
    for root, dirs, files in os.walk(OUT):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for name in files:
            if os.path.splitext(name)[1] not in CODE_EXT | {'.ipynb'}:
                continue
            p = os.path.join(root, name)
            try:
                text = open(p, encoding='utf-8', errors='ignore').read()
            except OSError:
                continue
            for pat in PATH_PATTERNS:
                if pat in text:
                    hits[pat].append(os.path.relpath(p, OUT))
    print('\n  path audit (files still containing cluster-absolute paths):', flush=True)
    for pat, files in hits.items():
        print(f'    {pat:20s} {len(files)} files', flush=True)
    return hits


GITIGNORE = """\
# data and derived objects -- never in the code repository
data/
results/
figures/
submission/
PerturbME_transfer/
*.h5ad
*.h5
*.npy
*.npz
*.pkl
*.bam
*.fastq.gz

# notebook outputs are stripped before release
.ipynb_checkpoints/

# python
__pycache__/
*.py[cod]
.venv/
env/

# os / editor
.DS_Store
.vscode/
"""


MIT_LICENSE = """\
MIT License

Copyright (c) 2026 {holder}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

# CONFIRM BEFORE PUBLISHING: most authors are at Genentech and correspondence is @gene.com,
# so the institution is the likely copyright holder rather than the individual authors.
COPYRIGHT_HOLDER = 'Genentech, Inc.'


def write_scaffolding(counts):
    with open(os.path.join(OUT, '.gitignore'), 'w') as f:
        f.write(GITIGNORE)
    with open(os.path.join(OUT, 'LICENSE'), 'w') as f:
        f.write(MIT_LICENSE.format(holder=COPYRIGHT_HOLDER))

    tmpl = os.path.join(REPO, 'scripts/code_release_README_template.md')
    if os.path.exists(tmpl):
        text = open(tmpl, encoding='utf-8').read()
        text = text.replace('{{N_PY}}', str(counts['src'] + counts['scripts']))
        text = text.replace('{{N_NB}}', str(counts['notebook']))
        with open(os.path.join(OUT, 'README.md'), 'w', encoding='utf-8') as f:
            f.write(text)
    else:
        print(f'  WARNING: no README template at {tmpl}', file=sys.stderr)

    setup = os.path.join(REPO, 'setup.sh')
    if os.path.exists(setup):
        shutil.copy2(setup, os.path.join(OUT, 'setup.sh'))


def main():
    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)

    n_code = sum(copy_code_tree(t) for t in INCLUDE_TREES)
    if INCLUDE_CROP_TOPLEVEL:
        n_code += copy_crop_toplevel()
    n_single = sum(copy_single(f) for f in INCLUDE_FILES)
    print(f'  individually selected: {n_single} files', flush=True)

    counts = {'src': n_code + n_single, 'scripts': 0}
    counts['notebook'] = copy_notebooks()

    print('\n  excluded as not part of the final paper:', flush=True)
    for what, why in EXCLUDED_WITH_REASON.items():
        print(f'    {what}\n        {why}', flush=True)

    write_scaffolding(counts)
    audit_paths()

    total = sum(os.path.getsize(os.path.join(r, f))
                for r, _, fs in os.walk(OUT) for f in fs)
    n = sum(len(fs) for _, _, fs in os.walk(OUT))
    print(f'\n  {n} files, {total / 2**20:.1f} MB total', flush=True)
    biggest = max((os.path.getsize(os.path.join(r, f)), os.path.join(r, f))
                  for r, _, fs in os.walk(OUT) for f in fs)
    print(f'  largest file: {biggest[0] / 2**20:.1f} MB  '
          f'{os.path.relpath(biggest[1], OUT)} '
          f'({"OK" if biggest[0] < 100 * 2**20 else "OVER GitHub 100 MB limit"})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
