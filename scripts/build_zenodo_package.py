"""
Assemble the Zenodo deposit for Perturb-ME (processed data, reviewer-accessible).

This is the record that carries the paper at submission time, while the ~850 GB of raw
FASTQ goes to GEO on its own (much slower) track. Everything here is processed or derived
data: no raw reads, so the whole package stays well inside Zenodo's default 50 GB quota
and needs no quota request.

Depends on scripts/build_deposit_h5ad.py having run first (it produces data/deposit/*.h5ad).

Run:
    LD_LIBRARY_PATH=/home/wangh256/miniforge3/envs/perturbme/lib \
    /home/wangh256/miniforge3/envs/perturbme/bin/python scripts/build_zenodo_package.py
"""
import hashlib
import os
import shutil
import sys

import numpy as np
import pandas as pd

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
OUT = os.path.join(REPO, 'submission/zenodo_package')

SUPP = os.path.join(REPO, 'doc/Supplementary Files')
MAGECK = os.path.join(SUPP, 'Supplementary Table 4 - MAGeCK bulk screen')
FIG2 = os.path.join(REPO, 'figures/nature_figures/Fig2')
MODEL = os.path.join(REPO, 'doc/co-scientist/gene_program_discovery')
DEPOSIT = os.path.join(REPO, 'data/deposit')

# (source, destination) pairs copied verbatim.
COPY = [
    # --- single-cell objects ---
    (os.path.join(DEPOSIT, 'PerturbME_A375_HLA_all_cells.h5ad'),
     '01_single_cell/PerturbME_A375_HLA_all_cells.h5ad'),
    (os.path.join(DEPOSIT, 'PerturbME_A375_HLA_MOI1.h5ad'),
     '01_single_cell/PerturbME_A375_HLA_MOI1.h5ad'),

    # --- regulatory model: canonical partition only (see
    #     doc/notes/20260612_Fig2_module_genelist_provenance.md for why the archived
    #     module runs must NOT be used) ---
    # target_modules.csv is NOT copied verbatim -- see annotate_target_modules()
    (os.path.join(FIG2, 'Fig2B_gene_programs.csv'),
     '02_regulatory_model/gene_programs.csv'),
    (os.path.join(FIG2, 'Fig2C_module_program_edges.csv'),
     '02_regulatory_model/module_program_edges.csv'),
    (os.path.join(FIG2, 'Fig2C_target_module_top8.csv'),
     '02_regulatory_model/target_module_representatives.csv'),
    (os.path.join(FIG2, 'Fig2C_gene_program_top8.csv'),
     '02_regulatory_model/gene_program_representatives.csv'),
    (os.path.join(FIG2, 'Fig2d_gene_reconciliation.csv'),
     '02_regulatory_model/Fig2d_curated_gene_reconciliation.csv'),

    # --- bulk CRISPR screen ---
    (os.path.join(MAGECK, 'all_samples.count.txt'),
     '03_bulk_screen/mageck_sgRNA_counts.txt'),
    (os.path.join(MAGECK, 'all_samples.count_normalized.txt'),
     '03_bulk_screen/mageck_sgRNA_counts_normalized.txt'),

    # --- agentic interpretation ---
    (os.path.join(SUPP, 'Supplementary File 1 - CoScientist Prompt.docx'),
     '04_cosci/CoScientist_prompt.docx'),
    (os.path.join(SUPP, 'Supplementary File 2 - CoScientist Response.pdf'),
     '04_cosci/CoScientist_response.pdf'),
    # NB: the former "Supplementary Code 1" is deliberately NOT shipped. It was a copy of
    # notebook/20260220_ElasticNet_Beta_Module_Analysis.ipynb, a superseded run whose feature
    # filter yields 2,579 features and a target partition agreeing with the published one for
    # only 193/221 targets. The full codebase is released instead, where
    # notebook/20260519_Figure2B_Rerun_BetaCorr.ipynb reproduces the published partition
    # exactly (verified: 224/224 targets, identical labels).

    # --- reference ---
    (os.path.join(REPO, 'doc/Supplementary Files/TableS1_CITE_antibody_panel.csv'),
     '05_reference/CITE_antibody_panel.csv'),
]

MAGECK_COMPARISONS = ['Low_vs_Input', 'High_vs_Input', 'CTL_vs_Input',
                      'Low_vs_CTL', 'High_vs_CTL', 'High_vs_Low']


def annotate_target_modules(dest_dir):
    """Ship the module partition with an explicit covariate_type column.

    Fig2B_target_modules.csv has 224 rows, not the 221 gene-target perturbations the
    manuscript describes: K-means was run over the beta columns, which also include the two
    pooled control categories (HIGH_CTL, LOW_CTL) and -- less expectedly -- the G1 cell-cycle
    covariate. Those three are genuinely part of what was clustered, so the partition is not
    altered here; instead each row is labelled so a reviewer can see at a glance that exactly
    221 rows are gene targets.
    """
    df = pd.read_csv(os.path.join(FIG2, 'Fig2B_target_modules.csv'))
    controls = {'HIGH_CTL', 'LOW_CTL'}
    technical = {'G1', 'G2M', 'S', 'UMI', 'cells_per_target'}
    df['covariate_type'] = np.where(
        df['target'].isin(controls), 'pooled_control',
        np.where(df['target'].isin(technical) | df['target'].str.startswith('TENX_'),
                 'technical_covariate', 'gene_target'))
    counts = df['covariate_type'].value_counts().to_dict()
    print(f'  target modules: {counts}', flush=True)
    assert counts.get('gene_target') == 221, \
        f"expected 221 gene targets, got {counts.get('gene_target')}"
    df.to_csv(os.path.join(dest_dir, 'target_modules.csv'), index=False)


def write_beta_matrix(dest_dir):
    """EM-refined ElasticNet coefficients, 18,683 features x 231 covariates."""
    B = np.load(os.path.join(MODEL, 'EN_B_EM.npy'))
    features = np.load(os.path.join(MODEL, 'feature_names.npy'), allow_pickle=True).astype(str)
    covs = np.load(os.path.join(MODEL, 'cov_names.npy'), allow_pickle=True).astype(str)
    assert B.shape == (features.size, covs.size), f'shape mismatch {B.shape}'
    df = pd.DataFrame(B, index=pd.Index(features, name='feature'), columns=covs)
    out = os.path.join(dest_dir, 'elasticnet_beta_matrix.csv.gz')
    df.to_csv(out, float_format='%.6g', compression='gzip')
    print(f'  beta matrix {df.shape} -> {os.path.basename(out)} '
          f'({os.path.getsize(out) / 2**20:.1f} MB)', flush=True)


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def main():
    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)

    missing = []
    for src, rel in COPY:
        dest = os.path.join(OUT, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if not os.path.exists(src):
            missing.append(src)
            print(f'  MISSING: {src}', file=sys.stderr, flush=True)
            continue
        shutil.copy2(src, dest)
        print(f'  {rel}  ({os.path.getsize(dest) / 2**20:.1f} MB)', flush=True)

    # MAGeCK per-comparison summaries
    for comp in MAGECK_COMPARISONS:
        for kind in ('gene_summary', 'sgrna_summary'):
            src = os.path.join(MAGECK, f'{comp}.{kind}.txt')
            if os.path.exists(src):
                dest = os.path.join(OUT, f'03_bulk_screen/mageck_{comp}.{kind}.txt')
                shutil.copy2(src, dest)

    shutil.copy2(os.path.join(REPO, 'scripts/zenodo_README_template.md'),
                 os.path.join(OUT, 'README.md'))

    annotate_target_modules(os.path.join(OUT, '02_regulatory_model'))
    write_beta_matrix(os.path.join(OUT, '02_regulatory_model'))

    # checksums
    rows = []
    for root, _, files in os.walk(OUT):
        for name in sorted(files):
            p = os.path.join(root, name)
            rows.append({'file': os.path.relpath(p, OUT),
                         'bytes': os.path.getsize(p),
                         'md5': md5(p)})
    df = pd.DataFrame(rows).sort_values('file')
    df.to_csv(os.path.join(OUT, 'MD5SUMS.csv'), index=False)

    total = df['bytes'].sum()
    print(f'\n{len(df)} files, {total / 2**30:.2f} GB total '
          f'({"OK" if total < 50 * 2**30 else "OVER"} vs Zenodo 50 GB default quota)')
    if missing:
        print(f'\n{len(missing)} MISSING inputs -- package incomplete', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
