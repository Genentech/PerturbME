"""
Build lean, deposit-ready AnnData objects for GEO / Zenodo / UCSC Cell Browser.

The working h5ad files are far too large to deposit:

    data/20251116_processed.h5ad       46 GB  (353,091 cells, all MOI)
    data/20251116_processed_1moi.h5ad  36 GB  (277,152 cells, MOI = 1)

Almost all of that is redundant:
  - layers/prev_counts is a DENSE (n_cells x 18,683) float array -> ~21-26 GB on its own
  - X and layers/counts hold byte-identical raw counts -> one of the two is pure duplication
  - obsp/{connectivities,distances} is the kNN graph, fully re-derivable from X_pca

This script keeps counts once (as X, float32, gzip), carries every annotation across, and
adds the two obs columns that the Methods describe but the working files never stored:

  cell_cycle_phase  Tirosh et al. G1/S/G2M call, from adata_arrs/cell_cycle_Tirosh.npy
  sequencing_run    the "TENX" model covariate, re-derived from the Sample name using the
                    exact rule in src/CROP/linear_model/lm_target_EM.py:224-232 (assign 2
                    first so that a sample containing both '1' and '2' lands in run 1)

Run (LD_LIBRARY_PATH is required - the conda pandas needs the env's libstdc++, same trap as
src/manuscript_figures/ExtFig_FACS_gating.py):

    LD_LIBRARY_PATH=/home/wangh256/miniforge3/envs/perturbme/lib \
    /home/wangh256/miniforge3/envs/perturbme/bin/python scripts/build_deposit_h5ad.py
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd

try:
    from anndata.io import read_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem
import anndata as ad

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
CELL_CYCLE = os.path.join(
    REPO, 'PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/adata_arrs/cell_cycle_Tirosh.npy')

FULL = os.path.join(REPO, 'data/20251116_processed.h5ad')
MOI1 = os.path.join(REPO, 'data/20251116_processed_1moi.h5ad')
OUTDIR = os.path.join(REPO, 'data/deposit')

JOBS = [
    (FULL, 'PerturbME_A375_HLA_all_cells.h5ad',
     'all 353,091 QC-passing cells, all sgRNA multiplicities'),
    (MOI1, 'PerturbME_A375_HLA_MOI1.h5ad',
     '277,152 sgRNA-singlet cells used for the regulatory model'),
]


def sequencing_run(sample_names):
    """Re-derive the 10x run covariate (TENX_1/2/3) from the sample name.

    Mirrors src/CROP/linear_model/lm_target_EM.py:224-232 exactly, including the
    assignment order: 2 first, then 1, then 3, so that a sample whose name contains
    both digits (e.g. Low_A21) is assigned to run 1.
    """
    s = np.asarray(sample_names).astype(str)
    run = np.full(s.size, '0', dtype='<U1')
    run[np.char.find(s, '2') != -1] = '2'
    run[np.char.find(s, '1') != -1] = '1'
    run[np.char.find(s, '3') != -1] = '3'
    assert not (run == '0').any(), 'sample name matched no 10x run digit'
    return run


def load_cell_cycle(obs_index, full_obs_index):
    """Tirosh phase calls, stored in the order of the 353K full object."""
    cc = np.load(CELL_CYCLE, allow_pickle=True).astype(str)
    assert cc.size == full_obs_index.size, (
        f'cell-cycle array is {cc.size} long but the full object has {full_obs_index.size} cells')
    if obs_index.size == full_obs_index.size:
        assert (obs_index == full_obs_index).all()
        return cc
    # subset object: align on barcode
    pos = pd.Index(full_obs_index).get_indexer(pd.Index(obs_index))
    assert (pos >= 0).all(), 'some cells are absent from the full object; cannot align'
    return cc[pos]


def build(src, out_name, description, full_obs_index):
    print(f'\n=== {os.path.basename(src)} -> {out_name} ===', flush=True)
    with h5py.File(src, 'r') as f:
        print('  reading counts ...', flush=True)
        X = read_elem(f['layers/counts'])
        obs = read_elem(f['obs'])
        var = read_elem(f['var'])
        obsm = {k: read_elem(f['obsm'][k]) for k in ('X_pca', 'X_umap') if k in f['obsm']}
        # Authoritative description of what X actually holds. NB the RNA and CITE blocks are
        # NOT on the same scale: RNA is raw integer UMI, CITE is already IgG-normalised.
        layer_info = {k: f['uns/counts_layer_info'][k][()].decode()
                      for k in f['uns/counts_layer_info']}

        # X and layers/counts are byte-identical in these files; verify before dropping one.
        n = min(1_000_000, f['X/data'].shape[0])
        assert np.array_equal(f['X/data'][:n], f['layers/counts/data'][:n]), \
            'X and layers/counts differ - do not drop X blindly'

    if X.dtype != np.float32:
        print(f'  downcasting counts {X.dtype} -> float32 (lossless for integer counts)',
              flush=True)
        X = X.astype(np.float32)

    obs = obs.copy()
    obs['cell_cycle_phase'] = pd.Categorical(load_cell_cycle(obs.index.values, full_obs_index))
    obs['sequencing_run'] = pd.Categorical(sequencing_run(obs['Sample'].values))

    var = var.copy()
    is_cite = var.index.str.startswith('CITE-')
    var['modality'] = pd.Categorical(np.where(is_cite, 'protein', 'RNA'))
    print(f'  features: {(~is_cite).sum()} RNA + {is_cite.sum()} CITE-seq protein', flush=True)

    adata = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm)
    adata.uns['title'] = 'Perturb-ME: genome-scale Perturb-CITE-seq of MHC-I regulation in A375'
    adata.uns['description'] = description
    adata.uns['X_contents'] = (
        'X is NOT on a single scale. The 18,670 RNA features (var.modality == "RNA") hold '
        f'{layer_info["rna_type"]}. The 13 CITE-seq features (var.modality == "protein", '
        f'prefixed "CITE-") hold {layer_info["cite_type"]}, i.e. they are already normalised '
        'and must not be re-normalised with the RNA counts.')
    adata.uns['counts_layer_info'] = layer_info

    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, out_name)
    print(f'  writing {out} ...', flush=True)
    adata.write_h5ad(out, compression='gzip', compression_opts=4)
    print(f'  done: {os.path.getsize(out) / 2**30:.2f} GB '
          f'(source {os.path.getsize(src) / 2**30:.1f} GB)', flush=True)
    print(adata, flush=True)
    return out


def main():
    with h5py.File(FULL, 'r') as f:
        full_obs_index = read_elem(f['obs/index']).astype(str)
    print(f'full object: {full_obs_index.size} cells; '
          f'barcodes unique: {pd.Index(full_obs_index).is_unique}', flush=True)

    for src, out_name, description in JOBS:
        if not os.path.exists(src):
            print(f'SKIP (missing): {src}', file=sys.stderr)
            continue
        build(src, out_name, description, full_obs_index)


if __name__ == '__main__':
    main()
