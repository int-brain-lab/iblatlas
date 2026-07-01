import logging
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter

import one.remote.aws as aws

from iblatlas import atlas
from iblatlas.genomics import agea

_logger = logging.getLogger(__name__)


def load(folder_cache=None):
    """
    Reads in the Allen gene expression experiments tables
    :param folder_cache:
    :return:
    df_cells: a dataframe of cells (8_879_868, 11), where each record corresponds to a single cell
    df_classes: a dataframe of classes (35, 3), where each record corresponds to a single class
    df_subclasses: a dataframe of subclasses (339, 4), where each record corresponds to a single subclass
    df_supertypes: a dataframe of supertypes (1202, 4), where each record corresponds to a single supertype
    df_clusters: a dataframe of clusters (5323, 5), where each record corresponds to a single cluster
    df_genes: a dataframe of genes (1672, 4), where each record corresponds to a single gene
    df_neurotransmitters: a dataframe of neurotransmitters (9, 2), where each record corresponds to a single
     neurotransmitter
    """
    OLD_VERSIONS = ['2024-02-02', '2023-06-12']
    folder_cache = Path(folder_cache or atlas.AllenAtlas._get_cache_dir().joinpath('merfish'))
    # check the AWS version and download the files if needed
    version_flag = next(folder_cache.glob('*.version'), None)
    if version_flag is None or version_flag.stem in OLD_VERSIONS:
        _logger.info(f'downloading gene expression data from {aws.S3_BUCKET_IBL} s3 bucket...')
        aws.s3_download_folder('atlas/merfish', folder_cache)
    # it is faster and more memory efficient to read the parquet files with dask, but we do
    # not want to require dask as a dependency so we provide the pandas alternative
    try:
        import dask.dataframe as dd
        df_cells = dd.read_parquet(list(folder_cache.rglob('*_cells.pqt')))
        df_cells = df_cells.compute()
    except Exception:  # there are more subtle errors than import errors if dask is intalled partially
        df_cells = pd.concat([pd.read_parquet(f) for f in folder_cache.rglob('*_cells.pqt')])
    # reads in the other tables
    df_classes = pd.read_parquet(folder_cache.joinpath('classes.pqt'))
    df_subclasses = pd.read_parquet(folder_cache.joinpath('subclasses.pqt'))
    df_supertypes = pd.read_parquet(folder_cache.joinpath('supertypes.pqt'))
    df_clusters = pd.read_parquet(folder_cache.joinpath('clusters.pqt'))
    df_genes = pd.read_parquet(folder_cache.joinpath('genes.pqt'))
    df_neurotransmitters = pd.read_parquet(folder_cache.joinpath('neurotransmitters.pqt'))
    return df_cells, df_classes, df_subclasses, df_supertypes, df_clusters, df_genes, df_neurotransmitters


def denoise_volume(volume, brain_mask, n_drop_non_neuronal=5, sigma=0.5, seed=42):
    """
    Denoise a raw MERFISH cell-type density volume.

    Raw volumes can contain NaN voxels with no nearby source data. Naively filling *every* NaN
    voxel with Dirichlet noise also fills the background outside the brain with random noise,
    which then bleeds a few voxels inward once Gaussian-smoothed. This restricts the
    fill/smooth/normalize steps to `brain_mask` and re-zeroes anything outside it after smoothing.

    :param volume: (n_types, dim0, dim1, dim2) raw density volume, as returned by
     `load_volume(..., label='')`
    :param brain_mask: (dim0, dim1, dim2) bool, True for voxels inside the brain (e.g.
     `atlas_agea.label != 0`)
    :param n_drop_non_neuronal: number of trailing non-neuronal types to drop (Astro-Epen,
     OPC-Oligo, OEC, Vascular, Immune are always the last 5 rows at the 'class' level); set to 0
     to keep all types
    :param sigma: Gaussian smoothing sigma in voxels, applied per type independently; set to 0 to
     disable
    :param seed: seed for the Dirichlet noise used to fill in-brain voxels that are NaN across
     every type
    :return: a (n_types - n_drop_non_neuronal, dim0, dim1, dim2) float32 array; in-brain voxels
     sum to ~1 over the type axis, out-of-brain voxels are 0
    """
    rng = np.random.RandomState(seed)
    vol = np.array(volume[: len(volume) - n_drop_non_neuronal], dtype=np.float32)
    all_nan = np.isnan(vol).all(axis=0) & brain_mask
    vol = np.nan_to_num(vol, nan=0.0)
    vol[:, all_nan] = rng.dirichlet(np.ones(vol.shape[0]), size=all_nan.sum()).T
    if sigma:
        vol = np.stack([gaussian_filter(channel, sigma=sigma) for channel in vol])
    vol *= brain_mask  # undo any smoothing bleed into the background
    vol /= vol.sum(axis=0, keepdims=True) + 1e-10
    return vol


def load_volume(level='class', label='processed', folder_cache=None):
    """
    Reads in a pre-computed MERFISH cell-type density volume and its type labels.

    These are dense 4-D arrays (one 3-D volume per cell type), built by
    `iblatlas/genomics/merfish_scrapping/02_create_volumes.py` on the same 200 um grid as
    `iblatlas.genomics.agea.load()` (not the default 25 um `AllenAtlas()` grid).

    :param level: taxonomy level to load, one of 'class', 'subclass', 'supertype', 'cluster'
    :param label: which volume to return
     - '': the raw, unprocessed volume (may contain NaN; returned memory-mapped)
     - 'processed': denoised via `denoise_volume()` (non-neuronal types dropped, NaNs filled,
       Gaussian-smoothed, renormalized to sum to 1 per in-brain voxel). Default -- unlike
       `agea.load()`, which defaults to the raw (`label=''`) volume.
    :param folder_cache:
    :return:
    volume: a (n_types, ml, dv, ap) array, one density volume per cell type, on the same grid as
     the `expression_volumes` returned by `agea.load()`. float16 memory-mapped if label='', float32
     in memory if label='processed'.
    labels: a (n_types,) array of type ids for each channel of `volume`, matching the index of the
     corresponding dataframe returned by `load()` (e.g. df_classes.index for level='class').
     Truncated to match `volume` when label='processed' drops non-neuronal types.
    atlas_agea: a brainatlas object with the labels and coordinates matching `volume` (same object
     as returned by `agea.load()` / `agea.load_atlas()`)
    """
    feasible_levels = ['class', 'subclass', 'supertype', 'cluster']
    assert level in feasible_levels, f'level must be one of {feasible_levels}'
    assert label in ('', 'processed'), "label must be '' (raw) or 'processed'"
    folder_cache = Path(folder_cache or atlas.AllenAtlas._get_cache_dir().joinpath('merfish'))
    # download only the 2 files needed for this level, not the whole `atlas/merfish` folder
    # (which also holds the multi-GB cell tables and the volumes for every other level)
    for filename in (f'merfish_{level}.npy', f'merfish_{level}_labels.npy'):
        file_path = folder_cache.joinpath(filename)
        if not file_path.exists():
            aws.s3_download_file(f'atlas/merfish/{filename}', file_path)
    volume = np.load(folder_cache.joinpath(f'merfish_{level}.npy'), mmap_mode='r')
    labels = np.load(folder_cache.joinpath(f'merfish_{level}_labels.npy'), allow_pickle=True)
    atlas_agea = agea.load_atlas()
    if label == 'processed':
        volume = denoise_volume(volume, atlas_agea.label != 0)
        labels = labels[:volume.shape[0]]
    return volume, labels, atlas_agea


def int2rgb(array, dtype=None):
    """
    One liner to convert rgba values stored as integer in dataframes
    :param array: rgba column of a dataframe or slice of the column
    :param dtype: optional, if int will return the uint8 view from 0-255 else will return floats from 0-1
    :return:
    """
    if dtype in (int, np.int8):
        return np.array(array).view('uint8').reshape(array.shape[0], 4)
    else:
        return np.array(array).view('uint8').reshape(array.shape[0], 4).astype(float) / 255
