"""Loads the lipid volumes derived from the Lipid Brain Atlas (Fusar Bassini et al. 2025)

MALDI mass spectrometry imaging of 173 lipids registered to the Allen CCF, interpolated on the
200 um AGEA / MERFISH grid. Generation scripts are in the `lipids_scrapping` folder.

The volumes are a derivative of the Lipid Brain Atlas data, distributed under CC-BY 4.0: credit
Fusar Bassini et al., Nature (2026), doi:10.1038/s41586-026-11050-0 and doi:10.5281/zenodo.15379565
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from one.remote import aws

from iblatlas import atlas
from iblatlas.genomics import agea

_logger = logging.getLogger(__name__)


def load_volume(folder_cache=None, min_reliability=None):
    """
    Reads in the pre-computed lipid volumes and the lipids table.

    The volumes are built by `iblatlas/genomics/lipids_scrapping/02_create_volumes.py` from the 2
    densely sampled naive male brains of the Lipid Brain Atlas, on the same 200 um grid as
    `iblatlas.genomics.agea.load()`.

    Parameters
    ----------
    folder_cache : str or Path, optional
        Folder containing `lipid_volumes.npy` and `lipids.pqt`, defaults to the `lipids` folder of the
        iblatlas cache. Missing files are downloaded from the IBL public S3 bucket (`atlas/lipids`).
    min_reliability : float, optional
        If set, keeps only the lipids whose `reliability` (voxel-wise correlation between the volumes
        of the 2 brains) is at least this value. The volume is then loaded in memory.

    Returns
    -------
    volume : np.ndarray
        (n_lipids, ml, dv, ap) float32 log intensities (uMAIA-normalised arbitrary units), NaN outside
        the brain and outside the sampled range (olfactory bulb and anterior frontal cortex).
        Memory-mapped unless `min_reliability` is set.
    df_lipids : pd.DataFrame
        (n_lipids, 4) one row per volume channel, columns:
        - lipid: lipid name
        - reliability: voxel-wise correlation between the volumes of the 2 dense brains
        - sex_log_ratio: median over Beryl regions of the female - male log intensity
        - sex_zscore: median over Beryl regions of the female - male difference over the pooled
          between-brain standard deviation (3 vs. 3 sparsely sampled brains, indicative only)
    atlas_agea : iblatlas.atlas.BrainAtlas
        Brain atlas with the labels and coordinates matching `volume` (same as `agea.load_atlas()`)
    """
    folder_cache = Path(folder_cache or atlas.AllenAtlas._get_cache_dir().joinpath('lipids'))
    for filename in ('lipid_volumes.npy', 'lipids.pqt'):
        file_path = folder_cache.joinpath(filename)
        if not file_path.exists():
            _logger.info(f'downloading {filename} from {aws.S3_BUCKET_IBL} s3 bucket...')
            aws.s3_download_file(f'atlas/lipids/{filename}', file_path)
    volume = np.load(folder_cache.joinpath('lipid_volumes.npy'), mmap_mode='r')
    df_lipids = pd.read_parquet(folder_cache.joinpath('lipids.pqt'))
    if min_reliability is not None:
        keep = df_lipids['reliability'].values >= min_reliability
        volume, df_lipids = volume[keep], df_lipids.loc[keep].reset_index(drop=True)
    return volume, df_lipids, agea.load_atlas()
