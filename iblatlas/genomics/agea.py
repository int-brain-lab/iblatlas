"""A package for loading 4345 formatted and registered gene expression volumes
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import one.remote.aws as aws

from iblatlas import atlas

_logger = logging.getLogger(__name__)

_, NML, NDV, NAP = DIM_EXP = (4345, 58, 41, 67)  # nexperiments, nml, ndv, nap


def load(folder_cache=None, expression_size=DIM_EXP):
    """
    Reads in the Allen gene expression experiments binary data.
    Generation scripts from the Allen Institute are available in the gene-expression-scrapping folder
    :param filename:
    :param folder_cache:
    :return:
        a dataframe of experiments (4345, 2), where each record corresponds to a single gene expression
        a memmap of all experiments brain volumes, size (4345, 58, 41, 67) corresponding to
    (nexperiments, ml, dv, ap). The spacing between slices is 200 um
        a brainatlas object with the labels and coordinates matching the gene expression volumes
    """
    OLD_VERSIONS = ['2023-06-12']
    folder_cache = Path(folder_cache or atlas.AllenAtlas._get_cache_dir().joinpath('agea'))
    # check the AWS version and download the files if needed
    version_flag = next(folder_cache.glob('*.version'), None)
    if version_flag is None or version_flag.stem in OLD_VERSIONS:
        _logger.info(f'downloading gene expression data from {aws.S3_BUCKET_IBL} s3 bucket...')
        aws.s3_download_folder('atlas/agea', folder_cache)

    # load the genes dataframe and the gene expression volumes
    file_parquet = Path(folder_cache).joinpath('gene-expression.pqt')
    file_expression = Path(folder_cache).joinpath('gene-expression.bin')
    df_genes = pd.read_parquet(file_parquet)
    expression_volumes = np.memmap(file_expression, dtype=np.float16, mode='r', offset=0, shape=expression_size)

    # create a brain atlas object with the gene expression volume geometry and pre-computed labels
    atlas_agea = atlas.BrainAtlas(
        image=np.load(Path(folder_cache).joinpath('image.npy')),
        label=np.load(Path(folder_cache).joinpath('label.npy')),
        dxyz=200 / 1e6 * np.array([1, -1, -1]),
        regions=atlas.BrainRegions(),
        iorigin=atlas.ALLEN_CCF_LANDMARKS_MLAPDV_UM['bregma'] / 200 + np.array([0, 0, 0]),
        dims2xyz=[0, 2, 1],
        xyz2dims=[0, 2, 1]
    )

    return df_genes, expression_volumes, atlas_agea
