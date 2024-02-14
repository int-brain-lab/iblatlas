import logging
from pathlib import Path

import pandas as pd
import numpy as np

import one.remote.aws as aws

from iblatlas import atlas

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
    OLD_VERSIONS = ['2023-06-12']
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
