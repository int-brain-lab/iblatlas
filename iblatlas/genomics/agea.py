import logging
from pathlib import Path

import numpy as np
import pandas as pd

from iblutil.io.hashfile import md5
import one.remote.aws as aws

from iblatlas.atlas import AllenAtlas

_logger = logging.getLogger(__name__)

_, NML, NDV, NAP = DIM_EXP = (4345, 58, 41, 67)  # nexperiments, nml, ndv, nap


def load(folder_cache=None):
    """
    Reads in the Allen gene expression experiments binary data.
    Generation scripts from the Allen Institute are available
    [1] E. S. Lein et al., “Genome-wide atlas of gene expression in the adult mouse brain,”
     Nature, vol. 445, no. 7124, Art. no. 7124, Jan. 2007, doi: 10.1038/nature05453.
    [2] L. Ng et al., “An anatomic gene expression atlas of the adult mouse brain,”
     Nat Neurosci, vol. 12, no. 3, Art. no. 3, Mar. 2009, doi: 10.1038/nn.2281.

    :param filename:
    :param folder_cache:
    :return:
        a dataframe of experiments, where each record corresponds to a single gene expression
        a memmap of all experiments brain volumes, size (4345, 58, 41, 67) corresponding to
    (nexperiments, ml, dv, ap). The spacing between slices is 200 um
    """
    OLD_MD5 = []
    folder_cache = folder_cache or AllenAtlas._get_cache_dir().joinpath('agea')
    file_parquet = Path(folder_cache).joinpath('gene-expression.pqt')
    file_bin = file_parquet.with_suffix(".bin")

    if not file_parquet.exists() or md5(file_parquet) in OLD_MD5:
        _logger.info(f'downloading gene expression data from {aws.S3_BUCKET_IBL} s3 bucket...')
        aws.s3_download_file(f'atlas/{file_parquet.name}', file_parquet)
        aws.s3_download_file(f'atlas/{file_bin.name}', file_bin)
    df_genes = pd.read_parquet(file_parquet)
    gexp_all = np.memmap(file_bin, dtype=np.float16, mode='r', offset=0, shape=DIM_EXP)
    return df_genes, gexp_all




