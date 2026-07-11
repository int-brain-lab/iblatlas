import logging
from pathlib import Path

import pandas as pd

import one.remote.aws as aws

from iblatlas import atlas

_logger = logging.getLogger(__name__)

FILENAMES = ('allen_mouse_regionalized_connectivity.pqt', 'allen_mouse_hierarchy_scores.pqt')
# bump this and append the previous value to OLD_VERSIONS whenever the files on S3 change
CURRENT_VERSION = '2026-07-11'
OLD_VERSIONS = []


def _sync(folder_cache):
    """
    Download the connectivity dataset if the local cache is missing or stale.

    Follows the same `<date>.version` sentinel convention as `iblatlas.genomics.agea.load` /
    `iblatlas.genomics.merfish.load`: a version flag file is dropped in the cache folder after
    a full download, and re-downloaded only if that flag is missing or in `OLD_VERSIONS`. Both
    files are small (a few MB total) so, unlike those two datasets, we always sync them
    together rather than checking per-file existence.
    """
    folder_cache.mkdir(parents=True, exist_ok=True)
    version_flag = next(folder_cache.glob('*.version'), None)
    if version_flag is None or version_flag.stem in OLD_VERSIONS:
        _logger.info(f'downloading mesoscale connectivity data from {aws.S3_BUCKET_IBL} s3 bucket...')
        for filename in FILENAMES:
            aws.s3_download_file(f'atlas/connectivity/{filename}', folder_cache.joinpath(filename))
        folder_cache.joinpath(f'{CURRENT_VERSION}.version').touch()


def load(folder_cache=None):
    """
    Reads in the regionalized mesoscale connectivity matrix.

    This is the Oh et al. (2014) / Knox et al. (2018) regularized-regression connectivity
    model -- fit jointly across all injection experiments to statistically correct for
    injection-site overlap contamination, rather than a naive per-experiment average --
    regionalized onto 291 Allen "summary structures".

    :param folder_cache:
    :return: a tidy long dataframe (671_628, 9), one row per (source structure, hemisphere,
     target structure, metric), with columns:
      - source_structure_id / source_acronym: source (injection-side) structure
      - target_structure_id / target_acronym: target structure
      - hemisphere: 'ipsi' or 'contra' (relative to the source)
      - metric: one of 'connection_strength', 'connection_density',
        'normalized_connection_strength', 'normalized_connection_density'
      - value: metric value for this (source, hemisphere, target) triple
      - source_volume_mm3 / target_volume_mm3: single-hemisphere volume (mm^3) of the
        source/target structure, including descendants (e.g. layers)

    Notes
    -----
    Only `connection_strength` is additive across a re-parcellation (it is proportional to
    integrated axon volume): summing it within each new group, on both the source and target
    side, gives the exact value for that group -- provided the new parcellation partitions the
    291 summary structures with no overlap or gaps (true of any Allen ontology grouping, e.g.
    Beryl or Cosmos, via `iblatlas.regions.BrainRegions.id2id`). `connection_density` and the
    normalized metrics are ratios: re-derive them from the aggregated `connection_strength` and
    volumes at the new parcellation, rather than summing or averaging the ratio columns
    directly.

    Example
    -------
    Reaggregate onto the 10 Cosmos regions:

    >>> from iblatlas.regions import BrainRegions
    >>> br = BrainRegions()
    >>> df = load()
    >>> strength = df[df['metric'] == 'connection_strength'].copy()
    >>> strength['source_cosmos'] = br.id2acronym(br.id2id(strength['source_structure_id'].values, mapping='Cosmos'))
    >>> strength['target_cosmos'] = br.id2acronym(br.id2id(strength['target_structure_id'].values, mapping='Cosmos'))
    >>> cosmos = strength.groupby(['source_cosmos', 'target_cosmos', 'hemisphere'])['value'].sum()
    """
    folder_cache = Path(folder_cache or atlas.AllenAtlas._get_cache_dir().joinpath('connectivity'))
    _sync(folder_cache)
    return pd.read_parquet(folder_cache.joinpath('allen_mouse_regionalized_connectivity.pqt'))


def load_hierarchy(folder_cache=None):
    """
    Reads in the Harris et al. (2019) precomputed cortical/thalamic hierarchy scores.

    :param folder_cache:
    :return: a dataframe (122, 9), one row per (area, correction scheme), with columns:
      - correction: 'cre_conf' (Cre-line confidence weighted) or 'no_conf' (unweighted)
      - area: structure acronym (e.g. 'VISp'), matching `source_acronym` / `target_acronym`
        returned by `load()`
      - region_type: 'C' = cortex, 'T' = thalamus
      - cc_before / cc_iter: hierarchy score using cortico-cortical connections only,
        before/after iterative refinement
      - cctc_before / cctc_iter: + thalamo-cortical connections
      - cctcct_before / cctcct_iter: + cortico-thalamic connections (full model used in the
        paper)

    Lower score = earlier/more feedforward in the hierarchy (e.g. VISp is near the bottom of
    the visual hierarchy); higher = later/more association-like.
    """
    folder_cache = Path(folder_cache or atlas.AllenAtlas._get_cache_dir().joinpath('connectivity'))
    _sync(folder_cache)
    return pd.read_parquet(folder_cache.joinpath('allen_mouse_hierarchy_scores.pqt'))
