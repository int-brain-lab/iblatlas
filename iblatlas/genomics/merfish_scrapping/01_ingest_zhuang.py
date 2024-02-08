"""
This script is used to ingest the MERFISH data from the Zhuang lab and the MERFISH subject
The purpose is to turn large csv files into parquet files for faster access.
Also the ~8M resulting cells table are not joined automatically to all the levels of hierarchy
"""

from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.colors

from iblutil.util import setup_logger
from iblatlas.atlas import AllenAtlas
from iblutil.numerical import ismember

ba = AllenAtlas()
logger = setup_logger('ibllib', level='INFO')
RELEASE = '20231215'
path_allen_cache = Path("/datadisk/Data/merfish_atlas/cache")
path_ibl_cache = Path("/datadisk/Data/merfish_atlas/ibl")
# path_allen_cache = Path("/Users/olivier/Documents/datadisk/Data/merfish_atlas/cache")
# path_ibl_cache = Path("/Users/olivier/Documents/datadisk/Data/merfish_atlas/ibl")

subjects = ['MERFISH-C57BL6J-638850', 'Zhuang-ABCA-1', 'Zhuang-ABCA-2', 'Zhuang-ABCA-3', 'Zhuang-ABCA-4']
df_subjects = pd.DataFrame(columns=['donor_label', 'donor_sex', 'donor_genotype']).set_index('donor_label')
df_neurotransmitters = pd.DataFrame(columns=['neurotransmitter', 'neurotransmitter_color']).set_index('neurotransmitter')
df_class = pd.DataFrame(columns=['class_id', 'class', 'class_color']).set_index('class_id')


def hex2rgba_dataframe(df, color_field, target_field):
    # stores rgba in a uint32 column for faster joins
    rgba = np.concatenate(df[color_field].apply(lambda x: np.array(matplotlib.colors.to_rgba(x))[np.newaxis, :]).values)
    df[target_field] = (rgba * 255).astype(np.uint8).view('uint32')
    return df


def int_after_colon(text):
    if isinstance(text, float):
        return text
    pattern = r'[a-zA-Z]+:(\d+)'
    match = re.search(pattern, text)
    if match:
        captured_integer = int(match.group(1))
        return float(captured_integer)
    else:
        return np.NaN


def fcn_groupby_check_unique(df_cells, group, fields):
    # fields = ['class', 'subclass_color']
    # group = 'subclass'
    df = df_cells.groupby(group).agg(
        **({k: pd.NamedAgg(column=k, aggfunc='first') for k in fields} | {
            f"{k}_unique": pd.NamedAgg(column=k, aggfunc='nunique') for k in fields}
           )
    )
    assert np.all(df.loc[:, [f'{f}_unique' for f in fields]].values == 1)
    df = df.drop(columns=[f'{f}_unique' for f in fields])
    return df


def fcn_integer_index(df, column, ndigits=2):
    df = df.reset_index()
    df = df.set_index(df[column].apply(lambda x: int(x[:ndigits])))
    df.index.rename(f'{column}_id', inplace=True)
    return df


def reindex_missing_rows(df):
    # add missing rows so that linear indexing matches pandas indexing for speed purposes
    # problem is to handle the types of inserted recorcds correctly so they export to parquet
    imiss = np.setxor1d(np.arange(np.max(df.index) + 1), df.index)
    missing_rec = df.iloc[df.index[0], :].copy()
    for k in missing_rec.keys():
        if isinstance(missing_rec[k], str):
            missing_rec[k] = ''
        elif isinstance(missing_rec[k], int | np.uint32 | np.int32 | np.uint64 | np.int64):
            missing_rec[k] = 0
        elif isinstance(missing_rec[k], float | np.float32 | np.float64):
            missing_rec[k] = np.NaN
    df = df.reindex(pd.Index(np.arange(np.max(df.index) + 1), name=df.index.name), fill_value=0)
    df.iloc[imiss, :] = missing_rec
    return df


ldf_neurotransmitters = []
ldf_classes = []
ldf_genes = []
ldf_subclasses = []
ldf_supertypes = []
ldf_clusters = []

for subject in subjects:

    path_subject = path_allen_cache.joinpath('metadata', subject, RELEASE)
    logger.warning(path_subject)
    # get the genes, as it happens the data is similar for all subjects
    df_genes = pd.read_csv(path_subject.joinpath('gene.csv'))
    df_genes = df_genes.set_index('gene_identifier')
    df_genes['mapped_ncbi_identifier'].map(int_after_colon)
    ldf_genes.append(df_genes)

    # export the csv in parquet
    logger.info('read cells')
    df_cells = pd.read_csv(path_subject.joinpath('views', 'cell_metadata_with_cluster_annotation.csv'))
    df_cells.set_index('cell_label', inplace=True)

    logger.info('read CCF coordinates and join on cells dataframe')
    # the CCF coordinates release is hard coded as there was no revision in 20231215 for this dataset
    file_csv_ccf = path_allen_cache.joinpath('metadata', f"{subject}-CCF", RELEASE, 'ccf_coordinates.csv')
    if not file_csv_ccf.exists():
        file_csv_ccf = path_allen_cache.joinpath('metadata', f"{subject}-CCF", '20230830', 'ccf_coordinates.csv')
    df_ccf = pd.read_csv(file_csv_ccf)
    logger.info(f'{df_cells.shape[0]:_} cells original, {df_ccf.shape[0]:_} registered in CCF')
    df_ccf.set_index('cell_label', inplace=True)
    df_ccf.rename(columns={'x': 'x_ccf', 'y': 'y_ccf', 'z': 'z_ccf'}, inplace=True)
    df_cells = df_cells.join(df_ccf, how='inner')
    logger.info(f'{df_cells.shape[0]:_} cells remaining after join with CCF')

    # update the subjects table
    logger.info('get subjects')
    assert np.all(df_cells['feature_matrix_label'] == df_cells['donor_label'])
    assert np.all(df_cells.groupby('donor_label').agg(
        donor_genotype=pd.NamedAgg(column="donor_genotype", aggfunc="nunique"),
        donor_sex=pd.NamedAgg(column="donor_sex", aggfunc="nunique"),
    ).values == 1)
    subject = df_cells['donor_label'].iloc[0]
    df_subjects.at[df_cells['donor_label'].iloc[0], 'donor_genotype'] = df_cells['donor_genotype'].iloc[0]
    df_subjects.at[df_cells['donor_label'].iloc[0], 'donor_sex'] = df_cells['donor_sex'].iloc[0]

    logger.info('neurotransmitters')
    # update the neurotransmitters table
    df_neurotransmitters = fcn_groupby_check_unique(df_cells, 'neurotransmitter', ['neurotransmitter_color'])
    ldf_neurotransmitters.append(df_neurotransmitters)

    logger.info('classes')
    # update the classes table, the neurotransmitters are not unique per class
    df_class = fcn_groupby_check_unique(df_cells, 'class', ['class_color'])
    df_class = fcn_integer_index(df_class, 'class', ndigits=2)
    ldf_classes.append(df_class)

    logger.info('subclasses')
    # update the subclass table, the neurotransmitters are not unique per subclass either
    df_subclass = fcn_groupby_check_unique(df_cells, 'subclass', ['subclass_color', 'class'])
    df_subclass = fcn_integer_index(df_subclass, 'subclass', ndigits=3)
    ldf_subclasses.append(df_subclass)

    logger.info('supertype')
    # update the supertype table, the neurotransmitters are not unique per supertype either
    df_sp = fcn_groupby_check_unique(df_cells, 'supertype', ['supertype_color', 'subclass'])
    df_sp = fcn_integer_index(df_sp, 'supertype', ndigits=4)
    ldf_supertypes.append(df_sp)

    logger.info('clusters')
    df_clusters = fcn_groupby_check_unique(df_cells, 'cluster', ['cluster_color', 'supertype', 'cluster_alias'])
    df_clusters = fcn_integer_index(df_clusters, 'cluster', ndigits=4)
    ldf_clusters.append(df_clusters)

    logger.info('cells')

    df_cells = df_cells.drop(columns=[
        'donor_genotype', 'donor_sex', 'feature_matrix_label', 'neurotransmitter_color', 'class_color', 'cluster_color',
        'cluster_alias', 'supertype_color', 'subclass_color'])
    df_cells['cluster'] = df_cells['cluster'].apply(lambda x: int(x[:4]))
    df_cells['supertype'] = df_cells['supertype'].apply(lambda x: int(x[:4]))
    df_cells['subclass'] = df_cells['subclass'].apply(lambda x: int(x[:3]))
    df_cells['class'] = df_cells['class'].apply(lambda x: int(x[:3]))
    # x: ml, y: dv, z: ap, x_ccf: ap, y_ccf: dv, z_ccf = ml,
    xyz = ba.ccf2xyz(df_cells.loc[:, ['z_ccf', 'x_ccf', 'y_ccf']].values * 1e3)
    aids = ba.get_labels(xyz, mode='clip', mapping='Allen')
    _, rids = ismember(aids, ba.regions.id)
    df_cells.drop(columns=['x', 'y', 'z', 'x_ccf', 'y_ccf', 'z_ccf', 'parcellation_index'], inplace=True)
    if 'cluster_confidence_score' in df_cells.columns:
        # the cluster_confidence_score and high_quality_transfer are not labeled for all subjects
        # so I decided to drop them as they are not complete
        df_cells.drop(columns=['cluster_confidence_score', 'high_quality_transfer'], inplace=True)
    if 'average_correlation_score' in df_cells.columns:
        df_cells.drop(columns=['average_correlation_score'], inplace=True)
    df_cells['x'] = xyz[:, 0]
    df_cells['y'] = xyz[:, 1]
    df_cells['z'] = xyz[:, 2]
    df_cells['Allen_id'] = aids
    df_cells.to_parquet(path_ibl_cache.joinpath(f"{subject}_cells.pqt"))


## %% writing the aggregate datasets for all subjects

df_neurotransmitters = pd.concat(ldf_neurotransmitters).drop_duplicates()
df_neurotransmitters = hex2rgba_dataframe(df_neurotransmitters, 'neurotransmitter_color', 'neurotransmitter_rgba')
df_neurotransmitters.to_parquet(path_ibl_cache.joinpath("neurotransmitters.pqt"))

df_genes = pd.concat(ldf_genes).drop_duplicates()
df_genes.to_parquet(path_ibl_cache.joinpath("genes.pqt"))

df_classes = pd.concat(ldf_classes).drop_duplicates()
df_classes = hex2rgba_dataframe(df_classes, 'class_color', 'class_rgba')
df_classes = reindex_missing_rows(df_classes)
df_classes.to_parquet(path_ibl_cache.joinpath("classes.pqt"))

df_subclasses = pd.concat(ldf_subclasses).drop_duplicates()
df_subclasses = hex2rgba_dataframe(df_subclasses, 'subclass_color', 'subclass_rgba')
df_subclasses = reindex_missing_rows(df_subclasses)
df_subclasses.to_parquet(path_ibl_cache.joinpath("subclasses.pqt"))

df_supertypes = pd.concat(ldf_supertypes).drop_duplicates()
df_supertypes = hex2rgba_dataframe(df_supertypes, 'supertype_color', 'supertype_rgba')
df_supertypes = reindex_missing_rows(df_supertypes)
df_supertypes.to_parquet(path_ibl_cache.joinpath("supertypes.pqt"))

df_clusters = pd.concat(ldf_clusters).drop_duplicates()
df_clusters = hex2rgba_dataframe(df_clusters, 'cluster_color', 'cluster_rgba')
df_clusters = reindex_missing_rows(df_clusters)
df_clusters.to_parquet(path_ibl_cache.joinpath("clusters.pqt"))
