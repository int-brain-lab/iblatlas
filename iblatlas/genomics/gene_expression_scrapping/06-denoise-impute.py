"""
This script will pre-process the raw AGEA volumes and write a newer version to disk.
The pre-processing steps include:
- merge hemispheres
- impute using PPCA
- remove curtaining effect trying to minimize the slice to slice amplitude variation per Cosmos regions
"""

from pathlib import Path
import tqdm
import joblib

import numpy as np
import pandas as pd
import scipy.ndimage
import seaborn as sns
import matplotlib.pyplot as plt

import scipy.optimize
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import pyppca  # pip install pyppca
from iblatlas.genomics import agea

df_genes, gene_expression_volumes, atlas_agea = agea.load()
mask_brain = atlas_agea.label == 0

folder_volumes = Path('/datadisk/Data/2025/denoise_agea')
folder_volumes = Path('/mnt/s1/2025/denoise_agea')


def convert_to_dual_volume(input_memmap, output_memmap):
    # Creates the dual volume
    for igene in tqdm.tqdm(np.arange(input_memmap.shape[0])):
        # igene = 625
        gvol = np.copy(input_memmap[igene]).astype(float)
        gvol[gvol < 0] = np.nan
        gvol_dual = np.nanmean(np.stack((gvol, np.flip(gvol, axis=0)), axis=3), axis=3)
        gvol_dual[np.isnan(gvol_dual)] = - 1
        output_memmap[igene] = gvol_dual.astype(np.float16)  # isl = 31


def compute_agea_pca(gene_expression_volumes, atlas_agea, pca_n_components=50):
    """
    Compute PCA embedding for AGEA and return the embedding and the indices of the voxels within the brain volume.

    :param gene_expression_volumes: numpy array of shape (n_genes, n_ml, n_dv, n_ap)
    :param atlas_agea: AGEA atlas object
    :param pca_n_components: Number of components for PCA
    :return: PCA embedding (n_voxels_embedding, n_pca_components)
     and flat indices of selected voxels within the brain volume np.array(int64) (n_voxels_embedding)
    :return:
    """
    inside_idx = np.where(atlas_agea.label.flatten() != 0)[0]
    ng = gene_expression_volumes.shape[0]
    # sel = atlas_agea.label.flatten() != 0  # remove void voxels
    # reshape in a big array nexp x nvoxels this takes a little while
    gexps = np.copy(gene_expression_volumes.reshape((ng, -1))[:, inside_idx].astype(np.float32).transpose())
    p_missing = np.mean(gexps < 0, axis=1)  # this is the proportion of missing genes
    gexps = gexps[p_missing < 0.6, :]  # we select only the voxels with less than 60% missing genes
    embedding_idx = inside_idx[p_missing < 0.6]  # we select only voxels within the brain volume
    # % run PCA on gexps
    scaler = StandardScaler()
    gexps = scaler.fit_transform(gexps)
    pca = PCA(n_components=pca_n_components)
    embedding = pca.fit_transform(gexps)
    return embedding, embedding_idx


def compute_agea_ppca(input_memmap, atlas_agea, output_memmap=None, ppca_n_components=50):
    inside_idx = np.where(atlas_agea.label.flatten() != 0)[0]
    outside_idx = np.unravel_index(np.where(atlas_agea.label.flatten() == 0), input_memmap.shape[1:])
    ng = input_memmap.shape[0]
    gexps = input_memmap.reshape((ng, -1))[:, inside_idx].astype(np.float32).transpose()
    gexps[gexps < 0] = np.nan
    C, ss, M, X, Ye = pyppca.ppca(gexps, d=ppca_n_components, dia=True)
    output_memmap = np.copy(input_memmap) if output_memmap is None else output_memmap
    output_memmap[:, *np.unravel_index(inside_idx, input_memmap.shape[1:])] = Ye.T  # noqa
    output_memmap[:, *outside_idx] = input_memmap[:, *outside_idx]
    return output_memmap


def train_region_predictor(atlas_agea, embedding, embedding_idx, mapping=None):
    """
    From the PCA embedding, predict the cosmos level label using a MLP.
    split the data in training and testing sets
    Accuracy: 0.7351404310907903 for allen regions
    Accuracy: 0.9329686479425212 for cosmos regions
    :param atlas_agea:
    :param embedding:
    :param embedding_idx:
    :param mapping:
    :return:
    """
    aids = np.abs(atlas_agea.regions.id[atlas_agea.label.flatten()[embedding_idx]])
    if mapping is not None:
        labels = atlas_agea.regions.remap(aids, source_map='Allen', target_map='Cosmos')
    else:
        labels = aids
    X_train, X_test, y_train, y_test = train_test_split(embedding, labels, test_size=0.2, random_state=42)
    mlp = MLPClassifier(hidden_layer_sizes=(50,), max_iter=300)
    mlp.fit(X_train, y_train)
    # y_pred = mlp.predict(X_test)
    accuracy = mlp.score(X_test, y_test)
    print("Accuracy:", accuracy)
    return accuracy


def curtaining(volume_agea, max_ratio=2):
    """
    Removes curtaining effect from a single gene expression volume
    :param volume_agea:
    :param max_ratio:
    :return:
    """

    def objective_function(weights, expression, counts):
        weighted_X = np.log(expression * weights[:, np.newaxis])
        std_x = np.abs(weighted_X - np.nanmedian(weighted_X, axis=0))  # Compute the weighted mean of the expression
        return np.nansum(std_x * counts)

    df_vol = pd.DataFrame({
        'iy': np.tile(np.arange(volume_agea.shape[2]), np.prod(volume_agea.shape[:2])),
        'gene': volume_agea.flatten(),
        'rindex': atlas_agea.label.flatten()}
    )
    df_vol = df_vol.loc[np.logical_and(df_vol['gene'] >= 0, df_vol['rindex'] != 0)]
    df_vol['allen_id'] = np.abs(atlas_agea.regions.id[df_vol['rindex']])
    df_vol['cosmos_id'] = atlas_agea.regions.remap(df_vol['allen_id'].to_numpy(), source_map='Allen',
                                                   target_map='Cosmos')
    # here maybe median is more robust than mean for outliers
    df_slices = df_vol.pivot_table(index='iy', values='gene', aggfunc=['median', 'count'], columns='cosmos_id')
    expression = df_slices.loc[:, 'median'].to_numpy()
    counts = df_slices.loc[:, 'count'].to_numpy()
    n_slices = expression.shape[0]
    # Optimization
    result = scipy.optimize.minimize(
        objective_function,
        np.ones(n_slices),
        args=(expression, counts),
        method='L-BFGS-B',
        bounds=[(1 / max_ratio, max_ratio) for _ in range(n_slices)]
    )
    coronal_weights = np.ones(atlas_agea.bc.ny)
    coronal_weights[df_slices.index] = result['x']
    return volume_agea * (coronal_weights[np.newaxis, np.newaxis, :])


def curtaining_parallel(input_memmap, output_memmap=None):
    """
    Try to remove the coronal slice curtaining effect in the AGEA volume.
    The idea is to find the coronal slices weights that minimize the standard deviation of the
     expression accross voxels within a cosmos regions
    :return:
    """
    output_memmap = np.copy(input_memmap) if output_memmap is None else output_memmap

    def compute_single_gene(igene):
        output_memmap[igene] = curtaining(input_memmap[igene])
    pass

    ng = input_memmap.shape[0]
    jobs = (joblib.delayed(compute_single_gene)(igene) for igene in np.arange(ng))
    list(tqdm.tqdm(joblib.Parallel(n_jobs=joblib.cpu_count() - 1, return_as='generator')(jobs), total=ng))


# %%
RECOMPUTE = False
files = {
    'dual': folder_volumes.joinpath('dual_memmap.bin'),
    'ppca': folder_volumes.joinpath('ppca.bin'),
    'dual_ppca': folder_volumes.joinpath('dual_ppca.bin'),
    'ppca_dual': folder_volumes.joinpath('ppca_dual.bin'),
    'dual_ppca_curtain': folder_volumes.joinpath('dual_ppca_curtain.bin'),
    'ppca_dual_curtain': folder_volumes.joinpath('ppca_dual_curtain.bin')
}
if RECOMPUTE:
    print('Create memmaps...')
    memmaps = {k: np.memmap(v, dtype=np.float16, mode='w+', offset=0,
                            shape=gene_expression_volumes.shape) for k, v in files.items()}
    convert_to_dual_volume(gene_expression_volumes, memmaps['dual'])
    compute_agea_ppca(memmaps['dual'], atlas_agea, output_memmap=memmaps['dual_ppca'])
    compute_agea_ppca(gene_expression_volumes, atlas_agea, output_memmap=memmaps['ppca'])
    convert_to_dual_volume(memmaps['ppca'], memmaps['ppca_dual'])
    curtaining_parallel(memmaps['dual_ppca'], output_memmap=memmaps['dual_ppca_curtain'])
    curtaining_parallel(memmaps['ppca_dual'], output_memmap=memmaps['ppca_dual_curtain'])
else:
    memmaps = {k: np.memmap(v, dtype=np.float16, mode='r+', offset=0,
                            shape=gene_expression_volumes.shape) for k, v in files.items()}
memmaps['baseline'] = gene_expression_volumes

# %%
file_results = folder_volumes.joinpath('results.csv')
if file_results.exists() and not RECOMPUTE:
    df_results = pd.read_csv(file_results)
    to_skip = list(df_results.volume.unique())
else:
    df_results = pd.DataFrame(columns=['volume', 'region', 'accuracy'])
    to_skip = []

results = []
for k, vol in memmaps.items():
    if k in to_skip:
        continue
    print(f'Compute PCA embeddings for {k}...')
    embedding, embedding_idx = compute_agea_pca(memmaps[k], atlas_agea)
    print(f'Compute Allen region prediction for {k}...')
    accuracy = train_region_predictor(atlas_agea, embedding, embedding_idx, mapping=None)
    results.append(dict(volume=k, region='Allen', accuracy=accuracy))
    print(f'Compute Cosmos region prediction for {k}...')
    accuracy = train_region_predictor(atlas_agea, embedding, embedding_idx, mapping='Cosmos')
    results.append(dict(volume=k, region='Cosmos', accuracy=accuracy))

if len(results) > 0:
    df_results = pd.concat([df_results, pd.DataFrame(results)], ignore_index=True, axis=0)
    df_results.to_csv(file_results, index=False)


# %% plot the accuracy for each method and region
volume_order = ['baseline', 'ppca', 'dual', 'ppca_dual', 'ppca_dual_curtain']
sns.catplot(x='region', y='accuracy', hue='volume', data=df_results, kind='bar', hue_order=volume_order)
plt.show()

# aws s3 cp /mnt/s1/2025/denoise_agea/ppca_dual_curtain.bin s3://ibl-brain-wide-map-public/atlas/agea/gene-expression-processed.bin  # NOQA
# touch /mnt/s1/2025/denoise_agea/2025-03-18.version
# aws s3 cp /mnt/s1/2025/denoise_agea/2025-03-18.version s3://ibl-brain-wide-map-public/atlas/agea/2025-03-18.version  # NOQA