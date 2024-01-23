"""
When looking at one or a few genes it is possible to upsample the gene expression volume
to the 25um or even the 10um resolution of the Allen CCF.
However when looking at thousand of genes, it is more feasible to keep the 200um resolution
Here we are downsampling a version of the CCF to 200um resolution that matches the gene expression atlas to
use as a label volume and to compute the brain coordinates of the gene expression voxels
"""
# %%
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats

from iblatlas import atlas

DIM_EXP = (4345, 58, 41, 67)  # nexperiments, nml, ndv, nap
#  path_agea_cache = atlas.AllenAtlas._get_cache_dir().joinpath('agea')
path_agea_cache = Path('/Users/olivier/Downloads/ONE/alyx.internationalbrainlab.org/histology/ATLAS/Needles/Allen/agea')
df_genes = pd.read_parquet(Path(path_agea_cache).joinpath('gene-expression.pqt'))
gexp_all = np.memmap(Path(path_agea_cache).joinpath('gene-expression.bin'), dtype=np.float16, mode='r', offset=0, shape=DIM_EXP)

# create a brain atlas object with the gene expression volume geometry
bg = atlas.BrainAtlas(
    image=gexp_all[0],
    label=gexp_all[0] != -1,
    dxyz=200 / 1e6 * np.array([1, -1, -1]),
    regions=atlas.BrainRegions(),
    iorigin=atlas.ALLEN_CCF_LANDMARKS_MLAPDV_UM['bregma'] / 200 + np.array([0, 0, 0]),
    dims2xyz=[0, 2, 1],
    xyz2dims=[0, 2, 1]
)
# instantiate the allen atlas at 25um resolution and compute the transform to the agea volume
ba = atlas.AllenAtlas(res_um=25)
xia2g = bg.bc.x2i(ba.bc.xscale)
yia2g = bg.bc.y2i(ba.bc.yscale)
zia2g = bg.bc.z2i(ba.bc.zscale)
# meshgrid the coordinates and aggregate per agea voxel, mean for image, and mode for the label
# (this takes less than a minute)
Xia2g, Yia2g, Zia2g = [aa.flatten() for aa in np.meshgrid(xia2g, yia2g, zia2g)]
df = pd.DataFrame(
    np.c_[Xia2g, Yia2g, Zia2g, ba.image.flatten(), ba.label.flatten()],
    columns=['ix', 'iy', 'iz', 'image', 'label'])
gene_expression_labels = df.groupby(['ix', 'iy', 'iz']).agg(
    image=pd.NamedAgg(column="image", aggfunc="mean"),
    label=pd.NamedAgg(column="label", aggfunc=lambda x: scipy.stats.mode(x)[0][0]),
).reset_index()

# create the image and label volume for the agea brain atlas
image = np.zeros(gexp_all[0].shape, dtype=np.float32)
label = np.zeros(gexp_all[0].shape, dtype=np.int16)
dim_columns = [['ix', 'iy', 'iz'][i] for i in bg.xyz2dims]
linear_indices = np.ravel_multi_index([gene_expression_labels[c].values for c in dim_columns], image.shape)
np.put(image, linear_indices, gene_expression_labels['image'].values)
np.put(label, linear_indices, gene_expression_labels['label'].values)

# make sure the transform went well by plotting a coronal and a sagittal slice from
# the same coordinate in both atlases
np.save(Path(path_agea_cache).joinpath('image.npy'), image)
np.save(Path(path_agea_cache).joinpath('label.npy'), label)
bg.image = image
bg.label = label

igenes = (0, 2325)
fig, axs = plt.subplots(3, 2, sharex=True, sharey=True)
ba.plot_cslice(0, ax=axs[0, 0])
bg.plot_cslice(0, ax=axs[1, 0])
ba.plot_cslice(0, ax=axs[0, 1], volume='annotation')
bg.plot_cslice(0, ax=axs[1, 1], volume='annotation')
bg.plot_cslice(0, ax=axs[2, 0], volume=gexp_all[igenes[0]], cmap='viridis')
bg.plot_cslice(0, ax=axs[2, 1], volume=gexp_all[igenes[1]], cmap='viridis')
fig.tight_layout()

fig, axs = plt.subplots(3, 2, sharex=True, sharey=True)
ba.plot_sslice(0, ax=axs[0, 0])
bg.plot_sslice(0, ax=axs[1, 0])
ba.plot_sslice(0, ax=axs[0, 1], volume='annotation')
bg.plot_sslice(0, ax=axs[1, 1], volume='annotation')
bg.plot_sslice(0, ax=axs[2, 0], volume=gexp_all[igenes[0]], cmap='viridis')
bg.plot_sslice(0, ax=axs[2, 1], volume=gexp_all[igenes[1]], cmap='viridis')
fig.tight_layout()

# %%
