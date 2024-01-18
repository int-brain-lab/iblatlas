from iblatlas.genomics import agea
df_genes, gene_expression_volumes = agea.load()

# this is the dataframe of 4345 genes (experiments) with their ids and acronyms
print(df_genes)

# for each gene, we have a 3D volume of gene expression values mapped at 200um resolution
print("n_volumes, ml, dv, ap", gene_expression_volumes.shape)

## %% relate the gene expression volumes to the Allen CCF
from iblatlas import atlas
import numpy as np
import matplotlib.pyplot as plt


igene = 456
gba = atlas.BrainAtlas(
    image=gene_expression_volumes[igene],
    label=gene_expression_volumes[igene] != -1,
    dxyz=200 / 1e6,
    regions=atlas.BrainRegions(),
    iorigin=atlas.ALLEN_CCF_LANDMARKS_MLAPDV_UM['bregma'] / 200 + np.array([0, 0, 0]),
    dims2xyz=[0, 2, 1],
    xyz2dims=[0, 2, 1]
)
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)

gba.plot_cslice(0.0, ax=axs[0], cmap='magma')
gba.plot_sslice(0.0, ax=axs[1], cmap='magma')

## %%

ccf_coords = np.meshgrid(*[np.arange(DIM_EXP[i]) * 200 for i in [1, 2, 3]])  # ml, dv, ap
xyzs = ba.ccf2xyz(np.c_[ccf_coords[0].flatten(), ccf_coords[2].flatten(), ccf_coords[1].flatten()])
aids = ba.get_labels(xyzs, mode='clip', mapping='Cosmos')


ba = atlas.AllenAtlas()

fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
ba.plot_cslice(0.0, ax=axs[0], cmap='magma')
ba.plot_sslice(0.0, ax=axs[1], cmap='magma')