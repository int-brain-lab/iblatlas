"""
This script was used to generate the following files
- depths_ind_10.npy
- depths_ind_25.npy
- depths_ind_50.npy
- depths_per_10.npy
- depths_per_25.npy
- depths_per_50.npy

These files are used by the function iblatlas.streamlines.utils.xyz_to_depth with the kwarg per=True to do a lookup of
xyz position to percentage depth from the surface of the cortex.

The original laplacian file that we work off (laplacian_10.nrrd) is provided by Allen and available to download
from here  https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/

(N.B we also have the IBL version of the Laplacian but have decided to use the Allen one for now)

Compared to the full laplacian volume provided by Allen we store
- The flattened index of voxels contained in the Isocortex (flattened index uses IBL order i.e AP, ML, DV)
    e.g depths_ind_10.npy
- The percentage depth from surface of cortex for each flattened index
    e.g depths_per_10.npy

We also subsample the laplacian volume to 25 um and 50 um resolution to generate the files for these resolution
atlases.
"""

import matplotlib.pyplot as plt
import numpy as np
import nrrd
import pandas as pd
from pathlib import Path

from iblatlas.atlas import AllenAtlas

file_path = Path('/Users/admin/Downloads/flatmap')

# -------------------------------------------------------------------------------------------------
# Prepare data for 10 um resolution
# -------------------------------------------------------------------------------------------------

# Load in the 10 um Allen Atlas
ba = AllenAtlas(10)
ba_shape = ba.image.shape
# Translate the label volume to cosmos mapping
idx = ba.regions.acronym2index('Isocortex')[1][0][0]
ba.label = ba.regions.mappings['Cosmos'][ba.label]
# Find the voxels that are contained within the isocortex, order AP, ML, DV
xyz = np.where(ba.label == idx)
# Compute the flattened idx that we will save
xyz_flat = np.ravel_multi_index(xyz, ba_shape).astype(np.int32)

del ba

# Load in the 10 um resolution laplacian, orientation AP, DV, ML
laplacian, _ = nrrd.read(file_path.joinpath('laplacian_10.nrrd'))
vals = laplacian[xyz[0], xyz[2], xyz[1]]
vals = vals.flatten().astype(np.float32)

del laplacian, xyz

# We just save the values for the voxels contained in isocortex
np.save(file_path.joinpath('depths_ind_10.npy'), xyz_flat)
np.save(file_path.joinpath('depths_per_10.npy'), vals)

# Check it makes sense
vol = np.full(np.prod(ba_shape), np.nan)
vol[xyz_flat] = vals
vol = vol.reshape(ba_shape)

ba = AllenAtlas(10)
ap = -2000 / 1e6
fig, ax = plt.subplots(figsize=(8, 5))
ax.imshow(np.moveaxis(ba.image[ba.bc.y2i(ap), :, :], 0, 1),
          extent=ba.extent(axis=1), cmap='gray')
im = ax.imshow(np.moveaxis(vol[ba.bc.y2i(ap), :, :], 0, 1),
               extent=ba.extent(axis=1))
plt.show()

del vol, xyz_flat, vals, ba

# -------------------------------------------------------------------------------------------------
# Downsample to 25 um and 50 um
# -------------------------------------------------------------------------------------------------

ba = AllenAtlas(10)
xyz_flat = np.load(file_path.joinpath('depths_ind_10.npy'))
vals = np.load(file_path.joinpath('depths_per_10.npy'))
# unravel so that they are in ixyz
ixyz = np.unravel_index(xyz_flat, ba.image.shape)
# Find the cartesian coordinates of all indices that we have depths for
xyz = ba.bc.i2xyz(np.c_[ixyz[1], ixyz[0], ixyz[2]])

del ba, ixyz

for res in [25, 50]:

    # Transform to 25 um atlas
    ba = AllenAtlas(res)
    df = pd.DataFrame()
    df['vals'] = vals

    # Find the index in the down-sampled volume for the cartesian coordinates
    ixyz = ba.bc.xyz2i(xyz)
    df['ind'] = np.ravel_multi_index((ixyz[:, 1], ixyz[:, 0], ixyz[:, 2]), ba.image.shape)

    # Group by index and find the median value
    df_grp = df.groupby('ind')
    df_med = df_grp.vals.median()
    vals_res = df_med.values
    xyz_flat_res = df_med.index.values
    np.save(file_path.joinpath(f'depths_ind_{res}.npy'), xyz_flat_res.astype(np.int32))
    np.save(file_path.joinpath(f'depths_per_{res}.npy'), vals_res.astype(np.float32))

    vol = np.full(np.prod(ba.image.shape), np.nan)
    vol[xyz_flat_res] = vals_res
    vol = vol.reshape(ba.image.shape)

    ap = -2000 / 1e6
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.imshow(np.moveaxis(ba.image[ba.bc.y2i(ap), :, :], 0, 1), extent=ba.extent(axis=1), cmap='gray')
    im = ax.imshow(np.moveaxis(vol[ba.bc.y2i(ap), :, :], 0, 1),
                   extent=ba.extent(axis=1))
    plt.show()


# -------------------------------------------------------------------------------------------------
# Check that the laplacian looks okay across the different resolutions
# -------------------------------------------------------------------------------------------------

fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
ap = -2000 / 1e6

for i, res in enumerate([10, 25, 50]):
    xyz_flat = np.load(file_path.joinpath(f'depths_ind_{res}.npy'))
    vals = np.load(file_path.joinpath(f'depths_per_{res}.npy'))
    ba = AllenAtlas(res)
    ba_shape = ba.image.shape

    vol = np.full(np.prod(ba_shape), np.nan)
    vol[xyz_flat] = vals
    vol = vol.reshape(ba_shape)

    ax = axs[i]
    ax.imshow(np.moveaxis(ba.image[ba.bc.y2i(ap), :, :], 0, 1), extent=ba.extent(axis=1), cmap='gray')
    im = ax.imshow(np.moveaxis(vol[ba.bc.y2i(ap), :, :], 0, 1),
                   extent=ba.extent(axis=1))

plt.show()
