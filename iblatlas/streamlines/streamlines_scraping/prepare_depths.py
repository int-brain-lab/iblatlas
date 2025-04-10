"""
This script was used to generate the following files
- depths_um_10.npy
- depths_um_25.npy
- depths_um_50.npy

These files are used by the function iblatlas.streamlines.utils.xyz_to_depth with the kwarg per=False to do a lookup of
xyz position to depth in um from the surface of the cortex.

The original streamlines file that we work off (streamlines_allen.npy) is available to download from here
https://drive.google.com/drive/u/0/folders/1MbO3n2LAPEocDTDwRKCCJuVyNSvCcaG2

(N.B we also have the IBL version of the streamlines but have decided to use the Allen one for now)

There are cases where multiple streamlines pass through the same voxel. In these cases the voxel has been
assigned the median depth value

For voxels that don't have any streamlines that pass through them the values are interpolated using nearest neighbour

We also subsample the depth volume to 25 um and 50 um resolution to generate the files for these resolution
atlases.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import NearestNDInterpolator
from iblatlas.atlas import AllenAtlas

file_path = Path('/Users/admin/int-brain-lab/flatmap')

# -------------------------------------------------------------------------------------------------
# Depth computation for 10 um
# -------------------------------------------------------------------------------------------------

ba = AllenAtlas(10)
# Order is AP, ML, DV
ba_shape = ba.image.shape

# Convert label volume to cosmos mapping
idx = ba.regions.acronym2index('Isocortex')[1][0][0]
ba.label = ba.regions.mappings['Cosmos'][ba.label]

# Find the voxels that lie in Isocortex AP, ML, DV
xyz = np.where(ba.label == idx)
# Compute the flattened idx for later
xyz_flat = np.ravel_multi_index(xyz, ba_shape)

streamlines = np.load(file_path.joinpath('streamlines_allen.npy'))
# Reorder to match atlas image so AP, ML, DV
streamlines = streamlines[:, :, [0, 2, 1]]

# Compute cartesian coordinate of top of each streamline
mlapdv_top = np.c_[ba.bc.i2x(streamlines[:, 0, 1]),
                   ba.bc.i2y(streamlines[:, 0, 0]),
                   ba.bc.i2z(streamlines[:, 0, 2])]
# Compute cartesian coordinate of bottom of each streamline
mlapdv_bottom = np.c_[ba.bc.i2x(streamlines[:, -1, 1]),
                      ba.bc.i2y(streamlines[:, -1, 0]),
                      ba.bc.i2z(streamlines[:, -1, 2])]

# Compute the distance in um of the interval along each streamline
steps_um = np.sqrt(np.sum((mlapdv_top - mlapdv_bottom) ** 2, axis=1)) * 1e6 / streamlines.shape[1]

del mlapdv_top, mlapdv_bottom, ba

# Get rid of streamlines that have no distance
kp_strm = np.where(steps_um > 0)[0]
streamlines = streamlines[kp_strm]
steps_um = steps_um[kp_strm]

del kp_strm

# Make an array of depths along each streamline
steps = np.arange(100) + 1
steps = np.tile(steps.reshape(1, 100), (streamlines.shape[0], 1))
depths = (steps * steps_um.reshape(steps_um.size, 1)).astype(np.float32)

del steps, steps_um

# Create a dataframe
df = pd.DataFrame()
df['inds'] = np.ravel_multi_index(
    (streamlines[:, :, 0], streamlines[:, :, 1], streamlines[:, :, 2]), ba_shape).flatten()
df['depths'] = depths.flatten()

del depths, streamlines

# Group by index in volume
df_grp = df.groupby('inds')

# Find the median value per voxel
agg_val = df_grp.median()

# Create a volume and set the values to the median values
vol = np.zeros(np.prod(ba_shape), dtype=np.float32)
vol[xyz_flat] = np.nan
vol[agg_val.index.values] = agg_val.depths.values
vol = vol.reshape(ba_shape)

# Find the indices where we have data (only consider Isocortex)
coords = np.where(np.bitwise_and(~np.isnan(vol), vol > 0))
vals = vol[coords[0], coords[1], coords[2]]
coords = np.array([coords[0], coords[1], coords[2]]).T

# These are the indices in the Isocortex where we don't have values for therefore need to interpolate
fill_coords = np.where(np.isnan(vol))
fill_coords = np.array([fill_coords[0], fill_coords[1], fill_coords[2]]).T

# Define a nearest neighbour interpolator
interpolator = NearestNDInterpolator(coords, vals)
fill_vals = interpolator(fill_coords)

# Fill the volume with the interpolated values
vol[fill_coords[:, 0], fill_coords[:, 1], fill_coords[:, 2]] = fill_vals
vals = vol.flat[xyz_flat].astype(np.float32)

np.save(file_path.joinpath('depths_um_10.npy'), vals)

# -------------------------------------------------------------------------------------------------
# Downsample to 25 and 50 um atlas
# -------------------------------------------------------------------------------------------------

# Load in 10um results
ba = AllenAtlas(10)

depths = np.load(file_path.joinpath('depths_um_10.npy'))

# unravel so that they are in ixyz
ixyz = np.unravel_index(xyz_flat, ba.image.shape)
# Find the cartesian coordinates of all indices that we have depths for
xyz = ba.bc.i2xyz(np.c_[ixyz[1], ixyz[0], ixyz[2]])

del ba, ixyz

for res in [25, 50]:

    # Transform to 25 um atlas
    ba = AllenAtlas(res)
    df = pd.DataFrame()
    df['depths'] = depths

    # Find the index in the downsampled volume for the cartesian coordinates
    ixyz = ba.bc.xyz2i(xyz)
    df['ind'] = np.ravel_multi_index((ixyz[:, 1], ixyz[:, 0], ixyz[:, 2]), ba.image.shape)

    # Group by index and find the median value
    df_grp = df.groupby('ind')
    df_med = df_grp.depths.median()
    depths_res = df_med.values
    xyz_flat_res = df_med.index.values

    np.save(file_path.joinpath(f'depths_um_{res}.npy'), depths_res.astype(np.float32))

    ibl_vol = np.zeros(np.prod(ba.image.shape))
    ibl_vol[xyz_flat_res] = depths_res
    ibl_vol = ibl_vol.reshape(ba.image.shape)
    ibl_vol[ibl_vol == 0] = np.nan

    ap = -2000 / 1e6
    cmap = 'tab20'
    vmin = 0
    vmax = 2000
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.imshow(np.moveaxis(ba.image[ba.bc.y2i(ap), :, :], 0, 1), extent=ba.extent(axis=1), cmap='gray')
    im = ax.imshow(np.moveaxis(ibl_vol[ba.bc.y2i(ap), :, :], 0, 1),
                   extent=ba.extent(axis=1), cmap=cmap, vmin=vmin, vmax=vmax)
    ba.plot_cslice(ap, volume='boundary', ax=ax)
    cbar = fig.colorbar(im, ax=ax, ticks=np.arange(vmin, vmax, 100))
    cbar.set_ticklabels(np.arange(vmin, vmax, 100))
    cbar.set_label('Depth from surface (um)')
    plt.show()
