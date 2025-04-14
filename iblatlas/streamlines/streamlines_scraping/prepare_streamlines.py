"""
This script was used to generate the following files
- paths_10.pqt
- paths_25.pqt
- paths_50.pqt
- dorsal_flatmap.npy

These files are used by the module iblatlas.streamlines to project 3D volumes or points onto the dorsal cortex
flatmap using the streamline paths

The original paths file that we work from (dorsal_flatmap_paths_10.h5) is provided by Allen and available to download
from here  https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/

Compared to Allen the paths are stored in a pandas Dataframes that has columns
- paths - indicates the streamline index
- lookup - indicates the voxel in the 3D Allen volume for each path

Compared to the streamlines provided by Allen we remove duplicate voxels per streamline, e.g compare
paths[677974] from Allen to paths[677974] from IBL

We have also subsampled the streamlines to the 25 um and 50 um volumes
"""

# -------------------------------------------------------------------------------------------------
# Process for 10 um
# -------------------------------------------------------------------------------------------------

import h5py
import pandas as pd
import numpy as np
from iblatlas.atlas import AllenAtlas, get_bc
from iblatlas.streamlines.utils import project_volume_onto_flatmap
from pathlib import Path
file_path = Path('/Users/admin/Downloads/flatmap')

ba = AllenAtlas(10)
ba_shape = ba.image.shape

del ba

f1 = h5py.File(file_path.joinpath('dorsal_flatmap_paths_10.h5'), 'r+')
paths = f1['paths'][:]
f1.close()

# Paths from Allen are in order AP, DV, ML
paths = np.unravel_index(paths, (1320, 800, 1140))
# Reorder so that the indices are in IBL order AP, ML, DV
paths = np.ravel_multi_index((paths[0], paths[2], paths[1]), ba_shape)
path_shape = paths.shape

# When the value = 0 in a path it is locations where the streamline has ended
# Find the longest streamline

path_df = pd.DataFrame()
paths_idx = np.repeat(np.arange(paths.shape[0]), paths.shape[1])
path_df['lookup'] = paths.flatten().astype(np.int32)
path_df['paths'] = paths_idx.astype(np.int32)

# Remove the rows where the lookup = 0. These indicate locations where there is no voxel for that depth
# e.g shorter streamlines
kp_idx = paths.flat != 0
path_df = path_df[kp_idx]
# Remove duplicate voxels within the same streamline
path_df = path_df.drop_duplicates().reset_index(drop=True)

path_df.to_parquet(file_path.joinpath('paths_10.pqt'))


# -------------------------------------------------------------------------------------------------
# Downsample for 25 um and 50 um
# -------------------------------------------------------------------------------------------------

res = 50  # 25 or 50

f1 = h5py.File(file_path.joinpath('dorsal_flatmap_paths_10.h5'), 'r+')
paths = f1['paths'][:]
f1.close()

bc = get_bc(10)
bc_res = get_bc(res)

# Paths are in order AP, DV, ML
path_shape = paths.shape
paths = paths.ravel()
paths = np.unravel_index(paths, (1320, 800, 1140))

xyz = bc.i2xyz(np.array([paths[2], paths[0], paths[1]]).T)
del paths
ixyz = bc_res.xyz2i(xyz)
del xyz
ixyz = np.ravel_multi_index((ixyz[:, 1], ixyz[:, 0], ixyz[:, 2]), (bc_res.ny, bc_res.nx, bc_res.nz))

paths = ixyz.reshape(path_shape)

del ixyz

# When the value = 0 in a path it is locations where the streamline has ended
# Find the longest streamline
path_df = pd.DataFrame()
paths_idx = np.repeat(np.arange(paths.shape[0]), paths.shape[1])
path_df['lookup'] = paths.flatten().astype(np.int32)
path_df['paths'] = paths_idx.astype(np.int32)

# Remove the rows where the lookup = 0. These indicate locations where there is no voxel for that depth
# e.g shorter streamlines
kp_idx = paths.flat != 0
path_df = path_df[kp_idx]
# Remove duplicate voxels within the same streamline
path_df = path_df.drop_duplicates()

path_df = path_df.drop_duplicates().reset_index(drop=True)

path_df.to_parquet(file_path.joinpath(f'paths_{res}.pqt'))


# -------------------------------------------------------------------------------------------------
# Extract the flatmap file
# -------------------------------------------------------------------------------------------------

f1 = h5py.File(file_path.joinpath('dorsal_flatmap_paths_10.h5'), 'r+')
flatmap = f1['flat'][:]
f1.close()

np.save(file_path.joinpath('dorsal_flatmap.npy'), flatmap)


# -------------------------------------------------------------------------------------------------
# Get the image and annotation files
# -------------------------------------------------------------------------------------------------
ba = AllenAtlas(10)

# Annotation
img = project_volume_onto_flatmap(ba.label, res_um=10, aggr='first', plot=False)
img[0, 0] = 0
np.save(file_path.joinpath('dorsal_annotation.npy'), img.astype(np.int16))

# Image
img = project_volume_onto_flatmap(ba.image, res_um=10, aggr='first', plot=False)
img[0, 0] = 0
np.save(file_path.joinpath('dorsal_image.npy'), img.astype(np.float16))
