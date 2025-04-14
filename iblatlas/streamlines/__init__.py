"""A package for working with Allen streamlines data.

This package provides a way to use the Allen streamline data to
1. Do a lookup of xyz position to depth from the surface of the cortex
2. Project a volume or xyz positions along the streamlines to display them on a flattened dorsal cortex flatmap

Streamlines are the paths that connect the surface of the isocortex to the white matter surface while following the
curvature of these surfaces.

For the xyz to depth position lookup the following files have been generated with 10, 25 and 50 um resolution

     depths/
     ├── depths_ind_{res}.npy
     ├── depths_per_{res}.npy
     ├── depths_um_{res}.npy

- The depths_ind_{res}.npy contains the flattened index for all voxels contained within the isocortex of the Allen
volume at the specified resolution. The ordering of the flattened index is AP, DV, ML
- The depths_per_{res}.npy contains the depth in percentage from the surface of the cortex of each voxel in the isocortex.
The percentage values were taken from the laplacian provided by the Allen institute. Voxels at the surface have a percentage
 depth of 0 while those at the white matter boundary have a percentage depth of 1.
- The depths_um_{res}.npy contains the depth in um from the surface of the cortex of each voxel in the isocortex. The depth
values have been computed using the streamlines provided by the Allen institute. For each streamline the length of
the streamline was computed and voxels assigned a depth based on the steps along the streamline. When multiple streamlines
pass through a given voxel the median depth has been taken.

The 25 and 50 um resolution datasets have been computed by subsampling the 10 um datasets.


To enable the projection onto the dorsal flatmap the following files have been generated with 10, 25 and 50 um resolution

     depths/
     ├── paths_{res}.pqt
     ├── dorsal_flatmap.npy
     ├── dorsal_annotation.npy
     ├── dorsal_image.npy


- The paths_{res}.pqt files contain information about the voxels that correspond to each streamline. They have been
derived from the dorsal_flatmap_paths_10.h5 paths provided by the Allen institute. We store the values as a pqt file and remove
values that correspond to no streamline (paths=0) and duplicate voxels per streamline. As we have removed duplicate
voxels per streamline this means that when taking the average or std along the streamline the behavior will differ
from the Allen code.
- The dorsal_flatmap.npy is a 2D numpy array with the dorsal cortex flatmap. The value in each voxel corresponds the
streamline index that corresponds to that voxel.
- The dorsal_annotation is 2D numpy array with the dorsal cortex flatmap. The value in each voxel corresponds to the
region id on the surface of the cortex that corresponds to that voxel
- The dorsal_image is a 2D numpy array with the dorsal cortex flatmap. The value in each voxel corresponds to the Allen
image volume projected onto the flatmap display

[1] Wang, Quanxin et al., “The Allen Mouse Brain Common Coordinate Framework: A 3D Reference Atlas”
 Cell, Volume 181, Issue 4, 936 - 953.e20 May. 2020, https://doi.org/10.1016/j.cell.2020.04.007

"""