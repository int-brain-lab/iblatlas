## [0.5.4]
### Modified
- `iblatlas.plots.plot_swanson_vector` option to show colorbar

## [0.5.3]
### Modified
- `iblatlas.plots` uses matplotlib.colormaps.get_cmap instead of deprecated cm.get_cmap

## [0.5.2]
### Modified
- `iblatlas.plots` uses colors.get_cmap instead of deprecated colors.get_cmap

## [0.5.1]
### Modified
- `iblatlas.genomics.merfish.load` does not return an error when dask is installed but
dask.dataframe is not.

## [0.5.0]
### Added
- `iblatlas.genomics.merfish` module for working with the Allen gene expression
 atlas in the `iblatlas.genomics.merfish` module

## [0.4.0]
### Added
- `iblatlas.genomics.agea` module for working with the Allen gene expression
 atlas in the `iblatlas.genomics.agea` module
### Modified
- slices of the atlas are now always returned with consistent sizes regardless of the volume layout on disk
- atlases now can have an extra dimension in the image volume, to allow for multiple layers

## [0.3.0]
### Modified
- Insertion._get_surface_intersection: can return None if and when mode is not set to raise and insertion is not intersecting with any surface
- BrainAtlas.mask(): method that computes a mask of the convex brain volume shape

## [0.2.1]

### Modified
- HOTFIX: Pass reg id into polygon functions in plot_scalar_on_slice

## [0.2.0]

### Modified
- Changed packaging from setup.py to pyproject.toml

## [0.1.1]

### Modified
- Change examples and docstrings to reference iblatlas instead of ibllib

## [0.1.0]

### Added
 - Atlas files that were located in the ibllib github repository. These include all files 
   in the ibllib.atlas module, in examples.atlas and ibllib.tests.test_atlas.py 
