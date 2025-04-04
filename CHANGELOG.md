## [0.8.0]
### Added
- `iblatlas.atlas.tilt_spherical` allows to recompute spherical coordinates after tilting the brain
- documented the creation of merfish volumes in ./iblatlas/genomics/merfish_scrapping/02_create_volumes.py
### Fixed
- `iblatlas.gui.atlasview` syntax is now compatible with Python 3.10
- `iblatlas.atlas.xyz_to_depth` allows for duplicate ixyz indices in the volume 

## [0.7.0]
### Added
- `iblatlas.atlas.xyz_to_depth` a lookup from xyz coord to depth in cortex
### Modified
- Moved `iblatlas.plots.get_bc10` to `iblatlas.atlas.get_bc`

## [0.6.0]
### Added
- `iblatlas.genomics.agea` option to load pre-processed volume
- `iblatlas.gui.atlasview` QT-based GUI to explore the Allen Atlas 

## [0.5.7]
### Modified
- `iblatlas.plots.plot_scalar_on_slice` option to add mask with specific color when using vector=True

## [0.5.6]
### Bugfix
- Cast list to numpy array for compatibility with numpy 2.0

## [0.5.5]
### Modified
- removed instances of `np.NaN` to replace with `np.nan` for numpy 2.0 compatibility
- set python 3.10 as minimum python version

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
