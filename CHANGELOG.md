## [1.0.0]

### Added
- `iblatlas.genomics` module for working with genomics data from Allen contains
  - the Allen gene expression atlas in the `iblatlas.genomics.agea` module
  - the Allen cell types atlas in the `iblatlas.genomics.merfish` module
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
