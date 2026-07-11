---
name: iblatlas-genomics-merfish
description: Load the MERFISH single-cell atlas (parquet dataframes) and the MERFISH cell-type density volumes using iblatlas, including how to denoise the volumes. Use when the user asks to load MERFISH, single-cell atlas data, cell-type data, MERFISH volumes, or how to denoise/preprocess them.
allowed-tools: [Read]
---
name: ea-load-merfish-dataframes

# Load MERFISH single-cell atlas

Read `sources/examples/05_load_merfish.py` for the canonical pattern and adapt it to context.

## Outputs

| Variable | Description |
|---|---|
| `df_cells` | One row per cell; columns `x`, `y`, `z` (MNI metres), `class`, `subclass`, `supertype`, `cluster` (integer keys) |
| `df_classes` | Cell-class metadata; key column: `class_rgba` (packed int, convert with `merfish.int2rgb()`) |
| `df_subclasses` | Subclass metadata |
| `df_supertypes` | Supertype metadata |
| `df_clusters` | Cluster metadata |
| `df_neurotransmitters` | Neurotransmitter metadata |

## Key points

- Cell coordinates are in **metres** (MNI space); multiply by `1e6` for µm when plotting.
- Get RGB colours for scatter plots: `rgb = merfish.int2rgb(df_classes.loc[df_cells['class'], 'class_rgba'])`.
- Filter a coronal slice: `mask = np.abs(df_cells['y'] - yslice) < (30 / 1e6)`.
- Use `df_cells['class']` as an index into `df_classes` (and similarly for other taxonomy levels).

---

# Load MERFISH cell-type density volumes

```python
from iblatlas.genomics import merfish
volume, labels, agea_atlas = merfish.load_volume(level='class')  # denoised by default; label='' for raw
```

| Variable | Shape | Description |
|---|---|---|
| `volume` | `(n_types, dim0, dim1, dim2)` | Per-voxel cell-type density |
| `labels` | `(n_types,)` | Type index per channel — same key space as `df_classes.index` above |
| `agea_atlas` | — | Atlas matching the volume grid |

`level` is one of `'class'`, `'subclass'`, `'supertype'`, `'cluster'`. The grid is the AGEA 200 µm
grid, not the default 25 µm `AllenAtlas()` — index it via `agea_atlas.bc`. See the
`merfish.load_volume` / `merfish.denoise_volume` docstrings for the `label` argument and denoising
details, and `sources/examples/08_load_merfish_volumes.py` for a worked loading + plotting example.
