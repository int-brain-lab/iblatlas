"""
MERFISH Volume Creation

This module processes MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization) data
to create volumetric representations of gene expression or cell type distributions in the mouse brain.

The script performs the following main operations:
1. Loads MERFISH and Allen Brain Atlas data
2. Aggregates cell data into slices (coronal and sagittal)
3. Combines and normalizes the data
4. Interpolates the data to create full 3D volumes
5. Saves the resulting volumes as memory-mapped numpy arrays

Key variables:
- LEVEL: Determines the level of cell classification ('class', 'subclass', or 'neurotransmitter')
- OVERWRITE: Boolean flag to control whether existing output files should be overwritten

Output:
- Memory-mapped numpy array containing the interpolated volumes
- Numpy array containing the labels for the volumes

Note: all cell types are interpolated from a single triangulation, it takes a few minutes per level
"""
from pathlib import Path
import tqdm

import numpy as np
import pandas as pd
import scipy.interpolate

from iblatlas.atlas import AllenAtlas
from iblatlas.genomics import agea, merfish

# 5322 clusters, 1201 supertypes, 338 subclasses, 34 classes
df_genes, expression_volumes, atlas_agea = agea.load()
ba = atlas_agea
# ba = AllenAtlas()
# av = atlasview.view(atlas=ba)
LEVEL = "subclass"  # 'class'  'neurotransmitter', or 'subclass'
OVERWRITE = False
ROOT_DIR = Path("/mnt/s1/2025/denoise_agea")

(
    df_cells,
    df_classes,
    df_subclasses,
    df_supertypes,
    df_clusters,
    _,
    df_neurotransmitters,
) = merfish.load()
df_cells["n_cells"] = 1.0
# The volumes are built for the right hemisphere and mirrored. The coronal donors cover both hemispheres, so the
# left cells are folded onto the right about the midline (44 um medial of bregma, see `BrainAtlas.midline_x`),
# in continuous coordinates. The hemisphere is kept so that the two halves of a coronal section remain two
# separate observations when binning the sections.
df_cells["hemisphere"] = np.where(df_cells["x"].values < ba.midline_x, -1, 1).astype(np.int8)
df_cells["x"] = ba.midline_x + np.abs(df_cells["x"].values - ba.midline_x)
# First ML voxel whose centre is right of the midline. A folded cell within a few um of the midline may round to
# the voxel left of it, so the ML indices are clipped to the right half.
IX_RIGHT = int(np.ceil(ba.imidline))
# Convert continuous CCF coordinates to integer voxel indices on the AGEA 200 µm grid
df_cells["ix"] = np.maximum(ba.bc.x2i(df_cells["x"].values), IX_RIGHT).astype(np.int16)
df_cells["iy"] = ba.bc.y2i(df_cells["y"].values).astype(np.int16)
df_cells["iz"] = ba.bc.z2i(df_cells["z"].values).astype(np.int16)


def aggregates_slices(
    df_cells, section_type="coronal", brain_atlas=None, pivot_column="class"
):
    """
    Aggregate per-cell MERFISH counts into sparse 3-D voxel bins, handling coronal
    and sagittal section geometries separately.

    Because MERFISH cells come from discrete physical sections, we first count cells
    within each (section, hemisphere, in-plane-voxel) bin, derive the mean out-of-plane coordinate
    for that bin, convert it to a voxel index, and finally sum across sections that
    land in the same 3-D voxel.  The column `n_slices` records how many sections
    contributed to each voxel so that downstream code can normalise accordingly.

    Parameters
    ----------
    df_cells : pd.DataFrame
        One row per cell, folded onto the right hemisphere.  Must contain columns: donor_label,
        brain_section_label, hemisphere, ix, iy, iz (voxel indices), x, y, z (CCF coordinates),
        n_cells, and the pivot_column.
    section_type : {'coronal', 'sagittal'}
        Coronal donors: C57BL6J-638850, Zhuang-ABCA-1, Zhuang-ABCA-2.
        Sagittal donors: Zhuang-ABCA-3, Zhuang-ABCA-4.
    brain_atlas : AllenAtlas, optional
    pivot_column : str
        Hierarchical level to count ('class', 'subclass', 'neurotransmitter').

    Returns
    -------
    df_vol : pd.DataFrame
        Sparse voxel table with columns [ix, iy, iz, <cell-type labels>, n_cells, n_slices].
    """
    brain_atlas = AllenAtlas() if brain_atlas is None else brain_atlas
    print(f"select {section_type} slices")
    # ['C57BL6J-638850',  'Zhuang-ABCA-1',  'Zhuang-ABCA-2',  'Zhuang-ABCA-3', 'Zhuang-ABCA-4']
    # Zhuang 3 and 4 are sagittal, all others coronal
    i_sagittal = df_cells["donor_label"].isin(["Zhuang-ABCA-3", "Zhuang-ABCA-4"])
    i_coronal = ~i_sagittal
    # In iblatlas x is mediolateral, y anteroposterior and z dorsoventral.
    # For coronal sections the out-of-plane axis is y (anteroposterior); for sagittal it is x (mediolateral).
    # ih: the in-plane voxel index used for grouping within a section, along with iz (dorsoventral)
    # h:  the continuous coordinate of the out-of-plane axis, averaged per bin then converted to a voxel index
    match section_type:
        case "coronal":
            (ih, h) = ("ix", "y")
            df_cells_slices = df_cells.loc[i_coronal, :].copy()
        case "sagittal":
            (ih, h) = ("iy", "x")
            df_cells_slices = df_cells.loc[i_sagittal, :].copy()
    # Count cells per (section, hemisphere, in-plane voxel, dv-voxel) bin, pivoted by cell-type label.
    # The hemisphere keeps the two folded halves of a coronal section as two observations.
    bins = ["brain_section_label", "hemisphere", ih, "iz"]
    df_slices = df_cells_slices.pivot_table(
        index=bins,
        values="n_cells",
        columns=pivot_column,
        aggfunc="sum",
        fill_value=0,
    )
    df_slices["n_cells"] = df_slices.sum(axis=1)
    df_slices["n_slices"] = 1
    # Compute the mean out-of-plane continuous coordinate within each bin,
    # which gives a representative position between sections for that voxel
    df_slices = df_slices.merge(
        df_cells_slices.loc[:, bins + [h]].groupby(bins).mean(),
        left_index=True,
        right_index=True,
        how="left",
    ).reset_index()
    # Convert the mean out-of-plane coordinate to a voxel index
    if section_type == "coronal":
        df_slices["iy"] = brain_atlas.bc.y2i(df_slices["y"].values, mode="clip").astype(
            np.int16
        )
    elif section_type == "sagittal":
        ix = brain_atlas.bc.x2i(df_slices["x"].values, mode="clip")
        df_slices["ix"] = np.maximum(ix, IX_RIGHT).astype(np.int16)
    # Sum across sections (and hemispheres) that map to the same 3-D voxel
    df_vol = (
        df_slices.drop(columns=["brain_section_label", "hemisphere", h])
        .groupby(["ix", "iy", "iz"])
        .sum()
        .reset_index()
    )
    return df_vol


# %%
level_dictionary = {
    "class": dict(
        df=df_classes, n=df_classes.shape[0] - 1, index=df_classes.index[1:].values
    ),
    "subclass": dict(
        df=df_subclasses,
        n=df_subclasses.shape[0] - 1,
        index=df_subclasses.index[1:].values,
    ),
    "neurotransmitter": dict(
        df=df_neurotransmitters,
        n=df_neurotransmitters.shape[0],
        index=df_neurotransmitters.index.values,
    ),
}

# Aggregate coronal and sagittal sections independently, then merge
df_vol_sagittal = aggregates_slices(
    df_cells, section_type="sagittal", brain_atlas=ba, pivot_column=LEVEL
)
df_vol_coronal = aggregates_slices(
    df_cells, section_type="coronal", brain_atlas=ba, pivot_column=LEVEL
)
df_vol = (
    pd.concat([df_vol_sagittal, df_vol_coronal], ignore_index=True)
    .groupby(["ix", "iy", "iz"])
    .sum()
    .reset_index()
)
# df_vol = df_vol_coronal
classes = level_dictionary[LEVEL]["index"]
# Normalise cell counts by the number of sections (and hemispheres) that overlapped each voxel,
# so voxels covered by several sections, both hemispheres, or both coronal and sagittal donors
# hold the mean count per section and are not double-counted
counts = list(classes) + ["n_cells"]
df_vol.loc[:, counts] = df_vol.loc[:, counts] / df_vol["n_slices"].values[:, np.newaxis]
linear_indices = ba._lookup_inds(df_vol.loc[:, ["ix", "iy", "iz"]].to_numpy())

# %% Now create the volumes via interpolation
file_memmap = ROOT_DIR.joinpath(f"merfish_{LEVEL}.npy")
file_labels = ROOT_DIR.joinpath(f"merfish_{LEVEL}_labels.npy")
np.save(file_labels, classes)
# Create an empty file of the correct size. The shape must be python ints: with numpy 2 the repr of numpy ints
# (np.int64(...)) would be written in the npy header and the file could not be read back
assert OVERWRITE or not file_memmap.exists(), f"{file_memmap} exists, set OVERWRITE = True"
memmap = np.lib.format.open_memmap(
    file_memmap,
    mode="w+",
    dtype=np.float16,
    shape=(int(classes.size), *map(int, ba.label.shape)),
)


# Interpolate the sparse voxel counts of all cell types onto the full bilateral 3-D grid.
# - The cells were folded onto the right hemisphere about the midline and both hemispheres
#   were pooled per voxel before normalisation, so each voxel appears once in `df_vol`.
#   Note that the Delaunay triangulation keeps a single value for duplicated points: mirroring
#   the points here instead would silently drop one hemisphere.
# - Linear scattered-data interpolation fills the right hemisphere; the left hemisphere
#   is then produced by mirroring the result about the midline.
# - All types share the same points, so the triangulation is computed once and all types are
#   interpolated in a single call (values of shape (n_points, n_types)).
# Grid covering the right hemisphere only (x from the midline to nx)
xxx, yyy, zzz = np.meshgrid(
    np.arange(IX_RIGHT, ba.bc.nx),
    np.arange(ba.bc.ny),
    np.arange(ba.bc.nz),
    indexing="ij",
)
points = df_vol.loc[:, ["ix", "iy", "iz"]].to_numpy().astype(float)
interpolator = scipy.interpolate.LinearNDInterpolator(
    points, df_vol.loc[:, classes].to_numpy().astype(float), fill_value=0, rescale=False
)
volumes_interp = interpolator(xxx, yyy, zzz)  # (nx_right, ny, nz, n_types)
# Mirror the right hemisphere about the midline: the left ML column i takes the right column
# nearest to 2 * imidline - i
ix_mirror = np.round(2 * ba.imidline - np.arange(IX_RIGHT)).astype(int)
ix_mirror = np.clip(ix_mirror, IX_RIGHT, ba.bc.nx - 1) - IX_RIGHT
for i in tqdm.tqdm(range(classes.size)):
    volume_interp = np.moveaxis(volumes_interp[..., i], ba.xyz2dims, np.arange(3))
    volume_left = np.take(volume_interp, ix_mirror, axis=ba.xyz2dims[0])
    memmap[i] = np.concatenate((volume_left, volume_interp), axis=ba.xyz2dims[0])
memmap.flush()


# %% Eventually copy the results to s3
print(
    f'aws s3 cp --recursive --exclude "*" --include "merfish_*.npy"'
    f' {ROOT_DIR}/ s3://ibl-brain-wide-map-public/atlas/merfish/ --profile ibl --dryrun'
)
# touch /mnt/s1/2025/denoise_agea/2025-04-04.version
# aws s3 cp /mnt/s1/2025/denoise_agea/2025-03-18.version s3://ibl-brain-wide-map-public/atlas/merfish/2025-04-04.version  # NOQA