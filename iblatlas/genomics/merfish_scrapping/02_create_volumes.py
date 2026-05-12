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

Note: This script may require significant computational resources and time to run, especially for the subclasses
"""
import joblib
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
# Convert continuous CCF coordinates to integer voxel indices on the AGEA 200 µm grid
df_cells["ix"] = ba.bc.x2i(df_cells["x"].values).astype(np.int16)
df_cells["iy"] = ba.bc.y2i(df_cells["y"].values).astype(np.int16)
df_cells["iz"] = ba.bc.z2i(df_cells["z"].values).astype(np.int16)


def aggregates_slices(
    df_cells, section_type="coronal", brain_atlas=None, pivot_column="class"
):
    """
    Aggregate per-cell MERFISH counts into sparse 3-D voxel bins, handling coronal
    and sagittal section geometries separately.

    Because MERFISH cells come from discrete physical sections, we first count cells
    within each (section, in-plane-voxel) bin, derive the mean out-of-plane coordinate
    for that bin, convert it to a voxel index, and finally sum across sections that
    land in the same 3-D voxel.  The column `n_slices` records how many sections
    contributed to each voxel so that downstream code can normalise accordingly.

    Parameters
    ----------
    df_cells : pd.DataFrame
        One row per cell.  Must contain columns: donor_label, brain_section_label,
        ix, iy, iz (voxel indices), x, y, z (CCF coordinates), n_cells, and the
        pivot_column.
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
    # For coronal sections the out-of-plane axis is y (dorsoventral); for sagittal it is x (mediolateral).
    # ih: the in-plane voxel index used for grouping within a section
    # h:  the continuous coordinate of the out-of-plane axis, averaged per bin then converted to a voxel index
    match section_type:
        case "coronal":
            (ih, h) = ("ix", "y")
            df_cells_slices = df_cells.loc[i_coronal, :].copy()
        case "sagittal":
            (ih, h) = ("iy", "x")
            df_cells_slices = df_cells.loc[i_sagittal, :].copy()
    # Count cells per (section, in-plane voxel, ap-voxel) bin, pivoted by cell-type label
    df_slices = df_cells_slices.pivot_table(
        index=["brain_section_label", ih, "iz"],
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
        df_cells_slices.loc[:, ["brain_section_label", ih, "iz", h]]
        .groupby(["brain_section_label", ih, "iz"])
        .mean(),
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
        df_slices["ix"] = brain_atlas.bc.x2i(df_slices["x"].values, mode="clip").astype(
            np.int16
        )
    # Sum across sections that map to the same 3-D voxel
    df_vol = df_slices.iloc[:, 1:].groupby(["ix", "iy", "iz"]).sum().reset_index()
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
# Normalise cell counts by the number of sections that overlapped each voxel,
# so voxels covered by both coronal and sagittal sections are not double-counted
df_vol.loc[:, classes] = (
    df_vol.loc[:, classes] / df_vol["n_slices"].values[:, np.newaxis]
)
linear_indices = ba._lookup_inds(df_vol.loc[:, ["ix", "iy", "iz"]].to_numpy())

# %% Now create the volumes via interpolation
file_memmap = ROOT_DIR.joinpath(f"merfish_{LEVEL}.npy")
file_labels = ROOT_DIR.joinpath(f"merfish_{LEVEL}_labels.npy")
np.save(file_labels, classes)
# Create an empty file of the correct size
if file_memmap.exists():
    assert OVERWRITE
    memmap = np.lib.format.open_memmap(
        file_memmap,
        mode="w+",
        dtype=np.float16,
        shape=tuple(np.r_[classes.size, ba.label.shape]),
    )
else:
    memmap = np.lib.format.open_memmap(
        file_memmap,
        mode="w+",
        dtype=np.float16,
        shape=tuple(np.r_[classes.size, ba.label.shape]),
    )


def interpolate_volume(i, iclass):
    """
    Interpolate sparse voxel counts for one cell-type label onto the full bilateral 3-D grid.

    Strategy:
    - MERFISH sections cover only the right hemisphere (positive mediolateral x).
    - All data-point x-coordinates are reflected to positive values so the interpolator
      sees a denser, right-hemisphere-only point cloud.
    - Linear scattered-data interpolation fills the right hemisphere; the left hemisphere
      is then produced by flipping the result along the mediolateral axis.
    """
    # Build a meshgrid covering the right hemisphere only (x from midline to nx)
    xxx, yyy, zzz = np.meshgrid(
        np.arange(ba.bc.x2i(0), ba.bc.nx),
        np.arange(ba.bc.ny),
        np.arange(ba.bc.nz),
        indexing="ij",
    )
    values = df_vol.loc[:, iclass].to_numpy()
    points = df_vol.loc[:, ["ix", "iy", "iz"]].to_numpy().astype(float)
    # Reflect all x-coordinates to the right hemisphere (positive mediolateral)
    # so that left-hemisphere observations also inform the interpolation
    points[:, 0] = ba.bc.x2i(np.abs(ba.bc.i2x(points[:, 0])))

    volume_interp = scipy.interpolate.griddata(
        points, values, (xxx, yyy, zzz), method="linear", fill_value=0, rescale=False
    )
    volume_interp = np.moveaxis(volume_interp, ba.xyz2dims, np.arange(3))

    # Mirror the right hemisphere to obtain the left, then concatenate into a full bilateral volume
    volume_interp_dual = np.concatenate(
        (np.flip(volume_interp, axis=ba.xyz2dims[0]), volume_interp), axis=0
    )
    memmap[i] = volume_interp_dual


# here it was taking too long for the 338 subclasses so I use multiprocessing to speed it up
jobs = []
for i, iclass in enumerate(classes):
    jobs.append(joblib.delayed(interpolate_volume)(i, iclass))
list(tqdm.tqdm(joblib.Parallel(return_as="generator", n_jobs=8)(jobs)))


# %% Eventually copy the results to s3
print(
    f'aws s3 cp --recursive --exclude "*" --include "merfish_*.npy"'
    f' {ROOT_DIR}/ s3://ibl-brain-wide-map-public/atlas/merfish/ --profile ibl --dryrun'
)
# touch /mnt/s1/2025/denoise_agea/2025-04-04.version
# aws s3 cp /mnt/s1/2025/denoise_agea/2025-03-18.version s3://ibl-brain-wide-map-public/atlas/merfish/2025-04-04.version  # NOQA