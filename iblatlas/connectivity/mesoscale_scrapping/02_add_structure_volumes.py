"""Add per-structure (single-hemisphere) volume columns to the regionalized connectivity
table produced by 01_build_regionalized_matrix.py, so it can be reaggregated to any
parcellation without a second data source.

Volumes are computed from the Allen CCF annotation volume via `iblatlas`
(`AllenAtlas` + `BrainRegions`), by cumulatively counting voxels labelled with each
structure or any of its descendants. The CCF annotation covers both hemispheres and is
bilaterally symmetric by construction, so the single-hemisphere volume used here is half
the whole-volume voxel count.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from iblatlas.atlas import AllenAtlas

OUTPUT_DIR = Path("./connectivity")
RES_UM = 25


def compute_structure_volumes(structure_ids: np.ndarray) -> pd.Series:
    """Compute single-hemisphere volume (mm^3) for each Allen structure id.

    Parameters
    ----------
    structure_ids : numpy.ndarray
        Allen CCF structure ids to compute volumes for.

    Returns
    -------
    pandas.Series
        Volume in mm^3, indexed by `structure_id`, including all descendant
        (e.g. layer) voxels.
    """
    ba = AllenAtlas(res_um=RES_UM)
    br = ba.regions
    br.compute_hierarchy()

    leaf_idx, counts = np.unique(ba.label, return_counts=True)
    voxel_counts = np.zeros(len(br.id), dtype=np.int64)
    for idx, cnt in zip(leaf_idx, counts):
        chain = br.hierarchy[:, idx]
        voxel_counts[chain[chain >= 0]] += cnt

    voxel_volume_mm3 = (RES_UM / 1000) ** 3
    bilateral_volume = voxel_counts * voxel_volume_mm3
    single_hemisphere_volume = bilateral_volume / 2

    id_to_idx = {sid: i for i, sid in enumerate(br.id)}
    volumes = {sid: single_hemisphere_volume[id_to_idx[sid]] for sid in structure_ids if sid in id_to_idx}
    return pd.Series(volumes, name="volume_mm3")


def main() -> None:
    connectivity_path = OUTPUT_DIR.joinpath("allen_mouse_regionalized_connectivity.pqt")
    connectivity = pd.read_parquet(connectivity_path)

    all_ids = np.union1d(connectivity["source_structure_id"].unique(), connectivity["target_structure_id"].unique())
    volumes = compute_structure_volumes(all_ids)

    connectivity["source_volume_mm3"] = connectivity["source_structure_id"].map(volumes)
    connectivity["target_volume_mm3"] = connectivity["target_structure_id"].map(volumes)

    missing = connectivity["source_volume_mm3"].isna().sum() + connectivity["target_volume_mm3"].isna().sum()
    print(f"missing volume values: {missing}")

    connectivity.to_parquet(connectivity_path, index=False)
    print(f"Updated {connectivity_path} with {len(volumes)} structure volumes (single-hemisphere, mm^3).")


if __name__ == "__main__":
    main()
