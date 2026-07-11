"""Build the tidy regionalized connectivity matrix.

Source: precomputed CSVs served by the Allen Institute
(`download.alleninstitute.org/publications/A_high_resolution_data-driven_model_of_the_mouse_connectome/`),
fetched via `mcmodels.core.VoxelModelCache`. These regionalize the Oh et al. (2014) / Knox et
al. (2018) voxel-scale regularized regression model -- fit across all injection experiments to
correct for injection-site overlap contamination -- onto 291 summary structures, for four
related metrics: `connection_strength`, `connection_density`,
`normalized_connection_strength`, `normalized_connection_density`.

`mcmodels` (github.com/AllenInstitute/mouse_connectivity_models) is an old, effectively
unmaintained package: it requires `allensdk` (which pins old numpy/pandas/scipy) and references
a private scikit-learn function removed in modern releases
(`sklearn.model_selection._search._check_param_grid`). Building this requires a dedicated
throwaway environment and a one-line patch to a local clone -- not something to add as a runtime
dependency of `iblatlas`. Run once to produce the shipped parquet; not part of the installed
package (see `pyproject.toml` excludes).

    uv venv .venv-mcmodels && uv pip install --python .venv-mcmodels/bin/python \
        allensdk six "scikit-learn>=1.3,<1.4" pyarrow
    git clone https://github.com/AllenInstitute/mouse_connectivity_models.git
    # patch mouse_connectivity_models/mcmodels/regressors/nonparametric/nadaraya_watson.py:
    #   wrap `from sklearn.model_selection._search import _check_param_grid` in a try/except
    #   ImportError, falling back to a no-op (only used by a grid-search fit path we never call)
    uv pip install --python .venv-mcmodels/bin/python --no-deps --no-build-isolation \
        -e mouse_connectivity_models
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from mcmodels.core import VoxelModelCache

CACHE_DIR = Path("./mcmodels_cache")  # VoxelModelCache download/cache location
OUTPUT_DIR = Path("./connectivity")

METRICS = {
    "connection_strength": "connection_strength.pkl",
    "connection_density": "connection_density.pkl",
    "normalized_connection_strength": "normalized_connection_strength.pkl",
    "normalized_connection_density": "normalized_connection_density.pkl",
}


def melt_matrix(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Melt a (source_id) x (hemisphere, target_id) matrix into a tidy frame.

    Parameters
    ----------
    df : pandas.DataFrame
        Index is source `structure_id` (int); columns are a
        `(hemisphere, target_id)` MultiIndex.
    metric : str
        Name of the metric this matrix represents.

    Returns
    -------
    pandas.DataFrame
        Columns: `source_structure_id`, `hemisphere`, `target_structure_id`,
        `metric`, `value`.
    """
    stacked = df.stack(level=[0, 1])
    stacked.index.set_names(["source_structure_id", "hemisphere", "target_structure_id"], inplace=True)
    tidy = stacked.rename("value").reset_index()
    tidy["target_structure_id"] = tidy["target_structure_id"].astype(int)
    tidy["metric"] = metric
    return tidy


def main() -> None:
    cache = VoxelModelCache(manifest_file=str(CACHE_DIR.joinpath("voxel_model_manifest.json")))
    structure_tree = cache.get_structure_tree()
    id_to_acronym = {
        s["id"]: s["acronym"] for s in structure_tree.get_structures_by_set_id([cache.SUMMARY_STRUCTURE_SET_ID])
    }

    frames = []
    for metric, filename in METRICS.items():
        df = pd.read_pickle(CACHE_DIR.joinpath(filename))
        frames.append(melt_matrix(df, metric))

    tidy = pd.concat(frames, ignore_index=True)
    tidy["source_acronym"] = tidy["source_structure_id"].map(id_to_acronym)
    tidy["target_acronym"] = tidy["target_structure_id"].map(id_to_acronym)

    tidy.to_parquet(OUTPUT_DIR.joinpath("allen_mouse_regionalized_connectivity.pqt"), index=False)
    print(
        f"Wrote {len(tidy)} rows, {tidy['source_structure_id'].nunique()} source structures, "
        f"metrics={sorted(tidy['metric'].unique())}"
    )


if __name__ == "__main__":
    main()
