"""Extract Harris et al. (2019) hierarchy scores into a tidy parquet table.

Source: precomputed Results tables from the AllenInstitute/MouseBrainHierarchy repository
(companion code to Harris, J.A. et al. "Hierarchical organization of cortical and thalamic
connectivity" Nature 575, 195-202 (2019)). The repo contains hierarchy scores under two
correction schemes -- with Cre-line confidence weighting ("CreConf") and without
("NoCreConf") -- for cortical (C) and thalamic (T) areas, at three levels of connectivity
considered (cortico-cortical only; + thalamo-cortical; + cortico-thalamic), before and after
the iterative self-consistency refinement described in the paper.

    git clone https://github.com/AllenInstitute/MouseBrainHierarchy.git
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_RESULTS = Path("./MouseBrainHierarchy/Results")
OUTPUT_DIR = Path("./connectivity")

COLUMN_RENAME = {
    "areas": "area",
    "CortexThalamus": "region_type",
    "CC before iteration": "cc_before",
    "CC iterated": "cc_iter",
    "CC+TC before iteration": "cctc_before",
    "CC+TC iterated": "cctc_iter",
    "CC+TC+CT before iteration": "cctcct_before",
    "CC+TC+CT iterated": "cctcct_iter",
}


def main() -> None:
    frames = []
    for correction, filename in [
        ("cre_conf", "hierarchy_summary_CreConf.xlsx"),
        ("no_conf", "hierarchy_summary_NoCreConf.xlsx"),
    ]:
        df = pd.read_excel(REPO_RESULTS.joinpath(filename), sheet_name="hierarchy_all_regions")
        df = df.rename(columns=COLUMN_RENAME)
        df.insert(0, "correction", correction)
        frames.append(df)

    hierarchy = pd.concat(frames, ignore_index=True)
    hierarchy.to_parquet(OUTPUT_DIR.joinpath("allen_mouse_hierarchy_scores.pqt"), index=False)
    print(f"Wrote {len(hierarchy)} rows ({hierarchy['area'].nunique()} unique areas).")


if __name__ == "__main__":
    main()
