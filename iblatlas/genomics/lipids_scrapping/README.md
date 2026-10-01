# Lipid volumes from the Lipid Brain Atlas

Lipid volumes on the AGEA / MERFISH 200 µm grid, from the Lipid Brain Atlas (Fusar Bassini et al. 2025):
- paper: [bioRxiv 2025.10.13.682018](https://www.biorxiv.org/content/10.1101/2025.10.13.682018v1), Nature 2026
- code: [lamanno-epfl/lipidbrainatlas](https://github.com/lamanno-epfl/lipidbrainatlas), [lamanno-epfl/EUCLID](https://github.com/lamanno-epfl/EUCLID)
- data: [Zenodo 15379565](https://zenodo.org/records/15379565)

## Usage

```python
from iblatlas.genomics import lipids
volume, df_lipids, atlas_agea = lipids.load_volume()  # (173, ml, dv, ap) log intensities, memory-mapped
volume, df_lipids, atlas_agea = lipids.load_volume(min_reliability=0.6)  # 158 reproducible lipids, in memory
```

`atlas_agea` is the same 200 µm atlas as the one returned by `agea.load_atlas()`, so the lipid volumes can be indexed like the gene expression and MERFISH volumes.

## Source data

- MALDI mass spectrometry imaging of 173 annotated lipids, 25 µm pixels, coronal sections registered to the Allen CCFv3 with STalign.
- One row per pixel in `maindata_2.parquet` (9.2 GB, 7.5M pixels): lipid intensities, CCF coordinates (`xccf`: ap, `yccf`: dv, `zccf`: ml, mm from the volume corner), subject metadata and lipizone labels.
- Intensities are uMAIA-normalised jointly across all brains: arbitrary units, comparable across sections and brains, with a 1e-4 floor.

| Samples | Sex / condition | Sections | Used for |
|---|---|---|---|
| ReferenceAtlas, SecondAtlas | male, naive | 27 + 35 | volumes and reliability |
| Male1-3, Female1-3 | male / female, naive | 4-6 each | sex effect |
| Pregnant1, 2, 4 | female, pregnant | 6 each | not used |
| GBA1 | male, knock-out | 5 | dropped |

## Scripts

1. `00_download_data.py`: downloads `maindata_2.parquet` from Zenodo and checks the md5.
2. `01_ingest_lba.py`: keeps the lipids (float32) and the subject / section metadata, converts the CCF coordinates to IBL xyz, adds the Allen region ids. Output: `lipid_pixels.pqt` (7.0M pixels, 5 GB), `lipids.pqt`.
3. `02_create_volumes.py`: builds the volumes following `merfish_scrapping/02_create_volumes.py`:
   - log intensities (clipped at the 1e-4 floor) are averaged within each (section, ml, dv) voxel bin. The ap voxel comes from the mean ap coordinate of the bin, then bins from different sections in the same voxel are averaged.
   - the sections are full coronal, so the left hemisphere is folded onto the right before interpolation.
   - linear interpolation of all lipids at once on the right hemisphere, then mirrored to the left.

   Output: `lipid_volumes.npy` (173, 58, 41, 67) float32, NaN outside the brain, and `lipids.pqt` with the columns below.

Runtime on a laptop: about 2 min per script, 8 GB peak memory. Disk: 15 GB in total.

## Lipids table

| Column | Description |
|---|---|
| `lipid` | lipid name, one row per volume channel |
| `reliability` | voxel-wise correlation between the volumes built separately from the 2 dense brains (median 0.74, range 0.35-0.87) |
| `sex_log_ratio` | median over Beryl regions of the female - male log intensity, from the 3 + 3 sparse naive brains |
| `sex_zscore` | same difference divided by the pooled between-brain standard deviation |

## Caveats

- **Coverage**: 91% of the in-brain voxels. The olfactory bulb and the most anterior frontal cortex (about 3 mm) were not sampled and are NaN.
- **Section effects**: lipids with a low `reliability` (e.g. SM 34:1;O2, r = 0.35) show intensity stripes from one section to the next, so filter on `reliability`.
- **Sex**: effects are small, |log ratio| < 0.07 for almost all lipids. Cer 40:2;O2 (+40% in females) and SM 34:1;O2 (+26%, least reliable) stand out. This is based on 3 vs. 3 brains with 4-6 sections each that are not matched, so use it as a flag, not a measurement. The paper's Bayesian male / female model outputs are not on Zenodo.
- **Source artefacts**: about 90 lipids contain a few small negative values, likely from the source's own imputation of lipids missing in some sections. They are clipped at the floor.
- **Multimodal data**: the Zenodo multimodal files do not assign lipids to the Allen MERFISH cells. They give MERFISH labels and genes to the lipid pixels, and lipids to Slide-seq beads (Langlieb et al. 2023). To get lipid values per MERFISH cell, sample `volume` at the cell coordinates.

## License

The derived lipid volumes and lipids table (`lipid_volumes.npy`, `lipids.pqt`) are distributed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/), as approved by the authors of the Lipid Brain Atlas. Their Zenodo records do not set a license. Any use must credit the original work:

- Fusar Bassini L. et al., "The lipidomic architecture of the mouse brain", Nature (2026), [doi:10.1038/s41586-026-11050-0](https://doi.org/10.1038/s41586-026-11050-0)
- source data: [doi:10.5281/zenodo.15379565](https://doi.org/10.5281/zenodo.15379565)

The volumes are a derivative of the source data: they are binned, log-transformed, averaged across 2 brains, interpolated and mirrored, as described above. The code in this folder follows the iblatlas license.
