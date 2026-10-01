"""
Lipid Volume Creation

Builds one volume per lipid on the AGEA / MERFISH 200 um grid, following the MERFISH procedure
(`merfish_scrapping/02_create_volumes.py`) adapted to continuous coronal MALDI data:
1. intensities are log-transformed and averaged within each (section, ml, dv) voxel bin, the ap voxel
   index comes from the mean ap coordinate of the bin, then bins from different sections landing in the
   same voxel are averaged
2. sections are full coronal: the left hemisphere is folded onto the right before interpolation
3. linear interpolation on the right hemisphere, mirrored to the left

Subjects:
- the volume uses the 2 densely sampled naive males (ReferenceAtlas, SecondAtlas)
- `reliability`: per-lipid voxel-wise correlation between the volumes of the 2 dense brains
- `sex_log_ratio`, `sex_zscore`: female vs. male effect per lipid from the sparse naive brains
  (Male1-3, Female1-3), computed on Beryl region averages. Pregnant brains are not used.

Output (log intensities, uMAIA-normalised arbitrary units, NaN outside the brain):
- lipid_volumes.npy: (n_lipids, ml, dv, ap) float32 on the grid of `iblatlas.genomics.agea.load_atlas()`
- lipids.pqt: the lipids table with the reliability and sex columns
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.interpolate

from iblatlas.genomics import agea

SAMPLES_DENSE = ['ReferenceAtlas', 'SecondAtlas']
FLOOR = 1e-4  # detection floor of the uMAIA-normalised intensities
path_ibl = Path.home().joinpath('Documents', 'datadisk', 'Data', 'lipid_brain_atlas', 'ibl')

ba = agea.load_atlas()
df_lipids = pd.read_parquet(path_ibl.joinpath('lipids.pqt'))
lipids = df_lipids['lipid'].tolist()
df_pixels = pd.read_parquet(path_ibl.joinpath('lipid_pixels.pqt'))
df_pixels = df_pixels.loc[df_pixels['Condition'] == 'naive']
# the intensities have a 1e-4 floor, the few negative values are imputation artefacts from the source data
df_pixels = pd.concat([np.log(df_pixels[lipids].clip(lower=FLOOR)), df_pixels.drop(columns=lipids)], axis=1)
df_pixels['x'] = np.abs(df_pixels['x'])  # fold the left hemisphere onto the right
df_pixels['ix'] = ba.bc.x2i(df_pixels['x'].values, mode='clip')
df_pixels['iz'] = ba.bc.z2i(df_pixels['z'].values, mode='clip')


def aggregate_voxels(df):
    """Averages pixels within each (section, ml, dv) bin, then bins from all sections within each voxel."""
    df_bins = df.groupby(['SectionID', 'ix', 'iz'])[lipids + ['y']].mean().reset_index()
    df_bins['iy'] = ba.bc.y2i(df_bins['y'].values, mode='clip')
    return df_bins.groupby(['ix', 'iy', 'iz'])[lipids].mean().reset_index()


def interpolate_volume(df_voxels):
    """Linearly interpolates the voxel averages on the right hemisphere and mirrors to the left.

    Returns a (n_lipids, ml, dv, ap) array, NaN outside the brain and outside the sampled range."""
    grid = np.meshgrid(np.arange(ba.bc.x2i(0), ba.bc.nx), np.arange(ba.bc.ny), np.arange(ba.bc.nz), indexing='ij')
    interpolator = scipy.interpolate.LinearNDInterpolator(df_voxels[['ix', 'iy', 'iz']].values, df_voxels[lipids].values)
    volume = np.moveaxis(interpolator(*grid), -1, 0)  # (n_lipids, x, y, z)
    volume = np.moveaxis(volume, np.array(ba.xyz2dims) + 1, np.arange(1, 4))  # (n_lipids, ml, dv, ap)
    volume = np.concatenate((np.flip(volume, axis=1), volume), axis=1)
    volume[:, ba.label == 0] = np.nan
    return volume.astype(np.float32)


# %% main volume from the 2 dense brains, and a per-lipid reliability from their agreement
volume = interpolate_volume(aggregate_voxels(df_pixels.loc[df_pixels['Sample'].isin(SAMPLES_DENSE)]))
vol_a, vol_b = (interpolate_volume(aggregate_voxels(df_pixels.loc[df_pixels['Sample'] == s])) for s in SAMPLES_DENSE)
both = np.all(np.isfinite(vol_a) & np.isfinite(vol_b), axis=0)
df_lipids['reliability'] = [np.corrcoef(a[both], b[both])[0, 1] for a, b in zip(vol_a, vol_b)]

# %% sex effect: female - male log intensity per Beryl region, from the sparse naive brains
df_sex = df_pixels.loc[~df_pixels['Sample'].isin(SAMPLES_DENSE)].copy()
df_sex['beryl'] = ba.regions.remap(df_sex['Allen_id'].values, source_map='Allen', target_map='Beryl')
df_sex = df_sex.loc[~df_sex['beryl'].isin(ba.regions.acronym2id(['void', 'root']))]
df_regions = df_sex.groupby(['beryl', 'Sex', 'Sample'])[lipids].mean()
n_brains = df_regions.groupby(level='beryl').size()
df_regions = df_regions.loc[n_brains.index[n_brains == df_sex['Sample'].nunique()]]  # regions seen in all brains
mean, std = (df_regions.groupby(level=['beryl', 'Sex']).agg(f) for f in ('mean', 'std'))
diff = mean.xs('female', level='Sex') - mean.xs('male', level='Sex')
pooled_std = np.sqrt((std.xs('female', level='Sex') ** 2 + std.xs('male', level='Sex') ** 2) / 2)
df_lipids['sex_log_ratio'] = diff.median().values
df_lipids['sex_zscore'] = (diff / pooled_std).median().values

# %%
np.save(path_ibl.joinpath('lipid_volumes.npy'), volume)
df_lipids.to_parquet(path_ibl.joinpath('lipids.pqt'))
# the derived volumes are distributed under CC-BY 4.0 (approved by the authors), see README.md
print(f'aws s3 cp --recursive --exclude "*" --include "lipid_volumes.npy" --include "lipids.pqt"'
      f' {path_ibl}/ s3://ibl-brain-wide-map-public/atlas/lipids/ --profile ibl --dryrun')
