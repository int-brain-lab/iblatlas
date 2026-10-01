"""
Converts the Lipid Brain Atlas pixel table to IBL coordinates

- keeps the 173 lipid intensities (float32) and the subject / section metadata
- converts the CCF coordinates (mm, ap/dv/ml from the volume corner) to IBL xyz (m, ml/ap/dv from bregma)
- drops the knock-out brain (GBA1) and the pixels without CCF coordinates

Samples: 2 densely sampled males (ReferenceAtlas, SecondAtlas, ~30 sections each) and 9 sparsely
sampled brains (Male1-3, Female1-3, Pregnant1,2,4; 4-6 sections each)
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from iblatlas.atlas import AllenAtlas

N_LIPIDS = 173  # the lipids are the first columns of the table
METADATA = ['Sample', 'Sex', 'Condition', 'SectionID', 'xccf', 'yccf', 'zccf', 'lipizone_names']
path_lba = Path.home().joinpath('Documents', 'datadisk', 'Data', 'lipid_brain_atlas')

ba = AllenAtlas()
pf = pq.ParquetFile(path_lba.joinpath('zenodo', 'maindata_2.parquet'))
lipids = pf.schema_arrow.names[:N_LIPIDS]
# reading row groups one by one and casting to float32 halves the memory footprint (7.5M pixels)
df_pixels = pd.concat([
    pf.read_row_group(i, columns=lipids + METADATA).to_pandas().astype({lip: np.float32 for lip in lipids})
    for i in range(pf.num_row_groups)
])
df_pixels = df_pixels.loc[(df_pixels['Condition'] != 'GBA1_disease') & df_pixels['xccf'].notna()]

# xccf: ap, yccf: dv, zccf: ml in mm -> IBL x: ml, y: ap, z: dv in m
xyz = ba.ccf2xyz(df_pixels.loc[:, ['zccf', 'xccf', 'yccf']].values.astype(float) * 1e3)
df_pixels = df_pixels.drop(columns=['xccf', 'yccf', 'zccf']).assign(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2])
df_pixels['Allen_id'] = ba.get_labels(xyz, mode='clip', mapping='Allen')
df_pixels = df_pixels.reset_index(drop=True)

path_lba.joinpath('ibl').mkdir(exist_ok=True)
df_pixels.to_parquet(path_lba.joinpath('ibl', 'lipid_pixels.pqt'))
pd.DataFrame({'lipid': lipids}).to_parquet(path_lba.joinpath('ibl', 'lipids.pqt'))
