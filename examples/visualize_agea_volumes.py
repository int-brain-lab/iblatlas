"""
Short example to show how to load the AGEA expression data both raw and processed,
and how to display it using the atlasview GUI
This requires having installed the package with optional dependencies: `pip install iblatlas[gui]`
"""
# %%
import numpy as np

from iblatlas.gui import atlasview
from iblatlas.genomics import agea

_, expression_volumes_raw, _ = agea.load()
df_genes, expression_volumes, atlas_agea = agea.load(label='processed')

av = atlasview.view(atlas=atlas_agea)

# %%
i = 690
levels = np.nanpercentile(expression_volumes[i, :], [0.5, 99.5])
av.set_volume(expression_volumes[i, :], levels=levels)
av.set_volume(expression_volumes_raw[i, :], levels=levels)
