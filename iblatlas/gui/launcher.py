from pathlib import Path
import numpy as np
from iblatlas.gui import atlasview
from iblatlas.genomics import agea, merfish
df_genes, gene_expression_volumes, atlas_agea = agea.load()

av = atlasview.view(atlas=atlas_agea)  # need to have an output argument here or the garbage collector will clean
#av = atlasview.view()  # need to have an output argument here or the garbage collector will clean

# %%
# gene idx 123 has a very decent quality volume
# gene 4111 is much poorer, with not too many holes
# 2451 lots of missing data
volume = gene_expression_volumes[2451].copy().astype(float)
av.set_volume(volume)


# %%
import pyqtgraph as pg

# %%

