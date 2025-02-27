from pathlib import Path
import numpy as np
from iblatlas.gui import atlasview
from iblatlas.genomics import agea, merfish
df_genes, gene_expression_volumes, atlas_agea = agea.load()
av = atlasview.view(atlas=atlas_agea)  # need to have an output argument here or the garbage collector will clean
#av = atlasview.view()  # need to have an output argument here or the garbage collector will clean


from iblatlas.atlas import AllenAtlas
ba = AllenAtlas()
ba.bc.ylim
ba.bc.zlim

# ba.bc.dxyz
# ba.dims2xyz
#
# atlas_agea.bc.dxyz
# atlas_agea.dims2xyz
