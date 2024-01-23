import matplotlib.pyplot as plt
from iblatlas.genomics import agea

df_genes, gene_expression_volumes, atlas_agea = agea.load()

igenes = (0,)
fig, axs = plt.subplots(3, 2, sharex=True, sharey=True)
atlas_agea.plot_cslice(0, ax=axs[0, 0])
atlas_agea.plot_cslice(0, ax=axs[1, 0], volume='annotation')
atlas_agea.plot_cslice(0, ax=axs[2, 0], volume=gene_expression_volumes[igenes[0]], cmap='viridis')
atlas_agea.plot_sslice(0, ax=axs[0, 1])
atlas_agea.plot_sslice(0, ax=axs[1, 1], volume='annotation')
atlas_agea.plot_sslice(0, ax=axs[2, 1], volume=gene_expression_volumes[igenes[0]], cmap='viridis')
fig.tight_layout()
