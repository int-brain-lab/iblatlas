---
name: iblatlas-genomics-agea
description: Load the Allen Gene Expression Atlas (AGEA) volumes using iblatlas. Use when the user asks to load AGEA, gene expression volumes, or the gene expression atlas.
allowed-tools: [Read]
---
name: ea-load-gene-expression-volumes

# Load AGEA (Allen Gene Expression Atlas)

Read `sources/examples/02_download_gene_expression_atlas.py` for the canonical pattern and adapt it to context.

## Outputs

| Variable | Description |
|---|---|
| `df_genes` | DataFrame with one row per gene; key column: `gene` (gene symbol string) |
| `expression_volumes` | Array of 3-D volumetric gene expression data, indexed by row position in `df_genes` |
| `agea_atlas` | `AllenAtlas`-like object with `.bc` coordinate converter and `.plot_cslice` / `.plot_sslice` methods |

## Key points

- `label='processed'` loads the denoised expression volumes; omit for raw.
- Select a gene volume by index: `igene = np.where(df_genes['gene'] == 'Slc17a7')[0][0]`, then `expression_volumes[igene]`.
- Coordinates via `agea_atlas.bc`: `.x2i()`, `.y2i()`, `.z2i()` convert MNI metres → voxel indices.
- Pass a gene volume directly to `agea_atlas.plot_cslice(y, volume=expression_volumes[igene], cmap='Blues')`.
