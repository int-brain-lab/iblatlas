---
name: iblatlas-connectivity-mesoscale
description: Load the Allen/Knox regularized-regression mesoscale mouse connectivity matrix and the Harris et al. cortico-thalamic hierarchy scores using iblatlas. Use when the user asks to load mesoscale connectivity, the mouse connectome, connection strength/density, or cortical/thalamic hierarchy scores.
allowed-tools: [Read]
---
name: ea-load-mesoscale-connectivity

# Load mesoscale mouse connectivity + cortico-thalamic hierarchy

```python
from iblatlas.connectivity import mesoscale
df = mesoscale.load()                   # regionalized connectivity matrix
hierarchy = mesoscale.load_hierarchy()  # Harris et al. hierarchy scores
```

## Outputs

| Variable | Description |
|---|---|
| `df` | Tidy long dataframe (671_628, 9): one row per (source structure, hemisphere, target structure, metric) |
| `hierarchy` | Dataframe (122, 9): one row per (area, correction scheme), 61 cortical/thalamic areas |

### `df` columns

- `source_structure_id` / `source_acronym`, `target_structure_id` / `target_acronym`: Allen CCF structure (291 "summary structures")
- `hemisphere`: `'ipsi'` or `'contra'` (relative to the source)
- `metric`: one of `'connection_strength'`, `'connection_density'`, `'normalized_connection_strength'`, `'normalized_connection_density'`
- `value`: the metric value for this (source, hemisphere, target) triple
- `source_volume_mm3` / `target_volume_mm3`: single-hemisphere structure volume (mm³), including descendants (e.g. layers)

### `hierarchy` columns

- `correction`: `'cre_conf'` (Cre-line confidence weighted) or `'no_conf'` (unweighted)
- `area`: acronym — matches `source_acronym` / `target_acronym` in `df` directly, same Allen CCF acronyms
- `region_type`: `'C'` cortex or `'T'` thalamus
- `cc_before` / `cc_iter`, `cctc_before` / `cctc_iter`, `cctcct_before` / `cctcct_iter`: hierarchy score at increasing model stages (cortico-cortical only → + thalamocortical → + corticothalamic), before/after iterative refinement. Lower = more feedforward/earlier (e.g. `VISp` ≈ −0.33 to −0.42, bottom of the visual hierarchy); higher = more association-like. Use `cctcct_iter` with `correction == 'cre_conf'` for the paper's main/full result.

## Key points

- **Only `connection_strength` is additive** across a re-parcellation (it's proportional to integrated axon volume). To reaggregate onto a coarser parcellation (e.g. Cosmos), map `source_structure_id` / `target_structure_id` via `iblatlas.regions.BrainRegions().id2id(ids, mapping='Cosmos')`, then `groupby([...])['value'].sum()`. This is only exact because a proper Allen ontology grouping (Beryl, Cosmos) partitions the 291 summary structures with no overlap or gaps — verify that if using a custom grouping.
- `connection_density` and the two normalized metrics are **ratios** — never sum or average them across a regrouping. Re-derive them from the aggregated `connection_strength` and `source_volume_mm3` / `target_volume_mm3` at the new grouping instead.
- This is anterograde **structural** connectivity (fit via regularized regression across injection experiments to correct for injection-site overlap contamination) — it carries no timing/conduction-delay information, and doesn't fully correct for fibers of passage.
- Empirically, mean outgoing `connection_strength` per area does **not** correlate significantly with the Harris hierarchy score (Spearman ρ≈0.14, p≈0.29, n=57 areas) — raw connectivity strength alone is a weak hierarchy proxy once injection-overlap/normalization artifacts are corrected for.
- If relating connectivity to response latency, prefer the Harris hierarchy score (layer-resolved feedforward/feedback classification) over raw connection strength as the proxy — Siegle et al. (2021, *Nature*, "Survey of spiking in the mouse visual system reveals functional hierarchy") correlated this same hierarchy score against measured visual response latency, a good template for this kind of analysis.
- Colour scatter/heatmap points or tick labels by Allen region: `iblatlas.regions.BrainRegions().get(br.acronym2id(acronyms)).rgb` gives per-region RGB (uint8, one row per acronym).
- See `examples/atlas_connectivity_mesoscale.ipynb` in this repo for a worked example (hierarchy-vs-strength scatter + Cosmos-level reaggregated heatmap, both colour-coded by Allen region).

## References

1. Oh SW et al. (2014) A mesoscale connectome of the mouse brain. *Nature* 508, 207–214. doi: 10.1038/nature13186
2. Knox JE et al. (2018) High-resolution data-driven model of the mouse connectome. *Network Neuroscience* 3, 217–236. doi: 10.1162/netn_a_00066
3. Harris JA et al. (2019) Hierarchical organization of cortical and thalamic connectivity. *Nature* 575, 195–202. doi: 10.1038/s41586-019-1716-z
