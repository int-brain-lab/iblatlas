"""A package for working with mesoscale mouse connectivity: a regionalized structural
connectivity matrix and a precomputed cortico-thalamic hierarchy.

Mesoscale connectivity matrix
------------------------------

This package provides a way to load the Allen Institute's regularized-regression mesoscale
connectivity matrix, regionalized onto 291 Allen "summary structures", together with
per-structure volumes so the matrix can be reaggregated onto any coarser parcellation
(e.g. Beryl or Cosmos, via `iblatlas.regions.BrainRegions.id2id`).

    connectivity/
    ├── allen_mouse_regionalized_connectivity.pqt
    └── allen_mouse_hierarchy_scores.pqt

-   allen_mouse_regionalized_connectivity.pqt is a tidy long table (671_628, 9), one row per
(source structure, hemisphere, target structure, metric): `connection_strength`,
`connection_density`, `normalized_connection_strength` or `normalized_connection_density`.
Only `connection_strength` is additive across a re-parcellation -- see
`iblatlas.connectivity.mesoscale.load` for details.
-   allen_mouse_hierarchy_scores.pqt is a dataframe (122, 9) of precomputed hierarchy scores for
61 cortical and thalamic areas, under two Cre-line-confidence correction schemes.

See the building scripts in ./connectivity/mesoscale_scrapping/.

[1] S. W. Oh et al., "A mesoscale connectome of the mouse brain," Nature, vol. 508, no. 7495,
 Art. no. 7495, Apr. 2014, doi: 10.1038/nature13186.
[2] J. E. Knox et al., "High-resolution data-driven model of the mouse connectome," Network
 Neuroscience, vol. 3, no. 1, pp. 217-236, Nov. 2018, doi: 10.1162/netn_a_00066.
[3] J. A. Harris et al., "Hierarchical organization of cortical and thalamic connectivity,"
 Nature, vol. 575, no. 7781, Art. no. 7781, Nov. 2019, doi: 10.1038/s41586-019-1716-z.

"""
