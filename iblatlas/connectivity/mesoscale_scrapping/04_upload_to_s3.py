"""Upload the built connectivity parquet files to the public S3 bucket.

Run manually with the `ibl` AWS CLI profile configured. `mesoscale.load()` / `load_hierarchy()`
use the same `<date>.version` sentinel convention as `iblatlas.genomics.agea.load` /
`iblatlas.genomics.merfish.load` (see `mesoscale._sync`): whenever these files are updated,
bump `mesoscale.CURRENT_VERSION`, append the previous value to `mesoscale.OLD_VERSIONS`, and
touch/upload a new `<CURRENT_VERSION>.version` marker here.

# touch ./connectivity/2026-07-11.version
# aws s3 cp ./connectivity/allen_mouse_regionalized_connectivity.pqt \
#     s3://ibl-brain-wide-map-public/atlas/connectivity/allen_mouse_regionalized_connectivity.pqt \
#     --profile ibl
# aws s3 cp ./connectivity/allen_mouse_hierarchy_scores.pqt \
#     s3://ibl-brain-wide-map-public/atlas/connectivity/allen_mouse_hierarchy_scores.pqt \
#     --profile ibl
# aws s3 cp ./connectivity/2026-07-11.version \
#     s3://ibl-brain-wide-map-public/atlas/connectivity/2026-07-11.version \
#     --profile ibl
"""
