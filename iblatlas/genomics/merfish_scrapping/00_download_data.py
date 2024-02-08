"""
Downloads the data from the Allen Brain Cell Atlas

Some useful resources as of 2024-01-31:
# https://alleninstitute.github.io/abc_atlas_access/notebooks/zhuang_merfish_tutorial.html
# https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html

# https://allen-brain-cell-atlas.s3.amazonaws.com/index.html
# https://ibl-brain-wide-map-public.s3.amazonaws.com/index.html

# the explorer setup is described by AWS here: https://github.com/awslabs/aws-js-s3-explorer/tree/master
"""
from pathlib import Path
import os
import json
import requests
from one.remote import aws

version = '20231215'
version = '20230830'
download_base = '/datadisk/Data/merfish_atlas/cache'

use_local_cache = False
manifest_path = 'releases/%s/manifest.json' % version

if not use_local_cache:
    url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/' + manifest_path
    manifest = json.loads(requests.get(url).text)
else:
    file = os.path.join(download_base, manifest_path)
    with open(file, 'rb') as f:
        manifest = json.load(f)

s3_allen, bucket_name = aws.get_s3_allen()
for r in manifest['directory_listing']:
    r_dict = manifest['directory_listing'][r]
    for d in r_dict['directories']:
        if d != 'metadata':
            continue
        d_dict = r_dict['directories'][d]
        local_path = Path(download_base).joinpath(d_dict['relative_path'])
        print(local_path)
        # !aws s3 ls s3://allen-brain-cell-atlas/metadata/Zhuang-ABCA-1/20231215/
        aws.s3_download_folder(d_dict['relative_path'], local_path, s3_allen, bucket_name)
