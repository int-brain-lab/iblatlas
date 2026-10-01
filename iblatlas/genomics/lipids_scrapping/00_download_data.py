"""
Downloads the Lipid Brain Atlas pixel table from Zenodo

Fusar Bassini et al. 2025, "The lipidomic architecture of the mouse brain"
https://www.biorxiv.org/content/10.1101/2025.10.13.682018v1
https://github.com/lamanno-epfl/lipidbrainatlas

MALDI-MSI of 173 lipids on 25 um pixels, coronal sections registered to the Allen CCFv3.
The only file we need is `maindata_2.parquet` (9.2 GB): one row per pixel with lipid intensities,
CCF coordinates (mm), subject metadata and lipizone labels.
"""
import hashlib
from pathlib import Path

import requests
from tqdm import tqdm

ZENODO_URL = 'https://zenodo.org/records/15379565/files'
FILES = {'maindata_2.parquet': '457977888a3e38b571740e2d79b0514d'}  # file name: md5
path_zenodo = Path.home().joinpath('Documents', 'datadisk', 'Data', 'lipid_brain_atlas', 'zenodo')
path_zenodo.mkdir(parents=True, exist_ok=True)

for file_name, md5 in FILES.items():
    file_local = path_zenodo.joinpath(file_name)
    if not file_local.exists():
        with requests.get(f'{ZENODO_URL}/{file_name}?download=1', stream=True) as r:
            r.raise_for_status()
            with open(file_local, 'wb') as f, tqdm(total=int(r.headers['content-length']), unit='B', unit_scale=True) as pbar:
                for chunk in r.iter_content(chunk_size=2 ** 24):
                    pbar.update(f.write(chunk))
    hasher = hashlib.md5()
    with open(file_local, 'rb') as f:
        while chunk := f.read(2 ** 24):
            hasher.update(chunk)
    assert hasher.hexdigest() == md5, f'{file_name} md5 mismatch, delete the file and download again'
