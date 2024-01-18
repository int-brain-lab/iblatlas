"""
Module to download and format gene expression atlas data from Allen Institute AGEA project
https://mouse.brain-map.org/agea
https://www.nature.com/articles/nature05453
http://help.brain-map.org/display/mousebrain/Documentation

volumes for all genes (size 51533 voxels x 4376 genes) ~ .84 Gb float32
correlation volumes (size 51533 voxels x 51533 genes) ~ 9.89 Gb float32

They are 200 micron voxel

ElementSpacing = 200 200 200
DimSize = 67 41 58 (159_326 voxels)

Which means you need to query:
X = 0 to (66*200)
Y = 0 to (40*200)
Z = 0 to (57*200)
"""
# %%
from pathlib import Path
import numpy as np
import pandas as pd


folder_download = Path("/datadisk/gdrive/2022/08_gene_expression_atlas")
file_parquet = folder_download.joinpath('experiments.pqt')


def make_experiments_parquet(file_parquet):
    """
    Query experiments and gene acronyms for the xml API and write dataframe
    with ids and corresponding acronyms
    :param file_parquet:
    :return:
    """
    import requests
    from xml.etree import ElementTree
    from copy import copy

    http_query = ("http://api.brain-map.org/api/v2/data/query.xml?criteria=model::SectionDataSet,rma::"
                  "criteria,[failed$eqfalse],products[id$eq1],plane_of_section[name$eqcoronal],"
                  "treatments[name$eqISH],rma::include,genes,rma::options[num_rows$eqall]")
    response = requests.get(http_query)
    tree = ElementTree.fromstring(response.text)

    def dictify(r, root=True, tags=None):
        if root:
            return {r.tag: dictify(r, False)}
        d = copy(r.attrib)
        if r.text and r.text.strip() != '':
            d = r.text.strip()
        for x in r.findall("./*"):
            if x.tag not in d:
                d[x.tag] = []
            d[x.tag].append(dictify(x, False))
        return d

    datasets = dictify(tree)['Response']['section-data-sets'][0]['section-data-set']
    ids = np.zeros(len(datasets), dtype=np.int64)
    genes = []
    for i, ds in enumerate(datasets):
        ids[i] = ds['id'][0]
        genes.append(ds['genes'][0]['gene'][0]['acronym'][0])

    df_genes = pd.DataFrame(data=np.c_[ids, np.array(genes)], columns=['id', 'gene'])
    df_genes.to_parquet(file_parquet)


if not file_parquet.exists():
    make_experiments_parquet(file_parquet)
df_experiments = pd.read_parquet(file_parquet)
