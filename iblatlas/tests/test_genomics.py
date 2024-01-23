import unittest
import tempfile
from pathlib import Path
import unittest.mock

import numpy as np
import pandas as pd

from iblatlas.genomics import agea


class TestLoadAgea(unittest.TestCase):

    def test_load(self):
        """
        Tests the mechanics of the load function by faking a 3 experiments gene expression dataset
        Here we just want to make sure the function returns the right dataframes and volumes,
        and that the atlas reads the right label and image volumes
        """
        with tempfile.TemporaryDirectory() as folder_cache:
            df_genes = pd.DataFrame({
                'id': np.array(['74658173', '71247618', '74511936']),
                'gene': ['Ndufa10', 'Grik2', 'Dmp1']
            })
            df_genes.to_parquet(Path(folder_cache).joinpath('gene-expression.pqt'))
            expression_matrices = np.random.rand(3, *agea.DIM_EXP[1:])
            expression_matrices.tofile(Path(folder_cache).joinpath('gene-expression.bin'))
            image = np.random.rand(*agea.DIM_EXP[1:])
            np.save(Path(folder_cache).joinpath('image.npy'), image)
            label = np.random.rand(*agea.DIM_EXP[1:])
            np.save(Path(folder_cache).joinpath('label.npy'), label)
            Path(folder_cache).joinpath('toto.version').touch()
            df_genes, gene_expression_volumes, atlas_agea = agea.load(
                folder_cache=folder_cache, expression_size=(3, 58, 41, 67))
            self.assertEqual(df_genes.shape, (3, 2))
            self.assertEqual(gene_expression_volumes.shape, (3, 58, 41, 67))
            self.assertEqual(atlas_agea.image.shape, (58, 41, 67))
            self.assertEqual(atlas_agea.label.shape, (58, 41, 67))
            # on windows we need to close the memmap for the tempdir to be deleted
            gene_expression_volumes._mmap.close()
