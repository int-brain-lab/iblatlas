import unittest
import tempfile
from pathlib import Path
import unittest.mock

import numpy as np
import pandas as pd

from iblatlas.genomics import agea
from iblatlas.genomics import merfish


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

    def test_load_atlas(self):
        """`load_atlas()` should build the same atlas as `load()`, from just image.npy/label.npy."""
        with tempfile.TemporaryDirectory() as folder_cache:
            np.save(Path(folder_cache).joinpath('image.npy'), np.random.rand(*agea.DIM_EXP[1:]))
            np.save(Path(folder_cache).joinpath('label.npy'), np.random.rand(*agea.DIM_EXP[1:]))
            atlas_agea = agea.load_atlas(folder_cache=folder_cache)
            self.assertEqual(atlas_agea.image.shape, (58, 41, 67))
            self.assertEqual(atlas_agea.label.shape, (58, 41, 67))


class TestMerfish(unittest.TestCase):

    def test_rgba(self):
        array = np.array([4278216447, 4293910672]).astype(np.uint32)
        assert np.all(np.isclose(
            merfish.int2rgb(array), np.array([[1., 0.4, 0., 1.], [0.56470588, 0.87843137, 0.9372549, 1.]])))
        assert np.all(np.isclose(
            merfish.int2rgb(array, dtype=int), np.array([[255, 102, 0, 255], [144, 224, 239, 255]], dtype=int)))


class TestMerfishDenoiseVolume(unittest.TestCase):

    def test_denoise_volume(self):
        """NaN-fill is restricted to in-brain voxels, and in-brain voxels sum to 1 afterwards."""
        n_types, shape = 5, (4, 4, 4)
        rng = np.random.RandomState(0)
        volume = rng.rand(n_types, *shape).astype(np.float32)
        brain_mask = np.ones(shape, dtype=bool)
        brain_mask[0, 0, 0] = False  # a background voxel, outside the brain
        volume[:, 1, 1, 1] = np.nan  # a voxel with no data at all, inside the brain

        denoised = merfish.denoise_volume(volume, brain_mask, n_drop_non_neuronal=2, sigma=0, seed=1)

        self.assertEqual(denoised.shape, (n_types - 2, *shape))
        np.testing.assert_allclose(denoised.sum(axis=0)[brain_mask], 1.0, atol=1e-5)
        np.testing.assert_array_equal(denoised[:, 0, 0, 0], 0.0)  # out-of-brain stays at 0
        self.assertTrue(np.all(denoised[:, 1, 1, 1] >= 0))  # NaN got a valid Dirichlet fill
        self.assertAlmostEqual(float(denoised[:, 1, 1, 1].sum()), 1.0, places=5)

    def test_denoise_volume_no_smoothing_no_drop(self):
        """With sigma=0, n_drop_non_neuronal=0 and no NaNs, this is just a per-voxel renormalize."""
        n_types, shape = 3, (2, 2, 2)
        volume = np.ones((n_types, *shape), dtype=np.float32)
        brain_mask = np.ones(shape, dtype=bool)
        denoised = merfish.denoise_volume(volume, brain_mask, n_drop_non_neuronal=0, sigma=0)
        self.assertEqual(denoised.shape, (n_types, *shape))
        np.testing.assert_allclose(denoised, 1.0 / n_types, atol=1e-6)


class TestMerfishLoadVolume(unittest.TestCase):

    def test_load_volume(self):
        with tempfile.TemporaryDirectory() as folder_cache:
            n_types, shape = 8, (4, 5, 6)
            volume = np.random.RandomState(2).rand(n_types, *shape).astype(np.float16)
            labels = np.arange(1, n_types + 1)
            np.save(Path(folder_cache, 'merfish_class.npy'), volume)
            np.save(Path(folder_cache, 'merfish_class_labels.npy'), labels)
            fake_atlas = unittest.mock.MagicMock(label=np.ones(shape))

            with unittest.mock.patch('iblatlas.genomics.agea.load_atlas', return_value=fake_atlas):
                vol_raw, lab_raw, ba = merfish.load_volume(level='class', label='', folder_cache=folder_cache)
                np.testing.assert_array_equal(np.asarray(vol_raw), volume)
                np.testing.assert_array_equal(lab_raw, labels)
                self.assertIs(ba, fake_atlas)
                # on windows we need to close the memmap for the tempdir to be deleted
                vol_raw._mmap.close()

                vol_proc, lab_proc, _ = merfish.load_volume(
                    level='class', label='processed', folder_cache=folder_cache)
                self.assertEqual(vol_proc.shape, (n_types - 5, *shape))  # default n_drop_non_neuronal=5
                np.testing.assert_array_equal(lab_proc, labels[:n_types - 5])

    def test_load_volume_bad_args(self):
        with self.assertRaises(AssertionError):
            merfish.load_volume(level='not_a_level')
        with self.assertRaises(AssertionError):
            merfish.load_volume(label='not_a_label')
