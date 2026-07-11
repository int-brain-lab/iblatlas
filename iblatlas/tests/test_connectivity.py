import unittest
import unittest.mock
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from iblatlas.connectivity import mesoscale


class TestMesoscale(unittest.TestCase):

    def test_load(self):
        """Tests the mechanics of load() by faking a small connectivity table."""
        with tempfile.TemporaryDirectory() as folder_cache:
            df = pd.DataFrame({
                'source_structure_id': [315, 315, 993],
                'source_acronym': ['Isocortex', 'Isocortex', 'MOs'],
                'target_structure_id': [993, 315, 315],
                'target_acronym': ['MOs', 'Isocortex', 'Isocortex'],
                'hemisphere': ['ipsi', 'ipsi', 'ipsi'],
                'metric': ['connection_strength'] * 3,
                'value': [1.0, 2.0, 3.0],
                'source_volume_mm3': [61.5, 61.5, 3.3],
                'target_volume_mm3': [3.3, 61.5, 61.5],
            })
            df.to_parquet(Path(folder_cache).joinpath('allen_mouse_regionalized_connectivity.pqt'))
            Path(folder_cache).joinpath(f'{mesoscale.CURRENT_VERSION}.version').touch()
            out = mesoscale.load(folder_cache=folder_cache)
            pd.testing.assert_frame_equal(out, df)

    def test_load_hierarchy(self):
        """Tests the mechanics of load_hierarchy() by faking a small hierarchy table."""
        with tempfile.TemporaryDirectory() as folder_cache:
            df = pd.DataFrame({
                'correction': ['cre_conf', 'no_conf'],
                'area': ['VISp', 'VISp'],
                'region_type': ['C', 'C'],
                'cctcct_iter': [-0.42, -0.38],
            })
            df.to_parquet(Path(folder_cache).joinpath('allen_mouse_hierarchy_scores.pqt'))
            Path(folder_cache).joinpath(f'{mesoscale.CURRENT_VERSION}.version').touch()
            out = mesoscale.load_hierarchy(folder_cache=folder_cache)
            pd.testing.assert_frame_equal(out, df)

    def test_sync_downloads_when_version_missing_or_old(self):
        """`_sync` should (re)download both files when there is no version flag, or when the
        flag matches a known-old version, and should leave an up-to-date cache alone."""
        with tempfile.TemporaryDirectory() as folder_cache:
            folder_cache = Path(folder_cache)
            with unittest.mock.patch('iblatlas.connectivity.mesoscale.aws.s3_download_file') as mock_dl:
                mesoscale._sync(folder_cache)  # no version flag yet -> downloads
                self.assertEqual(mock_dl.call_count, len(mesoscale.FILENAMES))
                self.assertTrue(folder_cache.joinpath(f'{mesoscale.CURRENT_VERSION}.version').exists())

            with unittest.mock.patch('iblatlas.connectivity.mesoscale.aws.s3_download_file') as mock_dl:
                mesoscale._sync(folder_cache)  # up-to-date flag present -> no download
                mock_dl.assert_not_called()

            folder_cache.joinpath(f'{mesoscale.CURRENT_VERSION}.version').unlink()
            folder_cache.joinpath('2020-01-01.version').touch()
            with unittest.mock.patch(
                'iblatlas.connectivity.mesoscale.OLD_VERSIONS', ['2020-01-01']
            ), unittest.mock.patch('iblatlas.connectivity.mesoscale.aws.s3_download_file') as mock_dl:
                mesoscale._sync(folder_cache)  # flag is a known-old version -> downloads again
                self.assertEqual(mock_dl.call_count, len(mesoscale.FILENAMES))

    def test_reaggregation_example(self):
        """The docstring's reaggregation recipe: connection_strength sums exactly across a
        parcellation that partitions source/target ids with no overlap or gaps."""
        df = pd.DataFrame({
            'source_structure_id': [1, 1, 2, 2],
            'target_structure_id': [10, 20, 10, 20],
            'hemisphere': ['ipsi'] * 4,
            'metric': ['connection_strength'] * 4,
            'value': [1.0, 2.0, 3.0, 4.0],
        })
        # groups {1, 2} -> 'A', {10, 20} -> 'B': the whole table collapses to a single A->A cell
        group = {1: 'A', 2: 'A', 10: 'A', 20: 'A'}
        df['source_group'] = df['source_structure_id'].map(group)
        df['target_group'] = df['target_structure_id'].map(group)
        total = df.groupby(['source_group', 'target_group'])['value'].sum()
        self.assertEqual(total.loc[('A', 'A')], df['value'].sum())
        np.testing.assert_allclose(total.loc[('A', 'A')], 10.0)
