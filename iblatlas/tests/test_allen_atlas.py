import numpy as np
import matplotlib.pyplot as plt
import unittest
from iblatlas.atlas import AllenAtlas, BrainRegions
from iblatlas.flatmaps import FlatMap
from iblatlas.plots import plot_swanson, annotate_swanson
from iblatlas.genomics import agea


class TestAtlasSlicesConversion(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.ba = AllenAtlas(25)
        self.ba.compute_surface()

    def test_compute_volume(self):
        self.ba.compute_regions_volume()
        # 'AAA' aid 23, index 594 1921, one sided volume should be around 1/4 of a sq. mm
        np.testing.assert_allclose(self.ba.regions.volume[594], 0.247453125)
        np.testing.assert_allclose(self.ba.regions.volume[31], 0)
        # test the cumulative hierarchy volume
        self.ba.compute_regions_volume(cumsum=True)
        np.testing.assert_allclose(self.ba.regions.volume[31], 16.666015625)

    def test_simple(self):
        ba = self.ba
        # extracts the top surface from the volume and make sure it's all populated
        ix, iy = np.meshgrid(np.arange(ba.bc.nx), np.arange(ba.bc.ny))
        iz = ba.bc.z2i(ba.top)
        inds = ba._lookup_inds(np.stack((ix, iy, iz), axis=-1))
        self.assertTrue(np.all(ba.label.flat[inds][~np.isnan(ba.top)] != 0))
        # one sample above, it's all zeros
        inds = ba._lookup_inds(np.stack((ix, iy, np.maximum(iz - 1, 0)), axis=-1))
        self.assertTrue(np.all(ba.label.flat[inds][~np.isnan(ba.top)] == 0))
        # plt.imshow(self._label2rgb(self.label.flat[inds]))  # show the surface
        # do the same for the bottom surface
        izb = ba.bc.z2i(ba.bottom)
        inds = ba._lookup_inds(np.stack((ix, iy, izb), axis=-1))
        self.assertTrue(np.all(ba.label.flat[inds][~np.isnan(ba.top)] != 0))
        # one sample below, it's all zeros
        inds = ba._lookup_inds(np.stack((ix, iy, np.maximum(izb + 1, 0)), axis=-1))
        self.assertTrue(np.all(ba.label.flat[inds][~np.isnan(ba.bottom)] == 0))

    def test_lookups(self):
        ba = self.ba
        # test the probabilistic indices lookup
        # radius_um = 200
        # mapping = 'Beryl'
        xyz = np.array([0, -.0058, -.0038])

        # from atlasview import atlasview  # mouais il va falloir changer ça
        # av = atlasview.view()  #

        aid = ba.get_labels(xyz, mapping='Beryl')
        aids, proportions = ba.get_labels(xyz, mapping='Beryl', radius_um=250)

        self.assertEqual(aid, 912)
        self.assertTrue(np.all(aids == np.array([997, 912, 976, 968])))
        expected = np.array([0.40709028, 0.35887036, 0.22757999, 0.00645937])
        np.testing.assert_allclose(proportions, expected, atol=1e-6)


class TestFlatMaps(unittest.TestCase):

    def test_flatmaps(self):
        fm = FlatMap(flatmap='dorsal_cortex')
        fm.plot_flatmap(depth=0)

        fm = FlatMap(flatmap='circles')
        fm.plot_flatmap()

        fm = FlatMap(flatmap='pyramid')
        fm.plot_flatmap(volume='image')


class TestSwanson(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.regions = BrainRegions()

    def test_swanson(self):
        br = self.regions

        # prepare array of acronyms
        acronyms = ['ACAd1', 'ACAv1', 'AId1', 'AIp1', 'AIv1', 'AUDd1', 'AUDp1', 'AUDpo1', 'AUDv1',
                    'SSp-m1', 'SSp-n1', 'SSp-tr1', 'SSp-ul1', 'SSp-un1', 'SSs1',
                    'VISC1', 'VISa1', 'VISal1', 'VISam1', 'VISl1', 'VISli1', 'VISp1', 'VISp2/3', 'VISpl1', 'VISpm1',
                    'SSp-n2/3', 'SSp-tr2/3', 'SSp-ul2/3', 'SSp-un2/3', 'SSs2/3',
                    'VISC2/3', 'VISa2/3', 'VISal2/3', 'VISam2/3', 'VISl2/3', 'VISli2/3', 'VISp2/3', 'VISpl1', 'VISpl2/3']

        # assign data to each acronym
        values = np.arange(len(acronyms))

        regions_rl = np.r_[br.acronym2id(acronyms), -br.acronym2id(acronyms)]
        values_rl = np.random.randn(regions_rl.size)

        fig, ax = plt.subplots()
        plot_swanson(hemisphere='both', br=br)
        fig, ax = plt.subplots()
        plot_swanson(br=br)
        fig, ax = plt.subplots()
        plot_swanson(acronyms, values, cmap='Blues', br=br)
        fig, ax = plt.subplots()
        plot_swanson(acronyms, values, cmap='Blues', hemisphere='both')
        fig, ax = plt.subplots()
        plot_swanson(regions_rl, values_rl, hemisphere='both', cmap='magma', br=br)
        annotate_swanson(ax, acronyms=acronyms)
        fig, ax = plt.subplots()
        plot_swanson(br=br, annotate=True)
        plt.close('all')


class TestGeneExpression(unittest.TestCase):

    def test_load(self):
        df_genes, gene_expression, atlas_agea = agea.load()
        self.assertEqual(df_genes.shape[0], gene_expression.shape[0])
        self.assertEqual(gene_expression.shape[1:], (58, 41, 67))
