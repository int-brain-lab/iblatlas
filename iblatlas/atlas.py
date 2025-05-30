"""
Classes for manipulating brain atlases, insertions, and coordinates.
"""
from pathlib import Path, PurePosixPath
from dataclasses import dataclass
import logging

import matplotlib.pyplot as plt
import numpy as np
import nrrd

from one.webclient import http_download_file
import one.params
import one.remote.aws as aws
from iblutil.numerical import ismember
from iblatlas.regions import BrainRegions, FranklinPaxinosRegions

ALLEN_CCF_LANDMARKS_MLAPDV_UM = {'bregma': np.array([5739, 5400, 332])}
"""dict: The ML AP DV voxel coordinates of brain landmarks in the Allen atlas."""

PAXINOS_CCF_LANDMARKS_MLAPDV_UM = {'bregma': np.array([5700, 4300 + 160, 330])}
"""dict: The ML AP DV voxel coordinates of brain landmarks in the Franklin & Paxinos atlas."""

S3_BUCKET_IBL = 'ibl-brain-wide-map-public'
"""str: The name of the public IBL S3 bucket containing atlas data."""

_logger = logging.getLogger(__name__)


def cart2sph(x, y, z):
    """
    Converts cartesian to spherical coordinates.

    Returns spherical coordinates (r, theta, phi).

    Parameters
    ----------
    x : numpy.array
        A 1D array of x-axis coordinates.
    y : numpy.array
        A 1D array of y-axis coordinates.
    z : numpy.array
        A 1D array of z-axis coordinates.

    Returns
    -------
    numpy.array
        The radial distance of each point.
    numpy.array
        The polar angle.
    numpy.array
        The azimuthal angle.

    See Also
    --------
    sph2cart
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    phi = np.arctan2(y, x) * 180 / np.pi
    theta = np.zeros_like(r)
    iok = r != 0
    theta[iok] = np.arccos(z[iok] / r[iok]) * 180 / np.pi
    if theta.size == 1:
        theta = float(theta)
    return r, theta, phi


def sph2cart(r, theta, phi):
    """
    Converts Spherical to Cartesian coordinates.

    Returns Cartesian coordinates (x, y, z).

    Parameters
    ----------
    r : numpy.array
        A 1D array of radial distances.
    theta : numpy.array
        A 1D array of polar angles.
    phi : numpy.array
        A 1D array of azimuthal angles.

    Returns
    -------
    x : numpy.array
        A 1D array of x-axis coordinates.
    y : numpy.array
        A 1D array of y-axis coordinates.
    z : numpy.array
        A 1D array of z-axis coordinates.

    See Also
    --------
    cart2sph
    """
    x = r * np.cos(phi / 180 * np.pi) * np.sin(theta / 180 * np.pi)
    y = r * np.sin(phi / 180 * np.pi) * np.sin(theta / 180 * np.pi)
    z = r * np.cos(theta / 180 * np.pi)
    return x, y, z


def rodrigues_rotation(v: np.ndarray, k: np.ndarray, theta: float) -> np.ndarray:
    """
    Rotate vector v around axis k by angle theta using Rodrigues' rotation formula.

    Parameters
    ----------
    v : np.ndarray
        The vector to be rotated, shape (3,)
    k : np.ndarray
        The axis of rotation (should be a unit vector), shape (3,)
    theta : float
        The angle of rotation in radians

    Returns
    -------
    np.ndarray
        The rotated vector, shape (3,)

    Notes
    -----
    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    """
    k = k / np.linalg.norm(k)  # Ensure k is a unit vector
    return (v * np.cos(theta) +
            np.cross(k, v) * np.sin(theta) +
            k * np.dot(k, v) * (1 - np.cos(theta)))


def tilt_spherical(theta: float, phi: float, tilt_angle: float = -5) -> tuple[float, float]:
    """
    Rotate the coordinates (theta, phi) around the x-axis by tilt_angle degrees.

    Parameters
    ----------
    theta : float
        Polar angle from z+ (up) in degrees
    phi : float
        Azimuthal angle from x+ (right), counter-clockwise, in degrees
    tilt_angle : float, optional
        The angle of rotation around the x-axis in degrees. Default is -5.
        A positive tilt_angle lifts the mouse chin upwards.

    Returns
    -------
    theta_ : float
        The new polar angle in degrees after tilt correction
    phi_ : float
        The new azimuthal angle in degrees after tilt correction,
        normalized to the range [0, 360)

    Notes
    -----
    This function uses the Rodrigues rotation formula to apply the tilt.
    """
    v = np.array(sph2cart(1, theta, phi))  # unit vector
    k = np.array([1, 0, 0])  # rotation axis
    vrot = rodrigues_rotation(v, k, np.deg2rad(tilt_angle))
    _, theta_, phi_ = cart2sph(*vrot)
    return theta_, (phi_ % 360 + 360) % 360


class BrainCoordinates:
    """
    Class for mapping and indexing a 3D array to real-world coordinates.

    * x = ml, right positive
    * y = ap, anterior positive
    * z = dv, dorsal positive

    The layout of the Atlas dimension is done according to the most used sections so they lay
    contiguous on disk assuming C-ordering: V[iap, iml, idv]

    Parameters
    ----------
    nxyz : array_like
        Number of elements along each Cartesian axis (nx, ny, nz) = (nml, nap, ndv).
    xyz0 : array_like
        Coordinates of the element volume[0, 0, 0] in the coordinate space.
    dxyz : array_like, float
        Spatial interval of the volume along the 3 dimensions.

    Attributes
    ----------
    xyz0 : numpy.array
        The Cartesian coordinates of the element volume[0, 0, 0], i.e. the origin.
    x0 : int
        The x-axis origin coordinate of the element volume.
    y0 : int
        The y-axis origin coordinate of the element volume.
    z0 : int
        The z-axis origin coordinate of the element volume.
    """

    def __init__(self, nxyz, xyz0=(0, 0, 0), dxyz=(1, 1, 1)):
        if np.isscalar(dxyz):
            dxyz = [dxyz] * 3
        self.x0, self.y0, self.z0 = list(xyz0)
        self.dx, self.dy, self.dz = list(dxyz)
        self.nx, self.ny, self.nz = list(nxyz)

    @property
    def dxyz(self):
        """numpy.array: Spatial interval of the volume along the 3 dimensions."""
        return np.array([self.dx, self.dy, self.dz])

    @property
    def nxyz(self):
        """numpy.array: Coordinates of the element volume[0, 0, 0] in the coordinate space."""
        return np.array([self.nx, self.ny, self.nz])

    """Methods distance to indices"""
    @staticmethod
    def _round(i, round=True):
        """
        Round an input value to the nearest integer, replacing NaN values with 0.

        Parameters
        ----------
        i : int, float, numpy.nan, numpy.array
            A value or array of values to round.
        round : bool
            If false this function is identity.

        Returns
        -------
        int, float, numpy.nan, numpy.array
            If round is true, returns the nearest integer, replacing NaN values with 0, otherwise
            returns the input unaffected.
        """
        nanval = 0
        if round:
            ii = np.array(np.round(i)).astype(int)
            ii[np.isnan(i)] = nanval
            return ii
        else:
            return i

    def x2i(self, x, round=True, mode='raise'):
        """
        Find the nearest volume image index to a given x-axis coordinate.

        Parameters
        ----------
        x : float, numpy.array
            One or more x-axis coordinates, relative to the origin, x0.
        round : bool
            If true, round to the nearest index, replacing NaN values with 0.
        mode : {'raise', 'clip', 'wrap'}, default='raise'
            How to behave if the coordinate lies outside of the volume: raise (default) will raise
            a ValueError; 'clip' will replace the index with the closest index inside the volume;
            'wrap' will return the index as is.

        Returns
        -------
        numpy.array
            The nearest indices of the image volume along the first dimension.

        Raises
        ------
        ValueError
            At least one x value lies outside of the atlas volume. Change 'mode' input to 'wrap' to
            keep these values unchanged, or 'clip' to return the nearest valid indices.
        """
        i = np.asarray(self._round((x - self.x0) / self.dx, round=round))
        if np.any(i < 0) or np.any(i >= self.nx):
            if mode == 'clip':
                i[i < 0] = 0
                i[i >= self.nx] = self.nx - 1
            elif mode == 'raise':
                raise ValueError("At least one x value lies outside of the atlas volume.")
            elif mode == 'wrap':  # This is only here for legacy reasons
                pass
        return i

    def y2i(self, y, round=True, mode='raise'):
        """
        Find the nearest volume image index to a given y-axis coordinate.

        Parameters
        ----------
        y : float, numpy.array
            One or more y-axis coordinates, relative to the origin, y0.
        round : bool
            If true, round to the nearest index, replacing NaN values with 0.
        mode : {'raise', 'clip', 'wrap'}
            How to behave if the coordinate lies outside of the volume: raise (default) will raise
            a ValueError; 'clip' will replace the index with the closest index inside the volume;
            'wrap' will return the index as is.

        Returns
        -------
        numpy.array
            The nearest indices of the image volume along the second dimension.

        Raises
        ------
        ValueError
            At least one y value lies outside of the atlas volume. Change 'mode' input to 'wrap' to
            keep these values unchanged, or 'clip' to return the nearest valid indices.
        """
        i = np.asarray(self._round((y - self.y0) / self.dy, round=round))
        if np.any(i < 0) or np.any(i >= self.ny):
            if mode == 'clip':
                i[i < 0] = 0
                i[i >= self.ny] = self.ny - 1
            elif mode == 'raise':
                raise ValueError("At least one y value lies outside of the atlas volume.")
            elif mode == 'wrap':  # This is only here for legacy reasons
                pass
        return i

    def z2i(self, z, round=True, mode='raise'):
        """
        Find the nearest volume image index to a given z-axis coordinate.

        Parameters
        ----------
        z : float, numpy.array
            One or more z-axis coordinates, relative to the origin, z0.
        round : bool
            If true, round to the nearest index, replacing NaN values with 0.
        mode : {'raise', 'clip', 'wrap'}
            How to behave if the coordinate lies outside of the volume: raise (default) will raise
            a ValueError; 'clip' will replace the index with the closest index inside the volume;
            'wrap' will return the index as is.

        Returns
        -------
        numpy.array
            The nearest indices of the image volume along the third dimension.

        Raises
        ------
        ValueError
            At least one z value lies outside of the atlas volume. Change 'mode' input to 'wrap' to
            keep these values unchanged, or 'clip' to return the nearest valid indices.
        """
        i = np.asarray(self._round((z - self.z0) / self.dz, round=round))
        if np.any(i < 0) or np.any(i >= self.nz):
            if mode == 'clip':
                i[i < 0] = 0
                i[i >= self.nz] = self.nz - 1
            elif mode == 'raise':
                raise ValueError("At least one z value lies outside of the atlas volume.")
            elif mode == 'wrap':  # This is only here for legacy reasons
                pass
        return i

    def xyz2i(self, xyz, round=True, mode='raise'):
        """
        Find the nearest volume image indices to the given Cartesian coordinates.

        Parameters
        ----------
        xyz : array_like
            One or more Cartesian coordinates, relative to the origin, xyz0.
        round : bool
            If true, round to the nearest index, replacing NaN values with 0.
        mode : {'raise', 'clip', 'wrap'}
            How to behave if any coordinate lies outside of the volume: raise (default) will raise
            a ValueError; 'clip' will replace the index with the closest index inside the volume;
            'wrap' will return the index as is.

        Returns
        -------
        numpy.array
            The nearest indices of the image volume.

        Raises
        ------
        ValueError
            At least one coordinate lies outside of the atlas volume. Change 'mode' input to 'wrap'
            to keep these values unchanged, or 'clip' to return the nearest valid indices.
        """
        xyz = np.array(xyz)
        dt = int if round else float
        out = np.zeros_like(xyz, dtype=dt)
        out[..., 0] = self.x2i(xyz[..., 0], round=round, mode=mode)
        out[..., 1] = self.y2i(xyz[..., 1], round=round, mode=mode)
        out[..., 2] = self.z2i(xyz[..., 2], round=round, mode=mode)
        return out

    """Methods indices to distance"""
    def i2x(self, ind):
        """
        Return the x-axis coordinate of a given index.

        Parameters
        ----------
        ind : int, numpy.array
            One or more indices along the first dimension of the image volume.

        Returns
        -------
        float, numpy.array
            The corresponding x-axis coordinate(s), relative to the origin, x0.
        """
        return ind * self.dx + self.x0

    def i2y(self, ind):
        """
        Return the y-axis coordinate of a given index.

        Parameters
        ----------
        ind : int, numpy.array
            One or more indices along the second dimension of the image volume.

        Returns
        -------
        float, numpy.array
            The corresponding y-axis coordinate(s), relative to the origin, y0.
        """
        return ind * self.dy + self.y0

    def i2z(self, ind):
        """
        Return the z-axis coordinate of a given index.

        Parameters
        ----------
        ind : int, numpy.array
            One or more indices along the third dimension of the image volume.

        Returns
        -------
        float, numpy.array
            The corresponding z-axis coordinate(s), relative to the origin, z0.
        """
        return ind * self.dz + self.z0

    def i2xyz(self, iii):
        """
        Return the Cartesian coordinates of a given index.

        Parameters
        ----------
        iii : array_like
            One or more image volume indices.

        Returns
        -------
        numpy.array
            The corresponding xyz coordinates, relative to the origin, xyz0.
        """

        iii = np.array(iii, dtype=float)
        out = np.zeros_like(iii)
        out[..., 0] = self.i2x(iii[..., 0])
        out[..., 1] = self.i2y(iii[..., 1])
        out[..., 2] = self.i2z(iii[..., 2])
        return out

    """Methods bounds"""
    @property
    def xlim(self):
        # FIXME Document
        return self.i2x(np.array([0, self.nx - 1]))

    @property
    def ylim(self):
        # FIXME Document
        return self.i2y(np.array([0, self.ny - 1]))

    @property
    def zlim(self):
        # FIXME Document
        return self.i2z(np.array([0, self.nz - 1]))

    def lim(self, axis):
        # FIXME Document
        if axis == 0:
            return self.xlim
        elif axis == 1:
            return self.ylim
        elif axis == 2:
            return self.zlim

    """returns scales"""
    @property
    def xscale(self):
        # FIXME Document
        return self.i2x(np.arange(self.nx))

    @property
    def yscale(self):
        # FIXME Document
        return self.i2y(np.arange(self.ny))

    @property
    def zscale(self):
        # FIXME Document
        return self.i2z(np.arange(self.nz))

    """returns the 3d mgrid used for 3d visualization"""
    @property
    def mgrid(self):
        # FIXME Document
        return np.meshgrid(self.xscale, self.yscale, self.zscale)


class BrainAtlas:
    """
    Objects that holds image, labels and coordinate transforms for a brain Atlas.
    Currently this is designed for the AllenCCF at several resolutions,
    yet this class can be used for other atlases arises.
    """

    """numpy.array: A 3D image volume."""
    image = None
    """numpy.array: A 3D annotation label volume."""
    label = None
    """numpy.array: One or several optional data volumes, the 3 last dimensions should match
    the image and the label volumes dimensions"""
    volumes = None

    def __init__(self, image, label, dxyz, regions, iorigin=[0, 0, 0],
                 dims2xyz=[0, 1, 2], xyz2dims=[0, 1, 2]):
        """
        self.image: image volume (ap, ml, dv)
        self.label: label volume (ap, ml, dv)
        self.bc: atlas.BrainCoordinate object
        self.regions: atlas.BrainRegions object
        self.top: 2d np array (ap, ml) containing the z-coordinate (m) of the surface of the brain
        self.dims2xyz and self.zyz2dims: map image axis order to xyz coordinates order
        """

        self.image = image
        self.label = label
        self.regions = regions
        self.dims2xyz = np.array(dims2xyz)
        self.xyz2dims = np.array(xyz2dims)
        assert np.all(self.dims2xyz[self.xyz2dims] == np.array([0, 1, 2]))
        assert np.all(self.xyz2dims[self.dims2xyz] == np.array([0, 1, 2]))
        # create the coordinate transform object that maps volume indices to real world coordinates
        nxyz = np.array(self.image.shape)[self.dims2xyz]
        bc = BrainCoordinates(nxyz=nxyz, xyz0=(0, 0, 0), dxyz=dxyz)
        self.bc = BrainCoordinates(nxyz=nxyz, xyz0=-bc.i2xyz(iorigin), dxyz=dxyz)

        self.surface = None
        self.boundary = None

    @staticmethod
    def _get_cache_dir():
        """
        ./histology/ATLAS/Needles/Allen
        Where . is the main ONE cache directory
        :return: pathlib.Path
        """
        par = one.params.get(silent=True)
        path_atlas = Path(par.CACHE_DIR).joinpath('histology', 'ATLAS', 'Needles', 'Allen')
        return path_atlas

    def mask(self):
        """
        Returns a Boolean mask volume of the brain atlas, where 1 is inside the convex brain and 0 is outside
        This returns an ovoid volume shaped like the brain and this will contain void values in the ventricules.
        :return: np.array Bool (nap, nml, ndv)
        """
        self.compute_surface()
        mask = np.logical_and(np.cumsum(self.label != 0, axis=-1) > 0,
                              np.flip(np.cumsum(np.flip(self.label, axis=-1) != 0, axis=-1), axis=-1) > 0)
        mask[np.isnan(self.top)] = 0
        return mask

    def compute_surface(self):
        """
        Get the volume top, bottom, left and right surfaces, and from these the outer surface of
        the image volume. This is needed to compute probe insertions intersections.

        NOTE: In places where the top or bottom surface touch the top or bottom of the atlas volume, the surface
        will be set to np.nan. If you encounter issues working with these surfaces check if this might be the cause.
        """
        if self.surface is None:  # only compute if it hasn't already been computed
            axz = self.xyz2dims[2]  # this is the dv axis
            _surface = (self.label == 0).astype(np.int8) * 2
            l0 = np.diff(_surface, axis=axz, append=2)
            _top = np.argmax(l0 == -2, axis=axz).astype(float)
            _top[_top == 0] = np.nan
            _bottom = self.bc.nz - np.argmax(np.flip(l0, axis=axz) == 2, axis=axz).astype(float)
            _bottom[_bottom == self.bc.nz] = np.nan
            self.top = self.bc.i2z(_top + 1)
            self.bottom = self.bc.i2z(_bottom - 1)
            self.surface = np.diff(_surface, axis=self.xyz2dims[0], append=2) + l0
            idx_srf = np.where(self.surface != 0)
            self.surface[idx_srf] = 1
            self.srf_xyz = self.bc.i2xyz(np.c_[idx_srf[self.xyz2dims[0]], idx_srf[self.xyz2dims[1]],
                                               idx_srf[self.xyz2dims[2]]].astype(float))

    def _lookup_inds(self, ixyz, mode='raise'):
        """
        Performs a 3D lookup from volume indices ixyz to the image volume
        :param ixyz: [n, 3] array of indices in the mlapdv order
        :return: n array of flat indices
        """
        idims = np.split(ixyz[..., self.xyz2dims], [1, 2], axis=-1)
        inds = np.ravel_multi_index(idims, self.bc.nxyz[self.xyz2dims], mode=mode)
        return inds.squeeze()

    def _lookup(self, xyz, mode='raise'):
        """
        Performs a 3D lookup from real world coordinates to the flat indices in the volume,
        defined in the BrainCoordinates object.

        Parameters
        ----------
        xyz : numpy.array
            An (n, 3) array of Cartesian coordinates.
        mode : {'raise', 'clip', 'wrap'}
            How to behave if any coordinate lies outside of the volume: raise (default) will raise
            a ValueError; 'clip' will replace the index with the closest index inside the volume;
            'wrap' will return the index as is.

        Returns
        -------
        numpy.array
            A 1D array of flat indices.
        """
        return self._lookup_inds(self.bc.xyz2i(xyz, mode=mode), mode=mode)

    def get_labels(self, xyz, mapping=None, radius_um=None, mode='raise'):
        """
        Performs a 3D lookup from real world coordinates to the volume labels
        and return the regions ids according to the mapping
        :param xyz: [n, 3] array of coordinates
        :param mapping: brain region mapping (defaults to original Allen mapping)
        :param radius_um: if not null, returns a regions ids array and an array of proportion
         of regions in a sphere of size radius around the coordinates.
        :param mode: {‘raise’, 'clip'} determines what to do when determined index lies outside the atlas volume
        'raise' will raise a ValueError (default)
        'clip' will replace the index with the closest index inside the volume
        :return: n array of region ids
        """
        mapping = mapping or self.regions.default_mapping

        if radius_um:
            nrx = int(np.ceil(radius_um / abs(self.bc.dx) / 1e6))
            nry = int(np.ceil(radius_um / abs(self.bc.dy) / 1e6))
            nrz = int(np.ceil(radius_um / abs(self.bc.dz) / 1e6))
            nr = [nrx, nry, nrz]
            iii = self.bc.xyz2i(xyz, mode=mode)
            # computing the cube radius and indices is more complicated as volume indices are not
            # necessarily in ml, ap, dv order so the indices order is dynamic
            rcube = np.meshgrid(*tuple((np.arange(
                -nr[i], nr[i] + 1) * self.bc.dxyz[i]) ** 2 for i in self.xyz2dims))
            rcube = np.sqrt(rcube[0] + rcube[1], rcube[2]) * 1e6
            icube = tuple(slice(-nr[i] + iii[i], nr[i] + iii[i] + 1) for i in self.xyz2dims)
            cube = self.regions.mappings[mapping][self.label[icube]]
            ilabs, counts = np.unique(cube[rcube <= radius_um], return_counts=True)
            return self.regions.id[ilabs], counts / np.sum(counts)
        else:
            regions_indices = self._get_mapping(mapping=mapping)[self.label.flat[self._lookup(xyz, mode=mode)]]
            return self.regions.id[regions_indices]

    def _get_mapping(self, mapping=None):
        """
        Safe way to get mappings if nothing defined in regions.
        A mapping transforms from the full allen brain Atlas ids to the remapped ids
        new_ids = ids[mapping]
        """
        mapping = mapping or self.regions.default_mapping
        if hasattr(self.regions, 'mappings'):
            return self.regions.mappings[mapping]
        else:
            return np.arange(self.regions.id.size)

    def _label2rgb(self, imlabel):
        """
        Converts a slice from the label volume to its RGB equivalent for display
        :param imlabel: 2D np-array containing label ids (slice of the label volume)
        :return: 3D np-array of the slice uint8 rgb values
        """
        if getattr(self.regions, 'rgb', None) is None:
            return self.regions.id[imlabel]
        else:  # if the regions exist and have the rgb attribute, do the rgb lookup
            return self.regions.rgb[imlabel]

    def tilted_slice(self, xyz, axis, volume='image'):
        """
        From line coordinates, extracts the tilted plane containing the line from the 3D volume
        :param xyz: np.array: points defining a probe trajectory in 3D space (xyz triplets)
        if more than 2 points are provided will take the best fit
        :param axis:
            0: along ml = sagittal-slice
            1: along ap = coronal-slice
            2: along dv = horizontal-slice
        :param volume: 'image' or 'annotation'
        :return: np.array, abscissa extent (width), ordinate extent (height),
        squeezed axis extent (depth)
        """
        if axis == 0:   # sagittal slice (squeeze/take along ml-axis)
            wdim, hdim, ddim = (1, 2, 0)
        elif axis == 1:  # coronal slice (squeeze/take along ap-axis)
            wdim, hdim, ddim = (0, 2, 1)
        elif axis == 2:  # horizontal slice (squeeze/take along dv-axis)
            wdim, hdim, ddim = (0, 1, 2)
        # get the best fit and find exit points of the volume along squeezed axis
        trj = Trajectory.fit(xyz)
        sub_volume = trj._eval(self.bc.lim(axis=hdim), axis=hdim)
        sub_volume[:, wdim] = self.bc.lim(axis=wdim)
        sub_volume_i = self.bc.xyz2i(sub_volume)
        tile_shape = np.array([np.diff(sub_volume_i[:, hdim])[0] + 1, self.bc.nxyz[wdim]])
        # get indices along each dimension
        indx = np.arange(tile_shape[1])
        indy = np.arange(tile_shape[0])
        inds = np.linspace(*sub_volume_i[:, ddim], tile_shape[0])
        # compute the slice indices and output the slice
        _, INDS = np.meshgrid(indx, np.int64(np.around(inds)))
        INDX, INDY = np.meshgrid(indx, indy)
        indsl = [[INDX, INDY, INDS][i] for i in np.argsort([wdim, hdim, ddim])[self.xyz2dims]]
        if isinstance(volume, np.ndarray):
            tslice = volume[indsl[0], indsl[1], indsl[2]]
        elif volume.lower() == 'annotation':
            tslice = self._label2rgb(self.label[indsl[0], indsl[1], indsl[2]])
        elif volume.lower() == 'image':
            tslice = self.image[indsl[0], indsl[1], indsl[2]]
        elif volume.lower() == 'surface':
            tslice = self.surface[indsl[0], indsl[1], indsl[2]]

        #  get extents with correct convention NB: matplotlib flips the y-axis on imshow !
        width = np.sort(sub_volume[:, wdim])[np.argsort(self.bc.lim(axis=wdim))]
        height = np.flipud(np.sort(sub_volume[:, hdim])[np.argsort(self.bc.lim(axis=hdim))])
        depth = np.flipud(np.sort(sub_volume[:, ddim])[np.argsort(self.bc.lim(axis=ddim))])
        return tslice, width, height, depth

    def plot_tilted_slice(self, xyz, axis, volume='image', cmap=None, ax=None, return_sec=False, **kwargs):
        """
        From line coordinates, extracts the tilted plane containing the line from the 3D volume
        :param xyz: np.array: points defining a probe trajectory in 3D space (xyz triplets)
        if more than 2 points are provided will take the best fit
        :param axis:
            0: along ml = sagittal-slice
            1: along ap = coronal-slice
            2: along dv = horizontal-slice
        :param volume: 'image' or 'annotation'
        :return: matplotlib axis
        """
        if axis == 0:
            axis_labels = np.array(['ap (um)', 'dv (um)', 'ml (um)'])
        elif axis == 1:
            axis_labels = np.array(['ml (um)', 'dv (um)', 'ap (um)'])
        elif axis == 2:
            axis_labels = np.array(['ml (um)', 'ap (um)', 'dv (um)'])

        tslice, width, height, depth = self.tilted_slice(xyz, axis, volume=volume)
        width = width * 1e6
        height = height * 1e6
        depth = depth * 1e6
        if not ax:
            plt.figure()
            ax = plt.gca()
            ax.axis('equal')
        if not cmap:
            cmap = plt.get_cmap('bone')
        # get the transfer function from y-axis to squeezed axis for second axe
        ab = np.linalg.solve(np.c_[height, height * 0 + 1], depth)
        height * ab[0] + ab[1]
        ax.imshow(tslice, extent=np.r_[width, height], cmap=cmap, **kwargs)
        sec_ax = ax.secondary_yaxis('right', functions=(
                                    lambda x: x * ab[0] + ab[1],
                                    lambda y: (y - ab[1]) / ab[0]))
        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        sec_ax.set_ylabel(axis_labels[2])
        if return_sec:
            return ax, sec_ax
        else:
            return ax

    @staticmethod
    def _plot_slice(im, extent, ax=None, cmap=None, volume=None, **kwargs):
        """
        Plot an atlas slice.

        Parameters
        ----------
        im : numpy.array
            A 2D image slice to plot.
        extent : array_like
            The bounding box in data coordinates that the image will fill specified as (left,
            right, bottom, top) in data coordinates.
        ax : matplotlib.pyplot.Axes
            An optional Axes object to plot to.
        cmap : str, matplotlib.colors.Colormap
            The Colormap instance or registered colormap name used to map scalar data to colors.
            Defaults to 'bone'.
        volume : str | np.array
            If 'boundary', assumes image is an outline of boundaries between all regions.
            FIXME How does this affect the plot?
        **kwargs
            See matplotlib.pyplot.imshow.

        Returns
        -------
        matplotlib.pyplot.Axes
            The image axes.
        """
        if not ax:
            ax = plt.gca()
            ax.axis('equal')
        if not cmap:
            cmap = plt.get_cmap('bone')

        if isinstance(volume, str):
            if volume == 'boundary':
                imb = np.zeros((*im.shape[:2], 4), dtype=np.uint8)
                imb[im == 1] = np.array([0, 0, 0, 255])
                im = imb

        ax.imshow(im, extent=extent, cmap=cmap, **kwargs)
        return ax

    def extent(self, axis):
        """
        :param axis: direction along which the volume is stacked:
         (2 = z for horizontal slice)
         (1 = y for coronal slice)
         (0 = x for sagittal slice)
        :return:
        """

        if axis == 0:
            extent = np.r_[self.bc.ylim, np.flip(self.bc.zlim)] * 1e6
        elif axis == 1:
            extent = np.r_[self.bc.xlim, np.flip(self.bc.zlim)] * 1e6
        elif axis == 2:
            extent = np.r_[self.bc.xlim, np.flip(self.bc.ylim)] * 1e6
        return extent

    def slice(self, coordinate, axis, volume='image', mode='raise', region_values=None,
              mapping=None, bc=None):
        """
        Get slice through atlas

        :param coordinate: coordinate to slice in metres, float
        :param axis: xyz convention:  0 for ml, 1 for ap, 2 for dv
            - 0: sagittal slice (along ml axis)
            - 1: coronal slice (along ap axis)
            - 2: horizontal slice (along dv axis)
        :param volume:
            - 'image' - allen image volume
            - 'rindex' - brain region index value
            - 'annotation' - allen annotation volume
            - 'surface' - outer surface of mesh
            - 'boundary' - outline of boundaries between all regions
            - 'volume' - custom volume, must pass in volume of shape ba.image.shape as regions_value argument
            - 'value' - custom value per allen region, must pass in array of shape ba.regions.id as regions_value argument
        :param mode: error mode for out of bounds coordinates
            -   'raise' raise an error
            -   'clip' gets the first or last index
        :param region_values: custom values to plot
            - if volume='volume', region_values must have shape ba.image.shape
            - if volume='value', region_values must have shape ba.regions.id
        :param mapping: mapping to use. Options can be found using ba.regions.mappings.keys()
        :return: 2d array or 3d RGB numpy int8 array of dimensions:
            - 0: nap x ndv (sagittal slice)
            - 1: nml x ndv (coronal slice)
            - 2: nap x nml (horizontal slice)
        """
        if axis == 0:
            index = self.bc.x2i(np.array(coordinate), mode=mode)
        elif axis == 1:
            index = self.bc.y2i(np.array(coordinate), mode=mode)
        elif axis == 2:
            index = self.bc.z2i(np.array(coordinate), mode=mode)

        def _take(vol, ind, axis):
            """
            This is a 2 steps process to get the slice along the correct axis
            1) slice the volume according to the mapped axis corresponding to the sclice
            we do this because np.take is 50 thousand times slower than straight slicing !
            2) reshape the output array according to the slice specifications
            """
            volume_axis = self.xyz2dims[axis]
            if mode == 'clip':
                ind = np.minimum(np.maximum(ind, 0), vol.shape[volume_axis] - 1)
            if volume_axis == 0:
                slic = vol[ind, :, :]
            elif volume_axis == 1:
                slic = vol[:, ind, :]
            elif volume_axis == 2:
                slic = vol[:, :, ind]
            output_sizes = [[1, 2], [0, 2], [1, 0]]  # we expect those sizes where index is the axis
            if np.diff(self.xyz2dims[output_sizes[axis]])[0] < 0:
                slic = slic.transpose().copy()
            return slic

        def _take_remap(vol, ind, axis, mapping):
            # For the labels, remap the regions indices according to the mapping
            return self._get_mapping(mapping=mapping)[_take(vol, ind, axis)]

        if isinstance(volume, np.ndarray):
            return _take(volume, index, axis=axis)
        elif volume == 'rindex':
            return _take_remap(self.label, index, axis=axis, mapping=mapping)
        elif volume in 'annotation':
            iregion = _take_remap(self.label, index, axis=axis, mapping=mapping)
            return self._label2rgb(iregion)
        elif volume == 'image':
            return _take(self.image, index, axis=axis)
        elif volume == 'value':
            return region_values[_take_remap(self.label, index, axis, mapping)]
        elif volume == 'image':
            return _take(self.image, index, axis=axis)
        elif volume in ['surface', 'edges']:
            self.compute_surface()
            return _take(self.surface, index, axis=axis)
        elif volume == 'boundary':
            iregion = _take_remap(self.label, index, axis, mapping)
            return self.compute_boundaries(iregion)

        elif volume == 'volume':
            if bc is not None:
                index = bc.xyz2i(np.array([coordinate] * 3))[axis]
            return _take(region_values, index, axis=axis)

    @staticmethod
    def compute_boundaries(values):
        """
        Compute the boundaries between regions on slice
        :param values:
        :return:
        """
        boundary = np.abs(np.diff(values, axis=0, prepend=0))
        boundary = boundary + np.abs(np.diff(values, axis=1, prepend=0))
        boundary = boundary + np.abs(np.diff(values, axis=1, append=0))
        boundary = boundary + np.abs(np.diff(values, axis=0, append=0))

        boundary[boundary != 0] = 1

        return boundary

    def plot_slices(self, xyz, *args, **kwargs):
        """
        From a single coordinate, plots the 3 slices that intersect at this point in a single
        matplotlib figure
        :param xyz: mlapdv coordinate in m
        :param args: arguments to be forwarded to plot slices
        :param kwargs: keyword arguments to be forwarded to plot slices
        :return: 2 by 2 array of axes
        """
        fig, axs = plt.subplots(2, 2)
        self.plot_cslice(xyz[1], *args, ax=axs[0, 0], **kwargs)
        self.plot_sslice(xyz[0], *args, ax=axs[0, 1], **kwargs)
        self.plot_hslice(xyz[2], *args, ax=axs[1, 0], **kwargs)
        xyz_um = xyz * 1e6
        axs[0, 0].plot(xyz_um[0], xyz_um[2], 'g*')
        axs[0, 1].plot(xyz_um[1], xyz_um[2], 'g*')
        axs[1, 0].plot(xyz_um[0], xyz_um[1], 'g*')
        return axs

    def plot_cslice(self, ap_coordinate, volume='image', mapping=None, region_values=None, **kwargs):
        """
        Plot coronal slice through atlas at given ap_coordinate

        :param: ap_coordinate (m)
        :param volume:
            - 'image' - allen image volume
            - 'annotation' - allen annotation volume
            - 'surface' - outer surface of mesh
            - 'boundary' - outline of boundaries between all regions
            - 'volume' - custom volume, must pass in volume of shape ba.image.shape as regions_value argument
            - 'value' - custom value per allen region, must pass in array of shape ba.regions.id as regions_value argument
        :param mapping: mapping to use. Options can be found using ba.regions.mappings.keys()
        :param region_values: custom values to plot
            - if volume='volume', region_values must have shape ba.image.shape
            - if volume='value', region_values must have shape ba.regions.id
        :param mapping: mapping to use. Options can be found using ba.regions.mappings.keys()
        :param **kwargs: matplotlib.pyplot.imshow kwarg arguments
        :return: matplotlib ax object
        """
        cslice = self.slice(ap_coordinate, axis=1, volume=volume, mapping=mapping, region_values=region_values)
        return self._plot_slice(np.moveaxis(cslice, 0, 1), extent=self.extent(axis=1), volume=volume, **kwargs)

    def plot_hslice(self, dv_coordinate, volume='image', mapping=None, region_values=None, **kwargs):
        """
        Plot horizontal slice through atlas at given dv_coordinate

        :param: dv_coordinate (m)
        :param volume:
            - 'image' - allen image volume
            - 'annotation' - allen annotation volume
            - 'surface' - outer surface of mesh
            - 'boundary' - outline of boundaries between all regions
            - 'volume' - custom volume, must pass in volume of shape ba.image.shape as regions_value argument
            - 'value' - custom value per allen region, must pass in array of shape ba.regions.id as regions_value argument
        :param mapping: mapping to use. Options can be found using ba.regions.mappings.keys()
        :param region_values: custom values to plot
            - if volume='volume', region_values must have shape ba.image.shape
            - if volume='value', region_values must have shape ba.regions.id
        :param mapping: mapping to use. Options can be found using ba.regions.mappings.keys()
        :param **kwargs: matplotlib.pyplot.imshow kwarg arguments
        :return: matplotlib ax object
        """

        hslice = self.slice(dv_coordinate, axis=2, volume=volume, mapping=mapping, region_values=region_values)
        return self._plot_slice(hslice, extent=self.extent(axis=2), volume=volume, **kwargs)

    def plot_sslice(self, ml_coordinate, volume='image', mapping=None, region_values=None, **kwargs):
        """
        Plot sagittal slice through atlas at given ml_coordinate

        :param: ml_coordinate (m)
        :param volume:
            - 'image' - allen image volume
            - 'annotation' - allen annotation volume
            - 'surface' - outer surface of mesh
            - 'boundary' - outline of boundaries between all regions
            - 'volume' - custom volume, must pass in volume of shape ba.image.shape as regions_value argument
            - 'value' - custom value per allen region, must pass in array of shape ba.regions.id as regions_value argument
        :param mapping: mapping to use. Options can be found using ba.regions.mappings.keys()
        :param region_values: custom values to plot
            - if volume='volume', region_values must have shape ba.image.shape
            - if volume='value', region_values must have shape ba.regions.id
        :param mapping: mapping to use. Options can be found using ba.regions.mappings.keys()
        :param **kwargs: matplotlib.pyplot.imshow kwarg arguments
        :return: matplotlib ax object
        """

        sslice = self.slice(ml_coordinate, axis=0, volume=volume, mapping=mapping, region_values=region_values)
        return self._plot_slice(np.swapaxes(sslice, 0, 1), extent=self.extent(axis=0), volume=volume, **kwargs)

    def plot_top(self, volume='annotation', mapping=None, region_values=None, ax=None, **kwargs):
        """
        Plot top view of atlas
        :param volume:
            - 'image' - allen image volume
            - 'annotation' - allen annotation volume
            - 'boundary' - outline of boundaries between all regions
            - 'volume' - custom volume, must pass in volume of shape ba.image.shape as regions_value argument
            - 'value' - custom value per allen region, must pass in array of shape ba.regions.id as regions_value argument

        :param mapping: mapping to use. Options can be found using ba.regions.mappings.keys()
        :param region_values:
        :param ax:
        :param kwargs:
        :return:
        """
        self.compute_surface()
        ix, iy = np.meshgrid(np.arange(self.bc.nx), np.arange(self.bc.ny))
        iz = self.bc.z2i(self.top)
        inds = self._lookup_inds(np.stack((ix, iy, iz), axis=-1))
        regions = self._get_mapping(mapping=mapping)[self.label.flat[inds]]
        if volume == 'annotation':
            im = self._label2rgb(regions)
        elif volume == 'image':
            im = self.top
        elif volume == 'value':
            im = region_values[regions]
        elif volume == 'volume':
            im = np.zeros((iz.shape))
            for x in range(im.shape[0]):
                for y in range(im.shape[1]):
                    im[x, y] = region_values[x, y, iz[x, y]]
        elif volume == 'boundary':
            im = self.compute_boundaries(regions)

        return self._plot_slice(im, self.extent(axis=2), ax=ax, volume=volume, **kwargs)


@dataclass
class Trajectory:
    """
    3D Trajectory (usually for a linear probe), minimally defined by a vector and a point.

    Examples
    --------
    Instantiate from a best fit from an n by 3 array containing xyz coordinates:

    >>> trj = Trajectory.fit(xyz)
    """
    vector: np.ndarray
    point: np.ndarray

    @staticmethod
    def fit(xyz):
        """
        Fits a line to a 3D cloud of points.

        Parameters
        ----------
        xyz : numpy.array
            An n by 3 array containing a cloud of points to fit a line to.

        Returns
        -------
        Trajectory
            A new trajectory object.
        """
        xyz_mean = np.mean(xyz, axis=0)
        return Trajectory(vector=np.linalg.svd(xyz - xyz_mean)[2][0], point=xyz_mean)

    def eval_x(self, x):
        """
        given an array of x coordinates, returns the xyz array of coordinates along the insertion
        :param x: n by 1 or numpy array containing x-coordinates
        :return: n by 3 numpy array containing xyz-coordinates
        """
        return self._eval(x, axis=0)

    def eval_y(self, y):
        """
        given an array of y coordinates, returns the xyz array of coordinates along the insertion
        :param y: n by 1 or numpy array containing y-coordinates
        :return: n by 3 numpy array containing xyz-coordinates
        """
        return self._eval(y, axis=1)

    def eval_z(self, z):
        """
        given an array of z coordinates, returns the xyz array of coordinates along the insertion
        :param z: n by 1 or numpy array containing z-coordinates
        :return: n by 3 numpy array containing xyz-coordinates
        """
        return self._eval(z, axis=2)

    def project(self, point):
        """
        projects a point onto the trajectory line
        :param point: np.array(x, y, z) coordinates
        :return:
        """
        # https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        if point.ndim == 1:
            return self.project(point[np.newaxis])[0]
        return (self.point + np.dot(point[:, np.newaxis] - self.point, self.vector) /
                np.dot(self.vector, self.vector) * self.vector)

    def mindist(self, xyz, bounds=None):
        """
        Computes the minimum distance to the trajectory line for one or a set of points.
        If bounds are provided, computes the minimum distance to the segment instead of an
        infinite line.
        :param xyz: [..., 3]
        :param bounds: defaults to None.  np.array [2, 3]: segment boundaries, inf line if None
        :return: minimum distance [...]
        """
        proj = self.project(xyz)
        d = np.sqrt(np.sum((proj - xyz) ** 2, axis=-1))
        if bounds is not None:
            # project the boundaries and the points along the traj
            b = np.dot(bounds, self.vector)
            ob = np.argsort(b)
            p = np.dot(xyz[:, np.newaxis], self.vector).squeeze()
            # for points below and above boundaries, compute cartesian distance to the boundary
            imin = p < np.min(b)
            d[imin] = np.sqrt(np.sum((xyz[imin, :] - bounds[ob[0], :]) ** 2, axis=-1))
            imax = p > np.max(b)
            d[imax] = np.sqrt(np.sum((xyz[imax, :] - bounds[ob[1], :]) ** 2, axis=-1))
        return d

    def _eval(self, c, axis):
        # uses symmetric form of 3d line equation to get xyz coordinates given one coordinate
        if not isinstance(c, np.ndarray):
            c = np.array(c)
        while c.ndim < 2:
            c = c[..., np.newaxis]
        # there are cases where it's impossible to project if a line is // to the axis
        if self.vector[axis] == 0:
            return np.nan * np.zeros((c.shape[0], 3))
        else:
            return (c - self.point[axis]) * self.vector / self.vector[axis] + self.point

    def exit_points(self, bc):
        """
        Given a Trajectory and a BrainCoordinates object, computes the intersection of the
        trajectory with the brain coordinates bounding box
        :param bc: BrainCoordinate objects
        :return: np.ndarray 2 y 3 corresponding to exit points xyz coordinates
        """
        bounds = np.c_[bc.xlim, bc.ylim, bc.zlim]
        epoints = np.r_[self.eval_x(bc.xlim), self.eval_y(bc.ylim), self.eval_z(bc.zlim)]
        epoints = epoints[~np.all(np.isnan(epoints), axis=1)]
        ind = np.all(np.bitwise_and(bounds[0, :] <= epoints, epoints <= bounds[1, :]), axis=1)
        return epoints[ind, :]


@dataclass
class Insertion:
    """
    Defines an ephys probe insertion in 3D coordinate. IBL conventions.

    To instantiate, use the static methods: `Insertion.from_track` and `Insertion.from_dict`.
    """
    x: float
    y: float
    z: float
    phi: float
    theta: float
    depth: float
    label: str = ''
    beta: float = 0

    @staticmethod
    def from_track(xyzs, brain_atlas=None):
        """
        Define an insersion from one or more trajectory.

        Parameters
        ----------
        xyzs : numpy.array
             An n by 3 array xyz coordinates representing an insertion trajectory.
        brain_atlas : BrainAtlas
            A brain atlas instance, used to attain the point of entry.

        Returns
        -------
        Insertion
        """
        assert brain_atlas, 'Input argument brain_atlas must be defined'
        traj = Trajectory.fit(xyzs)
        # project the deepest point into the vector to get the tip coordinate
        tip = traj.project(xyzs[np.argmin(xyzs[:, 2]), :])
        # get intersection with the brain surface as an entry point
        entry = Insertion.get_brain_entry(traj, brain_atlas)
        # convert to spherical system to store the insertion
        depth, theta, phi = cart2sph(*(entry - tip))
        insertion_dict = {
            'x': entry[0], 'y': entry[1], 'z': entry[2], 'phi': phi, 'theta': theta, 'depth': depth
        }
        return Insertion(**insertion_dict)

    @staticmethod
    def from_dict(d, brain_atlas=None):
        """
        Constructs an Insertion object from the json information stored in probes.description file.

        Parameters
        ----------
        d : dict
            A dictionary containing at least the following keys {'x', 'y', 'z', 'phi', 'theta',
            'depth'}.  The depth and xyz coordinates must be in um.
        brain_atlas : BrainAtlas, default=None
            If provided, disregards the z coordinate and locks the insertion point to the z of the
            brain surface.

        Returns
        -------
        Insertion

        Examples
        --------
        >>> tri = {'x': 544.0, 'y': 1285.0, 'z': 0.0, 'phi': 0.0, 'theta': 5.0, 'depth': 4501.0}
        >>> ins = Insertion.from_dict(tri)
        """
        assert brain_atlas, 'Input argument brain_atlas must be defined'
        z = d['z'] / 1e6
        if not hasattr(brain_atlas, 'top'):
            brain_atlas.compute_surface()
        iy = brain_atlas.bc.y2i(d['y'] / 1e6)
        ix = brain_atlas.bc.x2i(d['x'] / 1e6)
        # Only use the brain surface value as z if it isn't NaN (this happens when the surface touches the edges
        # of the atlas volume
        if not np.isnan(brain_atlas.top[iy, ix]):
            z = brain_atlas.top[iy, ix]
        return Insertion(x=d['x'] / 1e6, y=d['y'] / 1e6, z=z,
                         phi=d['phi'], theta=d['theta'], depth=d['depth'] / 1e6,
                         beta=d.get('beta', 0), label=d.get('label', ''))

    @property
    def trajectory(self):
        """
        Gets the trajectory object matching insertion coordinates
        :return: atlas.Trajectory
        """
        return Trajectory.fit(self.xyz)

    @property
    def xyz(self):
        return np.c_[self.entry, self.tip].transpose()

    @property
    def entry(self):
        return np.array((self.x, self.y, self.z))

    @property
    def tip(self):
        return sph2cart(- self.depth, self.theta, self.phi) + np.array((self.x, self.y, self.z))

    @staticmethod
    def _get_surface_intersection(traj, brain_atlas, surface='top', mode='raise'):
        """
        Computes the intersection of a trajectory with either the top or the bottom surface of an atlas.

        Parameters
        ----------
        traj: iblatlas.atlas.Trajectory object
        brain_atlas: iblatlas.atlas.BrainAtlas (or descendant) object
        surface: str, optional (defaults to 'top') 'top' or 'bottom'
        mode: str, optional (defaults to 'raise') 'raise' or 'none': raise an error if no intersection
         with the brain surface is found otherwise returns None

        Returns
        -------
        xyz: np.array, 3 elements, x, y, z coordinates of the intersection point with the surface
             None if no intersection is found and mode is not set to 'raise'
        """
        brain_atlas.compute_surface()
        distance = traj.mindist(brain_atlas.srf_xyz)
        dist_sort = np.argsort(distance)
        # In some cases the nearest two intersection points are not the top and bottom of brain
        # So we find all intersection points that fall within one voxel and take the one with
        # highest dV to be entry and lowest dV to be exit
        idx_lim = np.sum(distance[dist_sort] * 1e6 < np.max(brain_atlas.res_um))
        if idx_lim == 0:  # no intersection found
            if mode == 'raise':
                raise ValueError('No intersection found with brain surface')
            else:
                return
        dist_lim = dist_sort[0:idx_lim]
        z_val = brain_atlas.srf_xyz[dist_lim, 2]
        if surface == 'top':
            ma = np.argmax(z_val)
            _xyz = brain_atlas.srf_xyz[dist_lim[ma], :]
            _ixyz = brain_atlas.bc.xyz2i(_xyz)
            _ixyz[brain_atlas.xyz2dims[2]] += 1
        elif surface == 'bottom':
            ma = np.argmin(z_val)
            _xyz = brain_atlas.srf_xyz[dist_lim[ma], :]
            _ixyz = brain_atlas.bc.xyz2i(_xyz)

        xyz = brain_atlas.bc.i2xyz(_ixyz.astype(float))

        return xyz

    @staticmethod
    def get_brain_exit(traj, brain_atlas, mode='raise'):
        """
        Given a Trajectory and a BrainAtlas object, computes the brain exit coordinate as the
        intersection of the trajectory and the brain surface (brain_atlas.surface)
        :param brain_atlas:
        :return: 3 element array x,y,z
        """
        # Find point where trajectory intersects with bottom of brain
        return Insertion._get_surface_intersection(traj, brain_atlas, surface='bottom', mode=mode)

    @staticmethod
    def get_brain_entry(traj, brain_atlas, mode='raise'):
        """
        Given a Trajectory and a BrainAtlas object, computes the brain entry coordinate as the
        intersection of the trajectory and the brain surface (brain_atlas.surface)
        :param brain_atlas:
        :return: 3 element array x,y,z
        """
        # Find point where trajectory intersects with top of brain
        return Insertion._get_surface_intersection(traj, brain_atlas, surface='top', mode=mode)


class AllenAtlas(BrainAtlas):
    """
    The Allan Common Coordinate Framework (CCF) brain atlas.

    Instantiates an atlas.BrainAtlas corresponding to the Allen CCF at the given resolution
    using the IBL Bregma and coordinate system.
    """

    """pathlib.PurePosixPath: The default relative path of the Allen atlas file."""
    atlas_rel_path = PurePosixPath('histology', 'ATLAS', 'Needles', 'Allen')

    """numpy.array: A diffusion weighted imaging (DWI) image volume.

    This average template brain was created from images of 1,675 young adult C57BL/6J mouse brains
    acquired using serial two-photon tomography.
    This volume has a C-ordered shape of (ap, ml, dv) and contains uint16
    values.
    """
    image = None

    """numpy.array: An annotation label volume.

    The Allen atlas label volume has with the shape (ap, ml, dv) and contains uint16 indices
    of the Allen CCF brain regions to which each voxel belongs.
    """
    label = None

    def __init__(self, res_um=25, scaling=(1, 1, 1), mock=False, hist_path=None):
        """
        Instantiates an atlas.BrainAtlas corresponding to the Allen CCF at the given resolution
        using the IBL Bregma and coordinate system.

        Parameters
        ----------
        res_um : {10, 25, 50} int
            The Atlas resolution in micrometres; one of 10, 25 or 50um.
        scaling : float, numpy.array
            Scale factor along ml, ap, dv for squeeze and stretch (default: [1, 1, 1]).
        mock : bool
            For testing purposes, return atlas object with image comprising zeros.
        hist_path : str, pathlib.Path
            The location of the image volume. May be a full file path or a directory.

        Examples
        --------
        Instantiate Atlas from a non-default location, in this case the cache_dir of an ONE instance.
        >>> target_dir = one.cache_dir / AllenAtlas.atlas_rel_path
        ... ba = AllenAtlas(hist_path=target_dir)
        """
        LUT_VERSION = 'v01'  # version 01 is the lateralized version
        regions = BrainRegions()
        xyz2dims = np.array([1, 0, 2])  # this is the c-contiguous ordering
        dims2xyz = np.array([1, 0, 2])
        # we use Bregma as the origin
        self.res_um = res_um
        ibregma = (ALLEN_CCF_LANDMARKS_MLAPDV_UM['bregma'] / self.res_um)
        dxyz = self.res_um * 1e-6 * np.array([1, -1, -1]) * scaling
        if mock:
            image, label = [np.zeros((528, 456, 320), dtype=np.int16) for _ in range(2)]
            label[:, :, 100:105] = 1327  # lookup index for retina, id 304325711 (no id 1327)
        else:
            # Hist path may be a full path to an existing image file, or a path to a directory
            cache_dir = Path(one.params.get(silent=True).CACHE_DIR)
            hist_path = Path(hist_path or cache_dir.joinpath(self.atlas_rel_path))
            if not hist_path.suffix:  # check if folder
                hist_path /= f'average_template_{res_um}.nrrd'
            # get the image volume
            if not hist_path.exists():
                hist_path = _download_atlas_allen(hist_path)
            # get the remapped label volume
            file_label = hist_path.with_name(f'annotation_{res_um}.nrrd')
            if not file_label.exists():
                file_label = _download_atlas_allen(file_label)
            file_label_remap = hist_path.with_name(f'annotation_{res_um}_lut_{LUT_VERSION}.npz')
            if not file_label_remap.exists():
                label = self._read_volume(file_label).astype(dtype=np.int32)
                _logger.info("Computing brain atlas annotations lookup table")
                # lateralize atlas: for this the regions of the left hemisphere have primary
                # keys opposite to to the normal ones
                lateral = np.zeros(label.shape[xyz2dims[0]])
                lateral[int(np.floor(ibregma[0]))] = 1
                lateral = np.sign(np.cumsum(lateral)[np.newaxis, :, np.newaxis] - 0.5)
                label = label * lateral.astype(np.int32)
                # the 10 um atlas is too big to fit in memory so work by chunks instead
                if res_um == 10:
                    first, ncols = (0, 10)
                    while True:
                        last = np.minimum(first + ncols, label.shape[-1])
                        _logger.info(f"Computing... {last} on {label.shape[-1]}")
                        _, im = ismember(label[:, :, first:last], regions.id)
                        label[:, :, first:last] = np.reshape(im, label[:, :, first:last].shape)
                        if last == label.shape[-1]:
                            break
                        first += ncols
                    label = label.astype(dtype=np.uint16)
                    _logger.info("Saving npz, this can take a long time")
                else:
                    _, im = ismember(label, regions.id)
                    label = np.reshape(im.astype(np.uint16), label.shape)
                np.savez_compressed(file_label_remap, label)
                _logger.info(f"Cached remapping file {file_label_remap} ...")
            # loads the files
            label = self._read_volume(file_label_remap)
            image = self._read_volume(hist_path)

        super().__init__(image, label, dxyz, regions, ibregma, dims2xyz=dims2xyz, xyz2dims=xyz2dims)

    @staticmethod
    def _read_volume(file_volume):
        if file_volume.suffix == '.nrrd':
            volume, _ = nrrd.read(file_volume, index_order='C')  # ml, dv, ap
            # we want the coronal slice to be the most contiguous
            volume = np.transpose(volume, (2, 0, 1))  # image[iap, iml, idv]
        elif file_volume.suffix == '.npz':
            volume = np.load(file_volume)['arr_0']
        return volume

    def xyz2ccf(self, xyz, ccf_order='mlapdv', mode='raise'):
        """
        Converts anatomical coordinates to CCF coordinates.

        Anatomical coordinates are in meters, relative to bregma, which CFF coordinates are
        assumed to be the volume indices multiplied by the spacing in micormeters.

        Parameters
        ----------
        xyz : numpy.array
            An N by 3 array of anatomical coordinates in meters, relative to bregma.
        ccf_order : {'mlapdv', 'apdvml'}, default='mlapdv'
            The order of the CCF coordinates returned. For IBL (the default) this is (ML, AP, DV),
            for Allen MCC vertices, this is (AP, DV, ML).
        mode : {'raise', 'clip', 'wrap'}, default='raise'
            How to behave if the coordinate lies outside of the volume: raise (default) will raise
            a ValueError; 'clip' will replace the index with the closest index inside the volume;
            'wrap' will return the index as is.

        Returns
        -------
        numpy.array
            Coordinates in CCF space (um, origin is the front left top corner of the data
        volume, order determined by ccf_order
        """
        ordre = self._ccf_order(ccf_order)
        ccf = self.bc.xyz2i(xyz, round=False, mode=mode) * float(self.res_um)
        return ccf[..., ordre]

    def ccf2xyz(self, ccf, ccf_order='mlapdv'):
        """
        Convert anatomical coordinates from CCF coordinates.

        Anatomical coordinates are in meters, relative to bregma, which CFF coordinates are
        assumed to be the volume indices multiplied by the spacing in micormeters.

        Parameters
        ----------
        ccf : numpy.array
            An N by 3 array of coordinates in CCF space (atlas volume indices * um resolution). The
            origin is the front left top corner of the data volume.
        ccf_order : {'mlapdv', 'apdvml'}, default='mlapdv'
            The order of the CCF coordinates given. For IBL (the default) this is (ML, AP, DV),
            for Allen MCC vertices, this is (AP, DV, ML).

        Returns
        -------
        numpy.array
            The MLAPDV coordinates in meters, relative to bregma.
        """
        ordre = self._ccf_order(ccf_order, reverse=True)
        return self.bc.i2xyz((ccf[..., ordre] / float(self.res_um)))

    @staticmethod
    def _ccf_order(ccf_order, reverse=False):
        """
        Returns the mapping to go from CCF coordinates order to the brain atlas xyz
        :param ccf_order: 'mlapdv' or 'apdvml'
        :param reverse: defaults to False.
            If False, returns from CCF to brain atlas
            If True, returns from brain atlas to CCF
        :return:
        """
        if ccf_order == 'mlapdv':
            return [0, 1, 2]
        elif ccf_order == 'apdvml':
            if reverse:
                return [2, 0, 1]
            else:
                return [1, 2, 0]
        else:
            ValueError("ccf_order needs to be either 'mlapdv' or 'apdvml'")

    def compute_regions_volume(self, cumsum=False):
        """
        Sums the number of voxels in the labels volume for each region.
        Then compute volumes for all of the levels of hierarchy in cubic mm.
        :param: cumsum: computes the cumulative sum of the volume as per the hierarchy (defaults to False)
        :return:
        """
        nr = self.regions.id.shape[0]
        count = np.bincount(self.label.flatten(), minlength=nr)
        if not cumsum:
            self.regions.volume = count * (self.res_um / 1e3) ** 3
        else:
            self.regions.compute_hierarchy()
            self.regions.volume = np.zeros_like(count)
            for i in np.arange(nr):
                if count[i] == 0:
                    continue
                self.regions.volume[np.unique(self.regions.hierarchy[:, i])] += count[i]
            self.regions.volume = self.regions.volume * (self.res_um / 1e3) ** 3


def NeedlesAtlas(*args, **kwargs) -> AllenAtlas:
    """
    Instantiates an atlas.BrainAtlas corresponding to the Allen CCF at the given resolution
    using the IBL Bregma and coordinate system. The Needles atlas defines a stretch along AP
    axis and a squeeze along the DV axis.

    Parameters
    ----------
    res_um : {10, 25, 50} int
        The Atlas resolution in micrometres; one of 10, 25 or 50um.
    **kwargs
        See AllenAtlas.

    Returns
    -------
    AllenAtlas
        An Allen atlas object with MRI atlas scaling applied.

    Notes
    -----
    The scaling was determined by manually transforming the DSURQE atlas [1]_ onto the Allen CCF.
    The DSURQE atlas is an MRI atlas acquired from 40 C57BL/6J mice post-mortem, with 40um
    isometric resolution.  The alignment was performed by Mayo Faulkner.
    The atlas data can be found `here <http://repo.mouseimaging.ca/repo/DSURQE_40micron_nifti/>`__.
    More information on the dataset and segmentation can be found
    `here <http://repo.mouseimaging.ca/repo/DSURQE_40micron/notes_on_DSURQE_atlas>`__.

    References
    ----------
    .. [1] Dorr AE, Lerch JP, Spring S, Kabani N, Henkelman RM (2008). High resolution
       three-dimensional brain atlas using an average magnetic resonance image of 40 adult C57Bl/6J
       mice. Neuroimage 42(1):60-9. [doi 10.1016/j.neuroimage.2008.03.037]
    """
    DV_SCALE = 0.952  # multiplicative factor on DV dimension, determined from MRI->CCF transform
    AP_SCALE = 1.087  # multiplicative factor on AP dimension
    kwargs['scaling'] = np.array([1, AP_SCALE, DV_SCALE])
    return AllenAtlas(*args, **kwargs)


def MRITorontoAtlas(*args, **kwargs):
    """
    The MRI Toronto brain atlas.

    Instantiates an atlas.BrainAtlas corresponding to the Allen CCF at the given resolution
    using the IBL Bregma and coordinate system. The MRI Toronto atlas defines a stretch along AP
    a squeeze along DV *and* a squeeze along ML. These are based on 12 p65 mice MRIs averaged [1]_.

    Parameters
    ----------
    res_um : {10, 25, 50} int
        The Atlas resolution in micrometres; one of 10, 25 or 50um.
    **kwargs
        See AllenAtlas.

    Returns
    -------
    AllenAtlas
        An Allen atlas object with MRI atlas scaling applied.

    References
    ----------
    .. [1] Qiu, LR, Fernandes, DJ, Szulc-Lerch, KU et al. (2018) Mouse MRI shows brain areas
       relatively larger in males emerge before those larger in females. Nat Commun 9, 2615.
       [doi 10.1038/s41467-018-04921-2]
    """
    ML_SCALE = 0.952
    DV_SCALE = 0.885  # multiplicative factor on DV dimension, determined from MRI->CCF transform
    AP_SCALE = 1.031  # multiplicative factor on AP dimension
    kwargs['scaling'] = np.array([ML_SCALE, AP_SCALE, DV_SCALE])
    return AllenAtlas(*args, **kwargs)


def _download_atlas_allen(target_file_image):
    """
    Download the Allen Atlas from the alleninstitute.org Website.

    Parameters
    ----------
    target_file_image : str, pathlib.Path
        The full target file path to which to download the file. The name of the image file name
         must be either `average_template_<res>.nrrd` or `annotation_<res>.nrrd`, where <res> is
         one of {10, 25, 50}.

    Returns
    -------
    pathlib.Path
        The full path to the downloaded file.

    Notes
    -----
    - © 2015 Allen Institute for Brain Science. Allen Mouse Brain Atlas (2015) with region annotations (2017).
    - Available from: http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/
    - See Allen Mouse Common Coordinate Framework Technical White Paper for details
      http://help.brain-map.org/download/attachments/8323525/Mouse_Common_Coordinate_Framework.pdf

    """
    (target_file_image := Path(target_file_image)).parent.mkdir(exist_ok=True, parents=True)
    ROOT_URL = 'http://download.alleninstitute.org/informatics-archive/'

    if target_file_image.name.split('_')[0] == 'average':
        url = f'{ROOT_URL}current-release/mouse_ccf/average_template/'
    elif target_file_image.name.split('_')[0] == 'annotation':
        url = f'{ROOT_URL}current-release/mouse_ccf/annotation/ccf_2017/'
    else:
        raise ValueError('Unrecognized file image')
    url += target_file_image.name

    return Path(http_download_file(url, target_dir=target_file_image.parent))


class FranklinPaxinosAtlas(BrainAtlas):

    """pathlib.PurePosixPath: The default relative path of the atlas file."""
    atlas_rel_path = PurePosixPath('histology', 'ATLAS', 'Needles', 'FranklinPaxinos')

    def __init__(self, res_um=(10, 100, 10), scaling=(1, 1, 1), mock=False, hist_path=None):
        """The Franklin & Paxinos brain atlas.

        Instantiates an atlas.BrainAtlas corresponding to the Franklin & Paxinos atlas [1]_ at the
        given resolution, matched to the Allen coordinate Framework [2]_ and using the IBL Bregma
        and coordinate system. The Franklin Paxisnos volume has resolution of 10um in ML and DV
        axis and 100 um in AP direction.

        Parameters
        ----------
        res_um : list, numpy.array
            The Atlas resolution in micometres in each dimension.
        scaling : float, numpy.array
            Scale factor along ml, ap, dv for squeeze and stretch (default: [1, 1, 1]).
        mock : bool
            For testing purposes, return atlas object with image comprising zeros.
        hist_path : str, pathlib.Path
            The location of the image volume. May be a full file path or a directory.

        Examples
        --------
        Instantiate Atlas from a non-default location, in this case the cache_dir of an ONE instance.
        >>> target_dir = one.cache_dir / AllenAtlas.atlas_rel_path
        ... ba = FranklinPaxinosAtlas(hist_path=target_dir)

        References
        ----------
        .. [1] Paxinos G, and Franklin KBJ (2012) The Mouse Brain in Stereotaxic Coordinates, 4th
        edition (Elsevier Academic Press)
        .. [2] Chon U et al (2019) Enhanced and unified anatomical labeling for a common mouse
        brain atlas [doi 10.1038/s41467-019-13057-w]
        """
        # TODO interpolate?
        LUT_VERSION = 'v01'  # version 01 is the lateralized version
        regions = FranklinPaxinosRegions()
        xyz2dims = np.array([1, 0, 2])  # this is the c-contiguous ordering
        dims2xyz = np.array([1, 0, 2])
        # we use Bregma as the origin
        self.res_um = np.asarray(res_um)
        ibregma = (PAXINOS_CCF_LANDMARKS_MLAPDV_UM['bregma'] / self.res_um)
        dxyz = self.res_um * 1e-6 * np.array([1, -1, -1]) * scaling
        if mock:
            image, label = [np.zeros((528, 456, 320), dtype=np.int16) for _ in range(2)]
            label[:, :, 100:105] = 1327  # lookup index for retina, id 304325711 (no id 1327)
        else:
            # Hist path may be a full path to an existing image file, or a path to a directory
            cache_dir = Path(one.params.get(silent=True).CACHE_DIR)
            hist_path = Path(hist_path or cache_dir.joinpath(self.atlas_rel_path))
            if not hist_path.suffix:  # check if folder
                hist_path /= f'average_template_{res_um[0]}_{res_um[1]}_{res_um[2]}.npz'

            # get the image volume
            if not hist_path.exists():
                hist_path.parent.mkdir(exist_ok=True, parents=True)
                aws.s3_download_file(f'atlas/FranklinPaxinos/{hist_path.name}', str(hist_path))
            # get the remapped label volume
            file_label = hist_path.with_name(f'annotation_{res_um[0]}_{res_um[1]}_{res_um[2]}.npz')
            if not file_label.exists():
                file_label.parent.mkdir(exist_ok=True, parents=True)
                aws.s3_download_file(f'atlas/FranklinPaxinos/{file_label.name}', str(file_label))

            file_label_remap = hist_path.with_name(f'annotation_{res_um[0]}_{res_um[1]}_{res_um[2]}_lut_{LUT_VERSION}.npz')

            if not file_label_remap.exists():
                label = self._read_volume(file_label).astype(dtype=np.int32)
                _logger.info("computing brain atlas annotations lookup table")
                # lateralize atlas: for this the regions of the left hemisphere have primary
                # keys opposite to to the normal ones
                lateral = np.zeros(label.shape[xyz2dims[0]])
                lateral[int(np.floor(ibregma[0]))] = 1
                lateral = np.sign(np.cumsum(lateral)[np.newaxis, :, np.newaxis] - 0.5)
                label = label * lateral.astype(np.int32)
                _, im = ismember(label, regions.id)
                label = np.reshape(im.astype(np.uint16), label.shape)
                np.savez_compressed(file_label_remap, label)
                _logger.info(f"Cached remapping file {file_label_remap} ...")
            # loads the files
            label = self._read_volume(file_label_remap)
            image = self._read_volume(hist_path)

        super().__init__(image, label, dxyz, regions, ibregma, dims2xyz=dims2xyz, xyz2dims=xyz2dims)

    @staticmethod
    def _read_volume(file_volume):
        """
        Loads an atlas image volume given a file path.

        Parameters
        ----------
        file_volume : pathlib.Path
            The file path of an image volume.  Currently supports .nrrd and .npz files.

        Returns
        -------
        numpy.array
            The loaded image volume with dimensions (ap, ml, dv).

        Raises
        ------
        ValueError
            Unknown file extension, expects either '.nrrd' or '.npz'.
        """
        if file_volume.suffix == '.nrrd':
            volume, _ = nrrd.read(file_volume, index_order='C')  # ml, dv, ap
            # we want the coronal slice to be the most contiguous
            volume = np.transpose(volume, (2, 0, 1))  # image[iap, iml, idv]
        elif file_volume.suffix == '.npz':
            volume = np.load(file_volume)['arr_0']
        else:
            raise ValueError(
                f'"{file_volume.suffix}" files not supported, must be either ".nrrd" or ".npz"')
        return volume


def get_bc(res_um=10):
    """
    Get BrainCoordinates object for Allen Atlas with bregma origin without reading in volumes

    Parameters
    ----------
    res_um : float or int
        The resolution of the brain coordinates object to return in um. Options are 50, 25 or 10

    Returns
    -------
    BrainCoordinates object with bregma origin
    """

    dims2xyz = np.array([1, 0, 2])
    scaling = np.array([1, 1, 1])
    sf = 50 / res_um
    im_shape = np.array([int(264 * sf), int(228 * sf), int(160 * sf)])
    iorigin = (ALLEN_CCF_LANDMARKS_MLAPDV_UM['bregma'] / res_um)
    dxyz = res_um * 1e-6 * np.array([1, -1, -1]) * scaling
    nxyz = np.array(im_shape)[dims2xyz]
    bc = BrainCoordinates(nxyz=nxyz, xyz0=(0, 0, 0), dxyz=dxyz)
    bc = BrainCoordinates(nxyz=nxyz, xyz0=-bc.i2xyz(iorigin), dxyz=dxyz)

    return bc
