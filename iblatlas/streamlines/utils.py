import pandas as pd
import numpy as np
from iblutil.numerical import ismember
import matplotlib.pyplot as plt
from iblatlas.atlas import get_bc, BrainAtlas, aws, AllenAtlas
from iblatlas.regions import BrainRegions


def _download_depth_files(file_name):
    """
    Download and return path to relevant file
    :param file_name:
    :return:
    """
    file_path = BrainAtlas._get_cache_dir().joinpath('depths', file_name)
    if not file_path.exists():
        file_path.parent.mkdir(exist_ok=True, parents=True)
        aws.s3_download_file(f'atlas/depths/{file_path.name}', file_path)

    return file_path


def xyz_to_depth(xyz, per=True, res_um=25):
    """
    For a given xyz coordinates return the depth from the surface of the cortex. The depth is returned
    as a percentage if per=True and in um if per=False. Note the lookup will only work for xyz cooordinates
    that are in the Isocortex of the Allen volume. If coordinates outside of this region are given then
    the depth is returned as nan.

    Parameters
    ----------
    xyz : numpy.array
        An (n, 3) array of Cartesian coordinates. The order is ML, AP, DV and coordinates should be given in meters
        relative to bregma.

    per : bool
        Whether to do the lookup in percentage from the surface of the cortex or depth in um from the surface of the cortex.

    res_um : float or int
        The resolution of the brain atlas to do the depth lookup

    Returns
    -------
        numpy.array
        The depths from the surface of the cortex for each cartesian coordinate. If the coordinate does not lie within
        the Isocortex, depth value returned is nan
    """

    ind_flat = np.load(_download_depth_files(f'depths_ind_{res_um}.npy'))
    depth_file = f'depths_per_{res_um}.npy' if per else f'depths_um_{res_um}.npy'
    depths = np.load(_download_depth_files(depth_file))
    bc = get_bc(res_um=res_um)

    ixyz = bc.xyz2i(xyz, mode='clip')
    iravel = np.ravel_multi_index((ixyz[:, 1], ixyz[:, 0], ixyz[:, 2]), (bc.ny, bc.nx, bc.nz))
    a, b = ismember(iravel, ind_flat)

    lookup_depths = np.full(iravel.shape, np.nan, dtype=np.float32)
    lookup_depths[a] = depths[b]

    return lookup_depths


def get_mask(volume='annotation', br=None):
    """
    Generate a mask to plot results onto

    Parameters:
    -----------
    volume : str, optional
        The type of volume to project. Options are:
        - 'image': Projects the anatomical image using max intensity.
        - 'annotation': Projects the labeled regions and maps them to RGB colors.
        - 'boundary': Projects labeled regions and extracts anatomical boundaries.
        Default is 'annotation'.

    br : BrainRegions, optional
        An instance of the BrainRegions If None, a default BrainRegions is initialized.

    Returns:
    --------
    img : np.ndarray
        The resulting 2D flatmap projection image, either grayscale, RGB, or binary mask,
        depending on the selected volume type.
    """

    br = br or BrainRegions()
    if volume == 'image':
        img = np.load(_download_depth_files('dorsal_image.npy'))
    elif volume == 'annotation':
        img = np.load(_download_depth_files('dorsal_annotation.npy'))
        img = br.rgb[img]
    elif volume == 'boundary':
        img = np.load(_download_depth_files('dorsal_annotation.npy'))
        img = AllenAtlas.compute_boundaries(img)

    return img


def validate_aggr(aggr: str) -> None:
    """
    Validates if the provided aggregation type is valid.

    Parameters:
    ----------
    aggr : str
        The aggregation method to validate (e.g., 'mean', 'sum').

    Raises:
    ------
    AssertionError
        If the aggregation type is not one of the allowed values.
    """
    poss_aggrs = ['sum', 'count', 'mean', 'std', 'median', 'min', 'max', 'first', 'last']
    assert aggr in poss_aggrs, f"Aggregation must be one of {poss_aggrs}."


def project_volume_onto_flatmap(vol: np.ndarray, res_um: int = 25, aggr: str = 'mean', plot: bool = True,
                                cmap: str = 'viridis', clevels: tuple = None, ax: plt.Axes = None) -> np.ndarray:
    """
    Projects a 3D volume onto a 2D flatmap by aggregating values along streamline paths.

    Parameters:
    ----------
    vol : np.ndarray
        A 3D numpy array representing the volume data to be projected.

    res_um : int
        The resolution of the volume. Must be one of 10, 25 or 50.

    aggr : str, optional
        The aggregation method ('sum', 'count', 'mean', etc.), default is 'mean'.

    plot : bool, optional
        Whether to plot the resulting projection, default is True.

    cmap : str, optional
        The colormap to use for the plot, default is 'viridis'.

    clevels : tuple, optional
        The color limits to use for the plot, default is None.

    ax : matplotlib.axes.Axes, optional
        The axes on which to plot, default is None.

    Returns:
    -------
    np.ndarray
        The projected 2D array onto the flatmap.
    matplotlib.figure.Figure
        Matplotlib figure object if plot=True, otherwise None.
    matplotlib.axes.Axes
        Matplotlib axes object if plot=True, otherwise None.
    """
    bc = get_bc(res_um)
    assert vol.shape == (bc.ny, bc.nx, bc.nz), f"Volume does not have expected shape of {(bc.ny, bc.nx, bc.nz)}"

    # Validate the aggregation type
    validate_aggr(aggr)

    # Load the streamline paths
    path_df = pd.read_parquet(_download_depth_files(f'paths_{res_um}.pqt'))

    # Extract values from the volume using the path lookup
    path_df['vals'] = vol.flat[path_df['lookup'].values]

    # Aggregate the values along each path
    flat_df = path_df.groupby('paths').vals.agg(aggr)

    # Project the aggregated data onto the flatmap
    return _project_onto_flatmap(flat_df, plot=plot, cmap=cmap, clevels=clevels, ax=ax)


def project_points_onto_flatmap(xyz: np.ndarray, values: np.ndarray, res_um: int = 25, aggr: str = 'mean', plot: bool = True,
                                cmap: str = 'viridis', clevels: tuple = None, ax: plt.Axes = None) -> np.ndarray:
    """
    Projects 3D points with associated values onto a 2D flatmap.

    Parameters:
    ----------
    xyz : np.ndarray
        An array containing xyz coordinates of the points to be projected. xyz values should be given in metres

    values : np.ndarray
        A 1D array of values to associate with the points.

    res : int
        The resolution to load the corresponding streamline paths.

    aggr : str, optional
        The aggregation method ('sum', 'count', 'mean', etc.), default is 'mean'.

    plot : bool, optional
        Whether to plot the resulting projection, default is True.

    cmap : str, optional
        The colormap to use for the plot, default is 'viridis'.

    clevels : tuple, optional
        The color limits to use for the plot, default is None.

    ax : matplotlib.axes.Axes, optional
        The axes on which to plot, default is None.

    Returns:
    -------
    np.ndarray
        The projected 2D array onto the flatmap.
    matplotlib.figure.Figure
        Matplotlib figure object if plot=True, otherwise None.
    matplotlib.axes.Axes
        Matplotlib axes object if plot=True, otherwise None.
    """
    # Ensure that xyz and values have matching dimensions
    assert xyz.shape[0] == values.size, "xyz must have the same number of rows as values."

    # Validate the aggregation type
    validate_aggr(aggr)

    # Get the boundary coordinates for the given resolution
    bc = get_bc(res_um)

    # Convert coordinates xyz to indices in volume
    ixyz = bc.xyz2i(xyz, mode='clip')

    # Create DataFrame of values and their corresponding flattened indices
    val_df = pd.DataFrame()
    val_df['vals'] = values
    val_df['lookup'] = np.ravel_multi_index((ixyz[:, 1], ixyz[:, 0], ixyz[:, 2]), (bc.ny, bc.nx, bc.nz))

    # Remove entries with invalid lookups
    val_df = val_df[val_df['lookup'] != 0]

    # Load streamline paths
    path_df = pd.read_parquet(_download_depth_files(f'paths_{res_um}.pqt'))

    # Restrict dataframes to overlapping locations
    vals_in_paths, _ = ismember(val_df['lookup'].values, path_df['lookup'].values)
    val_df = val_df[vals_in_paths]

    # Keep only paths that match the val lookups
    paths_in_vals, _ = ismember(path_df['lookup'].values, val_df['lookup'].values)
    path_df = path_df[paths_in_vals]

    # Merge path data with values and aggregate by path
    flat_df = path_df.merge(val_df, on='lookup', how='left').groupby('paths').vals.agg(aggr)

    # Project the aggregated data onto the flatmap
    return _project_onto_flatmap(flat_df, plot=plot, cmap=cmap, clevels=clevels, ax=ax)


def _project_onto_flatmap(flat_df: pd.Series, plot: bool = True, cmap: str = 'viridis',
                          clevels: tuple = None, ax: plt.Axes = None) -> np.ndarray:
    """
    Function to project aggregated data onto a 2D flatmap.

    Parameters:
    ----------
    flat_df : pd.Series
        The data to project onto the flatmap.

    plot : bool, optional
        Whether to plot the resulting projection, default is True.

    cmap : str, optional
        The colormap to use for the plot, default is 'viridis'.

    clevels : tuple, optional
        The color limits to use for the plot, default is None.

    ax : matplotlib.axes.Axes, optional
        The axes on which to plot, default is None.

    Returns:
    -------
    np.ndarray
        The projected 2D array onto the flatmap.
    plt.figure
        Matplotlib figure object if plot=True, otherwise None.
    plt.axes
        Matplotlib axes object if plot=True, otherwise None.
    """
    # Load the flatmap
    flatmap = np.load(_download_depth_files('dorsal_flatmap.npy'))

    # Find the indices in the flatmap corresponding to the projected data
    _, b = ismember(flat_df.index, flatmap.flat)

    # Initialize the projection array with zeros
    proj = np.zeros(np.prod(flatmap.shape))

    # Assign the values from flat_df to the projection array
    proj[b] = flat_df.values

    # Reshape the 1D projection into the 2D flatmap shape
    proj = proj.reshape(flatmap.shape)

    # Plot the result if requested
    if plot:
        if ax:
            fig = ax.get_figure()
        else:
            fig, ax = plt.subplots()

        if clevels is None:
            clevels = (np.nanmin(proj), np.nanmax(proj))

        ax.imshow(proj, cmap=cmap, vmin=clevels[0], vmax=clevels[1])

        return proj, fig, ax
    else:
        return proj
