"""Regenerate the explanation-page figures in ``docs/explanation/images/`` and the landing-page
images in ``docs/images/landing/``.

Run by hand, not by the site build: ``python docs/figures/make_figures.py``. It downloads atlas
data on first use. ``frames-allen-vs-ibl.svg`` is drawn by hand; this script only refreshes the
sagittal template slice embedded in it.
"""
import base64
import io
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc
from PIL import Image

from iblatlas.atlas import AllenAtlas, MRITorontoAtlas, NeedlesAtlas
from iblatlas.genomics import agea, lipids, merfish
from iblatlas.plots import (plot_points_on_slice, plot_scalar_on_slice, plot_swanson,
                            plot_swanson_vector)

IMAGES_DIR = Path(__file__).resolve().parents[1] / 'explanation' / 'images'
LANDING_DIR = Path(__file__).resolve().parents[1] / 'images' / 'landing'
GUI_SCREENSHOT = Path(__file__).resolve().parents[1] / 'how-to' / 'images' / 'gui-hover.png'
LANDING_SIZE = (5, 3.75)  # inches; every landing image shares this 4:3 canvas so the cards line up
TILT_DEG = 5  # the pitch stated on the coordinates page; lambda is not defined in the package
OUTLINE_COLOURS = {'AllenAtlas': 'k', 'NeedlesAtlas': 'tab:blue', 'MRITorontoAtlas': 'tab:red'}
# Slices with no missing voxels: Gfap's raw volume is left/right asymmetric here, and Calb1's
# overall level changes in the processed volume.
AGEA_EXAMPLES = (('Gfap', -0.0016), ('Calb1', -0.0004))
ATLAS_FACTORIES = {'AllenAtlas': AllenAtlas, 'NeedlesAtlas': NeedlesAtlas,
                   'MRITorontoAtlas': MRITorontoAtlas}


def pitch_overlay(res_um: int = 25) -> plt.Figure:
    """
    The midline sagittal slice with the atlas horizontal and the flat-skull horizontal at bregma.

    The 5° angle between the two lines is marked with an arc. The package defines bregma but not
    lambda, so the flat-skull line is drawn from the documented angle (``TILT_DEG``) rather than
    through two landmarks.

    Parameters
    ----------
    res_um : int
        Atlas resolution in micrometres.

    Returns
    -------
    matplotlib.figure.Figure
        One axis, in micrometres relative to bregma, anterior on the left as in
        ``BrainAtlas.plot_sslice``.
    """
    atlas = AllenAtlas(res_um)
    fig, ax = plt.subplots(figsize=(8, 4.2), layout='constrained')
    atlas.plot_sslice(0, ax=ax, cmap='gray')
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    ap = np.array(xlim)
    ax.plot(ap, 0 * ap, color='tab:blue', lw=2, label='atlas horizontal')
    ax.plot(ap, ap * np.tan(np.deg2rad(TILT_DEG)), color='tab:red', lw=2,
            label='flat-skull horizontal (drawn from the angle)')
    ax.plot(0, 0, 'o', color='yellow', markeredgecolor='k', label='bregma')
    radius = 3200  # µm, the arc between the two lines, on the anterior side
    wedge_ap = np.array([0, xlim[0]])  # fill between the lines on the anterior side of bregma
    ax.fill_between(wedge_ap, 0, wedge_ap * np.tan(np.deg2rad(TILT_DEG)), color='gold', alpha=0.6)
    ax.add_patch(Arc((0, 0), 2 * radius, 2 * radius, theta1=0, theta2=TILT_DEG, color='k', lw=1.5))
    ax.text(radius + 200, 700, f'{TILT_DEG}°', fontsize=14, fontweight='bold', ha='center',
            bbox={'boxstyle': 'round', 'fc': 'white', 'ec': 'none', 'alpha': 0.85})
    ax.annotate('', xy=(radius, radius * np.tan(np.deg2rad(TILT_DEG / 2))), xytext=(radius + 200, 620),
                arrowprops={'arrowstyle': '-', 'color': 'k'})
    ax.set(xlim=xlim, ylim=(ylim[0], 1400), xlabel='AP (µm from bregma)', ylabel='DV (µm from bregma)',
           title=f'The atlas is pitched about {TILT_DEG}° nose-down relative to a flat skull')
    ax.legend(loc='lower left', fontsize=8)
    return fig


def scaling_overlay(res_um: int = 25) -> plt.Figure:
    """
    Coronal and sagittal brain outlines of the Allen, Needles and MRI Toronto atlases, overlaid.

    Parameters
    ----------
    res_um : int
        Atlas resolution in micrometres.

    Returns
    -------
    matplotlib.figure.Figure
        Two axes (coronal, sagittal) in metres relative to bregma.
    """
    fig, (cax, sax) = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    for name, factory in ATLAS_FACTORIES.items():
        atlas = factory(res_um)
        for ax, axis, coord in ((cax, 1, -0.0020), (sax, 0, 0.0)):
            brain = atlas.slice(coord, axis=axis, volume='rindex').T != 0  # rows = DV, as imshow
            extent = atlas.extent(axis) * 1e-6
            ax.contour(brain.astype(float), levels=[0.5], colors=OUTLINE_COLOURS[name],
                       extent=extent, origin='upper')
        sax.plot([], [], color=OUTLINE_COLOURS[name], label=name)
        cax.plot([], [], color=OUTLINE_COLOURS[name], label=name)
    cax.set(title='Coronal, AP = -2 mm', xlabel='ML (m)', ylabel='DV (m)')
    sax.set(title='Sagittal, ML = 0', xlabel='AP (m)', ylabel='DV (m)')
    sax.invert_xaxis()  # anterior on the left, as in plot_sslice
    cax.legend(fontsize=8)
    return fig


def mappings_coronal(res_um: int = 25, ap_coordinate: float = -0.0016) -> plt.Figure:
    """
    One coronal slice drawn with the Allen, Beryl and Cosmos mappings, beside the Swanson flat map.

    Each voxel takes the Allen colour of its mapped target, so the number of distinct patches
    falls as the mapping gets coarser. Swanson regions have no 3D position, so that mapping is
    shown as the flat map rather than as a slice.

    Parameters
    ----------
    res_um : int
        Atlas resolution in micrometres.
    ap_coordinate : float
        Anterior-posterior position of the slice, in metres from bregma.

    Returns
    -------
    matplotlib.figure.Figure
        Four axes in a row: three coronal slices, then the Swanson flat map.
    """
    atlas = AllenAtlas(res_um)
    fig, axs = plt.subplots(1, 4, figsize=(14, 3.6), layout='constrained',
                            width_ratios=(1, 1, 1, 1.3))
    for ax, name in zip(axs, ('Allen', 'Beryl', 'Cosmos')):
        indices = atlas.slice(ap_coordinate, axis=1, volume='rindex', mapping=name)
        n_targets = np.unique(indices[indices != 0]).size  # index 0 is the void outside the brain
        atlas.plot_cslice(ap_coordinate, volume='annotation', mapping=name, ax=ax)
        ax.set(title=f'{name}: {n_targets} regions on this slice', xticks=[], yticks=[],
               xlabel='', ylabel='')
    plot_swanson(ax=axs[3], br=atlas.regions)
    axs[3].set(title='Swanson: the flat map')
    axs[3].set_axis_off()
    for ax in axs:
        ax.set_anchor('N')  # equal-aspect axes differ in height; pin the tops so titles line up
    return fig


def agea_raw_vs_processed(examples: tuple = AGEA_EXAMPLES) -> plt.Figure:
    """
    Coronal slices of two genes' raw and processed AGEA volumes, side by side.

    The slices are fully measured in the raw volume, so the figure shows the everyday
    difference: the processed volume is ML-symmetric and its overall level can change. Each row
    shares a colour scale.

    Parameters
    ----------
    examples : tuple of (str, float)
        ``(gene symbol, AP coordinate in metres relative to bregma)``, one per row. If several
        experiments exist for a gene, the first is used.

    Returns
    -------
    matplotlib.figure.Figure
        Axes in row-major order: for each example, ``agea.load()`` then
        ``agea.load(label='processed')``.
    """
    fig, axs = plt.subplots(len(examples), 2, figsize=(8, 3.2 * len(examples)), layout='constrained')
    loaded = {label: agea.load(label=label) for label in ('', 'processed')}
    for row, (gene, ap_m) in zip(axs, examples):
        volumes = {label: np.asarray(vols[df.index[df.gene == gene][0]])
                   for label, (df, vols, _) in loaded.items()}
        in_brain = loaded['processed'][2].label != 0
        vmax = np.percentile(volumes['processed'][in_brain], 99.5)
        for ax, (label, volume) in zip(row, volumes.items()):
            loaded[label][2].plot_cslice(ap_m, ax=ax, volume=volume, cmap='viridis', vmin=0, vmax=vmax)
            title = 'agea.load()' if label == '' else f"agea.load(label='{label}')"
            ax.set_title(f'{title}: {gene}, AP {ap_m * 1e3:.1f} mm', fontsize=9)
    return fig


def sagittal_template_png(res_um: int = 25) -> bytes:
    """
    The midline sagittal template slice as an RGBA PNG, for the frames diagram.

    Anterior is on the left and dorsal at the top, as in ``BrainAtlas.plot_sslice``. The image
    spans the whole atlas volume, so its edges are the edges of the volume. Voxels outside the
    brain are transparent and the template is rescaled to mid-greys so it reads on a light or a
    dark page.

    Parameters
    ----------
    res_um : int
        Atlas resolution in micrometres. The image has one pixel per voxel.

    Returns
    -------
    bytes
        PNG data, DV rows by AP columns.
    """
    atlas = AllenAtlas(res_um)
    template = atlas.slice(0, axis=0, volume='image').T.astype(float)
    in_brain = atlas.slice(0, axis=0, volume='rindex').T != 0
    grey = 0.3 + 0.55 * (template - template.min()) / np.ptp(template)
    rgba = np.dstack([grey, grey, grey, in_brain.astype(float)])
    buffer = io.BytesIO()
    plt.imsave(buffer, rgba, format='png')
    return buffer.getvalue()


def embed_brain_in_frames_svg(svg_path: Path = IMAGES_DIR / 'frames-allen-vs-ibl.svg') -> None:
    """
    Refresh the sagittal template slice embedded in the hand-drawn frames diagram.

    Rewrites the ``href`` of the ``<image id="brain">`` element as a data URI, because an SVG
    shown through an ``<img>`` tag cannot load an external file. Nothing else in the SVG changes.

    Parameters
    ----------
    svg_path : Path
        The diagram, which must contain exactly one ``<image id="brain" ...>`` element.
    """
    uri = 'data:image/png;base64,' + base64.b64encode(sagittal_template_png()).decode()
    text = svg_path.read_text()
    text, count = re.subn(r'(<image id="brain"[^>]*?href=")[^"]*(")', lambda m: m[1] + uri + m[2], text)
    if count != 1:
        raise ValueError(f'{svg_path} must contain exactly one <image id="brain" href="...">')
    svg_path.write_text(text)


def _landing_axes(nrows: int = 1, height_ratios=None) -> tuple:
    """A figure on the shared landing canvas with ``nrows`` axes stacked, axes switched off."""
    fig, axs = plt.subplots(nrows, 1, figsize=LANDING_SIZE, layout='constrained',
                            height_ratios=height_ratios, squeeze=False)
    axs = axs[:, 0]
    for ax in axs:
        ax.axis('off')
    return fig, axs


def _clear_void(ax: plt.Axes, atlas: AllenAtlas, ap_coordinate: float) -> None:
    """Make the voxels outside the brain transparent in every image of a coronal-slice axis."""
    void = np.moveaxis(atlas.slice(ap_coordinate, axis=1, volume='rindex'), 0, 1) == 0
    for image in ax.images:
        rgba = image.to_rgba(image.get_array())
        rgba[void, 3] = 0
        image.set_data(rgba)


def landing_region_values(atlas: AllenAtlas) -> plt.Figure:
    """A coronal slice over the Swanson flat map, both coloured by one value per Beryl region.

    The values are random: the image illustrates the plotting calls, not a result.
    """
    br = atlas.regions
    acronyms = br.acronym[np.unique(br.mappings['Beryl'])]
    acronyms = acronyms[~np.isin(acronyms, ['void', 'root'])]
    values = np.random.default_rng(1).random(acronyms.size)
    fig, axs = _landing_axes(2, height_ratios=(1, 0.75))
    plot_scalar_on_slice(acronyms, values, coord=-2000, slice='coronal', mapping='Beryl',
                         hemisphere='both', background='image', cmap='Reds', brain_atlas=atlas,
                         ax=axs[0])
    _clear_void(axs[0], atlas, -0.002)
    plot_swanson_vector(acronyms, values, cmap='Reds', empty_color='silver', br=br, ax=axs[1])
    return fig


def landing_volumes() -> plt.Figure:
    """A coronal slice of the most reliable lipid volume, which is one value per voxel."""
    volume, df_lipids, atlas = lipids.load_volume(min_reliability=0.6)
    fig, axs = _landing_axes()
    atlas.plot_cslice(-0.002, ax=axs[0], volume=volume[df_lipids['reliability'].idxmax()],
                      cmap='viridis')
    return fig


def landing_probes(atlas: AllenAtlas) -> plt.Figure:
    """A synthetic probe track, 120 channels coloured by value, on a sagittal slice.

    The track is a straight line and the values a sine wave: the image illustrates the plotting
    call, not a recording.
    """
    t = np.linspace(0, 1, 120)
    xyz = np.c_[np.full(120, -0.0010), 0.0005 - 0.0006 * t, -0.0005 - 0.0035 * t]
    fig, axs = _landing_axes()
    plot_points_on_slice(xyz, values=1 + np.sin(t * 14) ** 2 * 4, coord=-1000, slice='sagittal',
                         mapping='Allen', background='boundary', cmap='Reds', aggr='mean',
                         fwhm=150, brain_atlas=atlas, ax=axs[0])
    return fig


def landing_mappings(atlas: AllenAtlas) -> plt.Figure:
    """One coronal slice with the Allen mapping above the Cosmos mapping."""
    fig, axs = _landing_axes(2)
    for ax, name in zip(axs, ('Allen', 'Cosmos')):
        atlas.plot_cslice(-0.0016, volume='annotation', mapping=name, ax=ax)
        _clear_void(ax, atlas, -0.0016)
    return fig


def landing_external_data() -> plt.Figure:
    """One MERFISH section, each cell coloured by its class."""
    df_cells, df_classes, *_ = merfish.load()
    in_section = df_cells['brain_section_label'] == 'Zhuang-ABCA-1.085'
    xyz = df_cells.loc[in_section, ['x', 'y', 'z']].values
    rgb = merfish.int2rgb(df_classes.loc[df_cells.loc[in_section, 'class'], 'class_rgba'])
    fig, axs = _landing_axes()
    axs[0].scatter(xyz[:, 0], xyz[:, 2], c=rgb, s=1)
    axs[0].set_aspect('equal', 'box')
    return fig


def landing_gui(out_dir: Path) -> None:
    """The GUI hover screenshot from the GUI how-to, cropped to the slice and its tooltip."""
    im = Image.open(GUI_SCREENSHOT).convert('RGB')
    width, height = im.size
    height = int(height * 0.93)  # drop the status line under the slice
    left = width - round(height * LANDING_SIZE[0] / LANDING_SIZE[1])  # keep the right edge, 4:3
    im.crop((left, 0, width, height)).save(out_dir / 'gui.png')


def landing_images(out_dir: Path = LANDING_DIR) -> None:
    """Write the six landing-page card images into ``out_dir``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    atlas = AllenAtlas(25)
    figures = {'region-values': landing_region_values(atlas), 'volumes': landing_volumes(),
               'probes': landing_probes(atlas), 'mappings': landing_mappings(atlas),
               'external-data': landing_external_data()}
    for name, fig in figures.items():
        fig.savefig(out_dir / f'{name}.png', dpi=110, facecolor='white')
        plt.close(fig)
    landing_gui(out_dir)


def main(out_dir: Path = IMAGES_DIR) -> None:
    """Write the four PNGs into ``out_dir``, overwriting any that exist."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, make in (('pitch-overlay', pitch_overlay), ('scaling-overlay', scaling_overlay),
                       ('agea-raw-vs-processed', agea_raw_vs_processed),
                       ('mappings-coronal', mappings_coronal)):
        fig = make()
        fig.savefig(out_dir / f'{name}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
    landing_images()


if __name__ == '__main__':
    main()
    embed_brain_in_frames_svg()
