"""Regenerate the explanation-page figures in ``docs/explanation/images/``.

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

from iblatlas.atlas import AllenAtlas, MRITorontoAtlas, NeedlesAtlas
from iblatlas.genomics import agea

IMAGES_DIR = Path(__file__).resolve().parents[1] / 'explanation' / 'images'
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


def main(out_dir: Path = IMAGES_DIR) -> None:
    """Write the three PNGs into ``out_dir``, overwriting any that exist."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, make in (('pitch-overlay', pitch_overlay), ('scaling-overlay', scaling_overlay),
                       ('agea-raw-vs-processed', agea_raw_vs_processed)):
        fig = make()
        fig.savefig(out_dir / f'{name}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)


if __name__ == '__main__':
    main()
    embed_brain_in_frames_svg()
