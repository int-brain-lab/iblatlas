"""
TopView is the main Widget with the related ControllerTopView Class
There are several SliceView windows (sagittal, coronal, possibly tilted etc...) that each have
a SliceController object
The underlying data model object is an iblatlas.atlas.AllenAtlas object

    TopView(QMainWindow)
    ControllerTopView(PgImageController)

    SliceView(QWidget)
    SliceController(PgImageController)

"""
import sys

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from qtpy import QtWidgets, uic, QtCore
from qtpy.QtGui import QTransform
import pyqtgraph as pg
import matplotlib

from iblatlas.atlas import AllenAtlas

from ibllib.misc import qt  # FIXME: remove ibllib dependency
from iblatlas.gui.braintree import BrainTree
from iblutil.numerical import ismember


class TopView(QtWidgets.QMainWindow):
    """
    Main Window of the application.
    This is a top view of the brain with 2 movable lines allowing to select sagittal and coronal
    slices.
    """
    @staticmethod
    def _instances():
        app = QtWidgets.QApplication.instance()
        return [w for w in app.topLevelWidgets() if isinstance(w, TopView)]

    @staticmethod
    def _get_or_create(title=None, **kwargs):
        av = next(filter(lambda e: e.isVisible() and e.windowTitle() == title,
                         TopView._instances()), None)
        if av is None:
            av = TopView(**kwargs)
            av.setWindowTitle(title)
        return av

    def __init__(self, **kwargs):
        super(TopView, self).__init__()
        self.ctrl = ControllerTopView(self, **kwargs)
        self.ctrl.image_layers = {'top': ImageLayer()}
        uic.loadUi(Path(__file__).parent.joinpath('topview.ui'), self)
        if (geom := self.ctrl.settings.value('geometry', None)) is not None:
            self.restoreGeometry(geom)
        self.plotItem_topview.setAspectLocked(True)
        self.plotItem_topview.addItem(self.ctrl.imageItem)
        # setup one horizontal and one vertical line that can be moved
        for line in self.ctrl.lines_coronal:
            line.sigDragged.connect(self._refresh_coronal)
        for line in self.ctrl.lines_sagittal:  # sigPositionChangeFinished
            line.sigDragged.connect(self._refresh_sagittal)
        for line in self.ctrl.lines_horizontal:  # sigPositionChangeFinished
            line.sigDragged.connect(self._refresh_horizontal)
        # set the horizontal slice start value in the middle of the volume
        self.ctrl.lines_horizontal[0].setValue(-.004)
        self._refresh_horizontal(self.ctrl.lines_horizontal[0])
        self.plotItem_topview.addItem(self.ctrl.line_coronal)
        self.plotItem_topview.addItem(self.ctrl.line_sagittal)
        # connect signals and slots: mouse moved
        s = self.plotItem_topview.getViewBox().scene()
        self.proxy = pg.SignalProxy(s.sigMouseMoved, rateLimit=60, slot=self.mouseMoveEvent)
        # combobox for the atlas remapping choices
        self.comboBox_mappings.addItems(self.ctrl.atlas.regions.mappings.keys())
        self.comboBox_mappings.currentIndexChanged.connect(self._refresh)
        # slider for transparency between image and labels
        self.slider_alpha.sliderMoved.connect(self.slider_alpha_move)
        self.ctrl.set_top()
        self.ctrl.fig_brain_tree.signal_region_selected.connect(self.on_brain_tree_selection)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)

    def focusInEvent(self, event):
        print('focusInEvent')
        for fig in self.ctrl.figures.values():
            fig.setWindowState(fig.windowState() & ~QtCore.Qt.WindowMinimized | QtCore.Qt.WindowActive)
            fig.activateWindow()

    def closeEvent(self, event):
        super(TopView, self).closeEvent(event)
        self.ctrl.settings.setValue('geometry', self.saveGeometry())
        for k, fig in self.ctrl.figures.items():
            self.ctrl.settings.setValue(f'geometry_{k}', fig.saveGeometry())
            fig.destroy()
        self.destroy()
        QtWidgets.QApplication.instance().quit()

    def add_scatter_feature(self, data):
        self.ctrl.scatter_data = data / 1e6
        self.ctrl.scatter_data_ind = self.ctrl.atlas.bc.xyz2i(self.ctrl.scatter_data)
        self.ctrl.fig_coronal.add_scatter()
        self.ctrl.fig_sagittal.add_scatter()
        self.line_coronal.sigDragged.connect(
            lambda: self.ctrl.set_scatter(self.ctrl.fig_coronal, self.line_coronal.value()))
        self.line_sagittal.sigDragged.connect(
            lambda: self.ctrl.set_scatter(self.ctrl.fig_sagittal, self.line_sagittal.value()))
        self.ctrl.set_scatter(self.ctrl.fig_coronal)
        self.ctrl.set_scatter(self.ctrl.fig_sagittal)

    def add_image_layer(self, **kwargs):
        """
        :param pg_kwargs: pyqtgraph setImage arguments: {'levels': None, 'lut': None,
        'opacity': 1.0}
        :param slice_kwargs: iblatlas.atlas.slice arguments: {'volume': 'image', 'mode': 'clip'}
        :return:
        """
        self.ctrl.fig_sagittal.add_image_layer(**kwargs)
        self.ctrl.fig_coronal.add_image_layer(**kwargs)

    def add_regions_feature(self, values, cmap, opacity=1.0):
        self.ctrl.values = values
        # creat cmap look up table
        colormap = matplotlib.cm.get_cmap(cmap)
        colormap._init()
        lut = (colormap._lut * 255).view(np.ndarray)
        lut = np.insert(lut, 0, [0, 0, 0, 0], axis=0)
        self.add_image_layer(pg_kwargs={'lut': lut, 'opacity': opacity}, slice_kwargs={
            'volume': 'value', 'region_values': values, 'mode': 'clip'})

    def slider_alpha_move(self):
        annotation_alpha = self.slider_alpha.value() / 100
        for _, fslice in self.ctrl.slices.items():
            fslice.ctrl.image_layers['image'].pg_kwargs['opacity'] = 1 - annotation_alpha
            fslice.ctrl.image_layers['annotation'].pg_kwargs['opacity'] = annotation_alpha
        self._refresh()

    def mouseMoveEvent(self, scenepos):
        if isinstance(scenepos, tuple):
            scenepos = scenepos[0]
        else:
            return
        pass
        # qpoint = self.imageItem.mapFromScene(scenepos)

    @QtCore.Slot(int)
    def on_brain_tree_selection(self, rid):
        self.ctrl.highlight_region = rid
        self._refresh()
        self.ctrl.set_top()

    def _refresh(self):
        self._refresh_sagittal()
        self._refresh_coronal()
        self._refresh_horizontal()

    def _refresh_coronal(self, line=None):
        line = self.ctrl.line_coronal if line is None else line
        self.ctrl.set_slice(self.ctrl.fig_coronal, val := line.value(),
                            mapping=self.comboBox_mappings.currentText())
        for line in self.ctrl.lines_coronal:
            line.setValue(val)

    def _refresh_sagittal(self, line=None):
        line = self.ctrl.line_sagittal if line is None else line
        self.ctrl.set_slice(self.ctrl.fig_sagittal, val := line.value(),
                            mapping=self.comboBox_mappings.currentText())
        for line in self.ctrl.lines_sagittal:
            line.setValue(val)

    def _refresh_horizontal(self, line=None):
        line = self.ctrl.lines_horizontal[0] if line is None else line
        self.ctrl.set_slice(self.ctrl.fig_horizontal, val := line.value(),
                            mapping=self.comboBox_mappings.currentText())
        for line in self.ctrl.lines_horizontal:
            line.setValue(val)

    def set_volume(self, volume: np.ndarray, colormap: str = 'magma', levels=None):
        self.ctrl.volume = volume
        cmap = pg.colormap.get(colormap)
        self.ctrl.levels = np.nanpercentile(volume, [0.5, 99.5]) if levels is None else levels
        for _, sl in self.ctrl.slices.items():
            sl.ctrl.image_layers['image'].image_item.setLookupTable(cmap.getLookupTable(alpha=True))
            sl.ctrl.image_layers['image'].pg_kwargs = {'mode': 'clip', 'levels': self.ctrl.levels}
        self._refresh()


class SliceView(QtWidgets.QWidget):
    """
    Window containing a volume slice
    """

    def __init__(self, topview: TopView, waxis, haxis, daxis, **kwargs):
        super(SliceView, self).__init__()
        self.topview = topview
        self.ctrl = SliceController(self, waxis, haxis, daxis, **kwargs)
        uic.loadUi(Path(__file__).parent.joinpath('sliceview.ui'), self)
        self.add_image_layer(slice_kwargs={'mode': 'clip'},
                             pg_kwargs={'opacity': 0.8}, name='image')
        self.add_image_layer(slice_kwargs={'volume': 'annotation', 'mode': 'clip'},
                             pg_kwargs={'opacity': 0.2}, name='annotation')
        self.add_image_layer(slice_kwargs={'volume': 'boundary', 'mode': 'clip'},
                             pg_kwargs={'opacity': 1}, name='boundary')
        # init the image display
        self.plotItem_slice.setAspectLocked(True)
        line_kwargs = {'movable': True, 'pen': pg.mkPen((0, 255, 0), width=3)}
        self.horizontal_line = pg.InfiniteLine(angle=0, pos=0, **line_kwargs)
        self.vertical_line = pg.InfiniteLine(angle=90, pos=0, **line_kwargs)
        self.plotItem_slice.addItem(self.horizontal_line)
        self.plotItem_slice.addItem(self.vertical_line)
        # connect signals and slots
        s = self.plotItem_slice.getViewBox().scene()
        self.proxy = pg.SignalProxy(s.sigMouseMoved, rateLimit=60, slot=self.mouseMoveEvent)
        s.sigMouseClicked.connect(self.mouseClick)

    def add_scatter(self):
        self.scatterItem = pg.ScatterPlotItem()
        self.plotItem_slice.addItem(self.scatterItem)

    def add_image_layer(self, name=None, **kwargs):
        """
        :param pg_kwargs: pyqtgraph setImage arguments: {'levels': None, 'lut': None,
        'opacity': 1.0}
        :param slice_kwargs: iblatlas.atlas.slice arguments: {'volume': 'image', 'mode': 'clip'}
        :return:
        """
        assert name is not None
        il = ImageLayer(**kwargs)
        self.ctrl.image_layers[name] = il
        self.plotItem_slice.addItem(il.image_item)

    def closeEvent(self, event):
        self.hide()

    def keyPressEvent(self, e):
        pass

    def mouseClick(self, event):
        if not event.double():
            return

    def mouseMoveEvent(self, scenepos):
        if isinstance(scenepos, tuple):
            scenepos = scenepos[0]
        else:
            return
        qpoint = self.ctrl.image_layers['image'].image_item.mapFromScene(scenepos)
        iw, ih, w, h, v, region = self.ctrl.cursor2xyamp(qpoint)
        self.label_v.setText(f"{v:.2f}")
        self.label_x.setText(f"{w * 1e6:.0f}")
        self.label_y.setText(f"{h * 1e6:.0f}")
        self.label_ix.setText(f"{iw:.0f}")
        self.label_iy.setText(f"{ih:.0f}")
        if region is None:
            self.label_region.setText("")
            self.label_acronym.setText("")
        else:
            self.label_region.setText(region['name'][0])
            self.label_acronym.setText(region['acronym'][0])

    def replace_image_layer(self, index, **kwargs):
        if index and len(self.imageItem) >= index:
            il = self.image_layers.pop(index)
            self.plotItem_slice.removeItem(il.image_item)
        self.add_image_layer(**kwargs)


class PgImageController:
    """
    Abstract class that implements mapping fr`om axes to voxels for any window.
    Not instantiated directly.
    """
    def __init__(self, win, res=25):
        self.qwidget = win
        self.transform = None  # affine transform image indices 2 data domain
        self.image_layers: dict = {}

    def cursor2xyamp(self, qpoint):
        """Used for the mouse hover function over image display"""
        iw, ih = self.cursor2ind(qpoint)
        v = self.im[iw, ih]
        w, h, _ = np.matmul(self.transform, np.array([iw, ih, 1]))
        return iw, ih, w, h, v

    def cursor2ind(self, qpoint):
        """ image coordinates over the image display"""
        iw = np.max((0, np.min((int(np.floor(qpoint.x())), self.nw - 1))))
        ih = np.max((0, np.min((int(np.round(qpoint.y())), self.nh - 1))))
        return iw, ih

    @property
    def imageItem(self):
        """returns the first image item"""
        return next((self.image_layers[k].image_item for k in self.image_layers))

    def set_image(self, pg_image_item, im, dw, dh, w0, h0, **pg_kwargs):
        """
        :param im:
        :param dw:
        :param dh:
        :param w0:
        :param h0:
        :param pgkwargs: og.ImageItem.setImage() parameters: level=None, lut=None, opacity=1
        :return:
        """
        self.im = im
        self.nw, self.nh = self.im.shape[0:2]
        pg_image_item.setImage(self.im, **pg_kwargs)
        transform = [dw, 0., 0., 0., dh, 0., w0, h0, 1.]
        self.transform = np.array(transform).reshape((3, 3)).T
        pg_image_item.setTransform(QTransform(*transform))

    def set_points(self, x=None, y=None):
        # at the moment brush and size are fixed! These need to be arguments
        # For the colour need to convert the colour to QtGui.QColor
        self.qwidget.scatterItem.setData(x=x, y=y, brush='b', size=5)


class ControllerTopView(PgImageController):
    """
    TopView ControllerTopView
    """
    def __init__(self, qmain: TopView, res: int = 25, volume='image', atlas=None, **kwargs):
        super(ControllerTopView, self).__init__(qmain)
        line_kwargs = {'movable': True, 'pen': pg.mkPen((0, 255, 0), width=3)}
        self.line_coronal = pg.InfiniteLine(angle=0, pos=0, **line_kwargs)
        self.line_sagittal = pg.InfiniteLine(angle=90, pos=0, **line_kwargs)
        self.settings = QtCore.QSettings('IBL', 'Atlas')
        self.highlight_region = None
        self.volume = volume
        self.atlas = AllenAtlas(res) if atlas is None else atlas
        self.atlas.regions.compute_hierarchy()
        self.fig_top = self.qwidget = qmain
        # Setup Coronal slice: width: ml, height: dv, depth: ap
        self.fig_coronal = SliceView(qmain, waxis=0, haxis=2, daxis=1)
        self.fig_coronal.setWindowTitle('Coronal Slice')
        self.set_slice(self.fig_coronal)
        self.fig_coronal.show()
        # Setup Sagittal slice: width: ap, height: dv, depth: ml
        self.fig_sagittal = SliceView(qmain, waxis=1, haxis=2, daxis=0)
        self.fig_sagittal.setWindowTitle('Sagittal Slice')
        self.set_slice(self.fig_sagittal)
        self.fig_sagittal.show()
        # Setup Horizontal slice: width: ml, height: ap, depth: dv
        self.fig_horizontal = SliceView(qmain, waxis=1, haxis=0, daxis=2)
        self.fig_horizontal.setWindowTitle('Horizontal Slice')
        self.set_slice(self.fig_horizontal)
        self.fig_horizontal.show()
        # The last figure is the brain tree architecture
        self.fig_brain_tree = BrainTree()
        for k, fig in self.figures.items():
            if (geom := self.settings.value(f'geometry_{k}', None)) is not None:
                fig.restoreGeometry(geom)

    @property
    def lines_sagittal(self):
        return [self.fig_coronal.vertical_line, self.fig_horizontal.horizontal_line, self.line_sagittal]

    @property
    def lines_coronal(self):
        return [self.fig_horizontal.vertical_line, self.fig_sagittal.vertical_line, self.line_coronal]

    @property
    def lines_horizontal(self):
        return [self.fig_coronal.horizontal_line, self.fig_sagittal.horizontal_line]

    @property
    def slices(self) -> dict:
        return {
            'coronal': self.fig_coronal,
            'sagittal': self.fig_sagittal,
            'horizontal': self.fig_horizontal,
        }

    @property
    def figures(self) -> dict:
        return self.slices | {'brain_tree': self.fig_brain_tree}

    def set_slice(self, fig, coord=0, mapping="Allen"):
        waxis, haxis, daxis = (fig.ctrl.waxis, fig.ctrl.haxis, fig.ctrl.daxis)
        # construct the transform matrix image 2 ibl coordinates
        dw = self.atlas.bc.dxyz[waxis]
        dh = self.atlas.bc.dxyz[haxis]
        wl = self.atlas.bc.lim(waxis)[0] - dw / 2
        hl = self.atlas.bc.lim(haxis)[0] - dh / 2
        # the ImageLayer object carries slice kwargs and pyqtgraph ImageSet kwargs
        # reversed order so the self.im is set with the base layer
        for layer_name, layer in fig.ctrl.image_layers.items():
            if layer_name == 'boundary':
                if self.highlight_region is None:
                    layer.image_item.setOpacity(0)
                    continue
                else:
                    ir = self.highlight_region
                    _, iir = self.atlas.regions.descendants(self.atlas.regions.id[ir], return_indices=True)
                    slice_labels = self.atlas.slice(coord, axis=daxis, mapping='Allen', volume='rindex', mode='clip')
                    _slice, _ = ismember(slice_labels, iir)
                    _slice = self.atlas.compute_boundaries(_slice)
                    _slice = np.tile(_slice.astype(np.uint8)[:, :, np.newaxis], (1, 1, 4)) * 255
                    layer.image_item.setOpacity(1)
            elif layer_name == 'image':
                _slice = self.atlas.slice(coord, axis=daxis, mapping=mapping, volume=self.volume, **layer.slice_kwargs)
            else:
                _slice = self.atlas.slice(coord, axis=daxis, mapping=mapping, **layer.slice_kwargs)
            fig.ctrl.set_image(layer.image_item, _slice, dw, dh, wl, hl, **layer.pg_kwargs)

        fig.ctrl.slice_coord = coord

    def set_top(self):
        self.atlas.compute_surface()
        img = self.atlas.top.T.copy()
        img[np.isnan(img)] = np.nanmin(img)  # img has dims ml, ap
        if (ir := self.highlight_region) is not None:
            _, iir = self.atlas.regions.descendants(self.atlas.regions.id[ir], return_indices=True)
            bounds = np.any(ismember(self.atlas.label, iir)[0], axis=self.atlas.xyz2dims[-1]).astype(bool).T
            img[bounds] = np.nan
        if np.diff(ismember(self.atlas.dims2xyz, [0, 1])[1])[0] > 0:
            img = img.T
        dw, dh = (self.atlas.bc.dxyz[0], self.atlas.bc.dxyz[1])
        wl, hl = (self.atlas.bc.xlim, self.atlas.bc.ylim)
        self.set_image(self.image_layers['top'].image_item, img, dw, dh, wl[0], hl[0])

    def set_scatter(self, fig, coord=0):
        waxis = fig.ctrl.waxis
        # dealing with coronal slice
        if waxis == 0:
            idx = np.where(self.scatter_data_ind[:, 1] == self.atlas.bc.y2i(coord))[0]
            x = self.scatter_data[idx, 0]
            y = self.scatter_data[idx, 2]
        else:
            idx = np.where(self.scatter_data_ind[:, 0] == self.atlas.bc.x2i(coord))[0]
            x = self.scatter_data[idx, 1]
            y = self.scatter_data[idx, 2]
        fig.ctrl.set_points(x, y)


class SliceController(PgImageController):

    def __init__(self, fig, waxis=None, haxis=None, daxis=None, wdir=1, hdir=1):
        """
        :param waxis: brain atlas axis corresponding to display abscissa (coronal: 0, sagittal: 1)
        :param haxis: brain atlas axis corresponding to display ordinate (coronal: 2, sagittal: 2)
        :param daxis: brain atlas axis corresponding to display abscissa (coronal: 1, sagittal: 0)
        """
        super(SliceController, self).__init__(fig)
        self.waxis = waxis
        self.haxis = haxis
        self.daxis = daxis
        self.wdir = wdir
        self.hdir = hdir

    def cursor2xyamp(self, qpoint):
        """
        Extends the superclass method to also get the brain region from the model
        :param qpoint:
        :return:
        """
        iw, ih, w, h, v = super(SliceController, self).cursor2xyamp(qpoint)
        ctrl = self.qwidget.topview.ctrl
        xyz = np.zeros(3)
        xyz[np.array([self.waxis, self.haxis, self.daxis])] = [w, h, self.slice_coord]
        mapping = self.qwidget.topview.comboBox_mappings.currentText()
        try:
            region = ctrl.atlas.regions.get(ctrl.atlas.get_labels(xyz, mapping=mapping))
        except ValueError:
            region = None
        i = ctrl.atlas._lookup(xyz, mode='clip')
        vol = ctrl.atlas.image if isinstance(ctrl.volume, str) else ctrl.volume
        v = np.take(vol, i)
        return iw, ih, w, h, v, region


@dataclass
class ImageLayer:
    """
    Class for keeping track of image layers.
    :param image_item
    :param pg_kwargs: pyqtgraph setImage arguments: {'levels': None, 'lut': None, 'opacity': 1.0}
    :param slice_kwargs: iblatlas.atlas.slice arguments: {'volume': 'image', 'mode': 'clip'}
    :param
    """
    image_item: pg.ImageItem = field(default_factory=pg.ImageItem)
    pg_kwargs: dict = field(default_factory=lambda: {})
    slice_kwargs: dict = field(default_factory=lambda: {'volume': 'image', 'mode': 'clip'})


def view(res=25, title=None, atlas=None):
    """ application entry point """
    qt.create_app()
    av = TopView._get_or_create(title=title, res=res, atlas=atlas)
    av.show()
    return av


def main():
    app = QtWidgets.QApplication([])
    w = TopView()
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
