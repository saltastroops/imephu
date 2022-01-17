from __future__ import annotations
import os
from typing import Any, BinaryIO, Union, Optional

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.mpl_normalize import simple_norm
from astropy.visualization.interval import AsymmetricPercentileInterval
from astropy.visualization.wcsaxes.core import WCSAxesSubplot

"""A finder chart."""


class FinderChart:
    """A finder chart for an astronomical observation.

    The finder chart is generated from a FITS file. By default the finder chart just
    shows the region with an overlaid coordinate grid. The axes show WCS coordinates.
    Use the `add_annotation` method to add more content to the finder chart.

    You can display the finder chart on the screen or save it as a file. If you save it
    as a pdf, metadata is added to the file if it has been added with the `add_metadata`
    method.

    This class uses Matplotlib for generating the finder chart.

    .. warning:: Be careful when working with several finder charts as all the finder
       charts use the same Matplotlib figure.

    Parameters
    ----------
    name: `str`, `path-like` or `binary file-like`
        FITS file to display.
    """

    def __init__(self, name: Union[str, BinaryIO, os.PathLike[Any]]):
        self._hdu = fits.open(name)[0]
        # The FITS data is read in only when it is needed. To avoid trying to read from
        # a closed stream later on, we thus force the data to be read in immediately.
        self._data = self._hdu.data
        self._wcs = WCS(self._hdu)

    def show(self) -> None:
        """Display the finder chart on the screen."""
        self._create_plot()
        plt.show()

    def save(
        self,
        name: Union[str, BinaryIO, os.PathLike[Any]],
        format: Optional[str] = None,
    ) -> None:
        """Save the finder chart in a file.

        If `format` is not set, the file extension of `name` is used to figure out the
        file or, if there is no extension, Matplotlib's default is used. See
        Matplotlib's `matplotlib.pyplot.savefig` function for more details regarding
        the format.

        In case of a pdf the finder chart's metadata is stored as pdf metadata in the
        generated file.

        Parameters
        ----------
        name: `str`, `path-like` or `binary file-like`
            The file to which the finder chart is saved. An existing file is replaced.
        format: `str`, optional
            The format in which to store the finder chart.
        """
        self._create_plot()
        plt.savefig(name, format=format)

    def _create_plot(self) -> None:
        ax = plt.subplot(projection=self._wcs)

        self._add_fits_content(ax)
        self._update_axes(ax)

        # ax.coords.grid(True)

    def _add_fits_content(self, ax: WCSAxesSubplot) -> None:
        """
         Add content of the FITS file to the plot.

        The code has been adapted from APLpy.
        """
        cmap = plt.cm.get_cmap("gist_yarg")

        # value range in the image
        pmin = 0.25
        pmax = 99.75
        interval = AsymmetricPercentileInterval(pmin, pmax, n_samples=10000)
        vmin, vmax = interval.get_limits(self._hdu.data)
        vmin = -0.1 * (vmax - vmin) + vmin
        vmax = 0.1 * (vmax - vmin) + vmax

        normalizer = simple_norm(self._data, power=2, min_cut=vmin, max_cut=vmax)

        ax.imshow(
            self._hdu.data,
            cmap=cmap,
            interpolation="nearest",
            origin="lower",
            norm=normalizer,
            aspect="equal",
        )

    def _update_axes(self, ax: WCSAxesSubplot):
        axis_type_names = {
            "pos.eq.dec": "Dec (ICRS)",
            "pos.eq.ra": "RA (ICRS)",
        }
        axis_types = self._wcs.world_axis_physical_types
        x_axis_type = axis_types[0]
        y_axis_type = axis_types[1]
        if x_axis_type != "pos.eq.ra":
            raise ValueError(
                "Only pos.eq.ra is supported for the physical type of "
                "the first world axis"
            )
        if y_axis_type != "pos.eq.dec":
            raise ValueError(
                "Only pos.eq.dec is supported for the physical type of "
                "the second world axis"
            )

        ax.coords[0].set_axislabel(axis_type_names[x_axis_type])
        ax.coords[1].set_axislabel(axis_type_names[y_axis_type])
