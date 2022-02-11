from __future__ import annotations

import os
from io import BytesIO
from typing import Any, BinaryIO, List, Optional, Union

import matplotlib.pyplot as plt
import pikepdf
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.visualization.interval import AsymmetricPercentileInterval
from astropy.visualization.mpl_normalize import simple_norm
from astropy.visualization.wcsaxes.core import WCSAxesSubplot
from astropy.wcs import WCS
from matplotlib.figure import Figure

import imephu
from imephu.annotation import Annotation
from imephu.service.survey import load_fits

"""A finder chart."""


class FinderChart:
    """A finder chart for an astronomical observation.

    The finder chart is generated from a FITS file. By default the finder chart just
    shows the region with an overlaid coordinate grid. The axes show WCS coordinates.
    Use the `add_annotation` method to add more content to the finder chart.

    You can display the finder chart on the screen or save it as a file.

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
        self._annotations: List[Annotation] = []

    @staticmethod
    def from_survey(survey: str, fits_center: SkyCoord, size: Angle) -> "FinderChart":
        """Create a finder chart from a sky survey FITS image.

        Parameters
        ----------
        survey: `str`
            The name of the survey to query for the FITS file.
        fits_center: `~astropy.coordinates.SkyCoord`
            The center position of the loaded FITS file.
        size: `~astropy.coordinates.Angle`
            The width and height of the FITS file image, as an angle on the sky.
        """
        fits_image = load_fits(survey, fits_center, size)
        return FinderChart(fits_image)

    @property
    def wcs(self) -> WCS:
        """Return a deep copy of the WCS object of the finder chart.

        Returns
        -------
        `~astropy.wcs.WCS`
            A deep copy of the WCS object.
        """
        return self._wcs.deepcopy()

    def add_annotation(self, annotation: Annotation) -> None:
        """Add an annotation to the finder chart.

        Annotations will be plotted onto the finder chart in the order they have been
        added. So, for example, the annotation added last will be output on top of all
        the other annotations.
        """
        self._annotations.append(annotation)

    def show(self) -> None:
        """Display the finder chart on the screen."""
        figure = self._create_plot()
        plt.show()
        plt.close(figure)

    def save(
        self,
        name: Union[str, BinaryIO, os.PathLike[str], os.PathLike[bytes]],
        format: Optional[str] = None,
    ) -> None:
        """Save the finder chart in a file.

        If `format` is not set, the file extension of `name` is used to figure out the
        file or, if there is no extension, Matplotlib's default is used. See
        Matplotlib's `matplotlib.pyplot.savefig` function for more details regarding
        the format.

        Parameters
        ----------
        name: `str`, `path-like` or `binary file-like`
            The file to which the finder chart is saved. An existing file is replaced.
        format: `str`, optional
            The format in which to store the finder chart.
        """
        figure = self._create_plot()
        if format and format.lower() == "pdf":
            pdf = BytesIO()
            plt.savefig(pdf, format=format, bbox_inches="tight")
            pdf_with_metadata = FinderChart._update_pdf_metadata(pdf.getvalue())
            if hasattr(name, "write"):
                name.write(pdf_with_metadata)  # type: ignore
            else:
                with open(name, "wb") as f:  # type: ignore
                    f.write(pdf_with_metadata)
        else:
            plt.savefig(name, format=format, bbox_inches="tight")
        plt.close(figure)

    def _create_plot(self) -> Figure:
        figure = plt.figure(figsize=(10, 9))
        ax = plt.subplot(projection=self._wcs)

        self._add_fits_content(ax)
        self._update_axes(ax)

        ax.grid(True, color="blue", alpha=0.2)

        for annotation in self._annotations:
            annotation.add_to(ax)

        return figure

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

    def _update_axes(self, ax: WCSAxesSubplot) -> None:
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

        x = ax.coords[0]
        y = ax.coords[1]
        x.set_axislabel(axis_type_names[x_axis_type])
        y.set_axislabel(axis_type_names[y_axis_type])

        x.display_minor_ticks(True)
        y.display_minor_ticks(True)

        x.set_ticks(size=7)
        y.set_ticks(size=7)

    @staticmethod
    def _update_pdf_metadata(pdf: bytes) -> bytes:
        """Update a pdf by setting some metadata values.

        The following metadata values are set:

        - ``dc:title``: The string ``imephu x.y.z`` is assigned (where ``x.y.z`` is the
          version).
        - ``xmp:CreatorTool``: A string of the form "imephu x.y.z" is assigned.

        The resulting pdf is returned; the original pdf remains unchanged.

        See the Core properties chapter of the `XMP specification
        <https://wwwimages2.adobe.com/content/dam/acom/en/devnet/xmp/pdfs/XMP%20SDK%20Release%20cc-2016-08/XMPSpecificationPart1.pdf>`_
        for more details about the metadata properties.

        Parameters
        ----------
        pdf: `bytes`
            The pdf whose ``xml:creatorTool`` metadata value should be set.

        Returns
        -------
        `bytes`
            The updated pdf content.
        """
        with pikepdf.open(BytesIO(pdf)) as document:
            with document.open_metadata() as meta:
                meta["dc:title"] = "Finder Chart"
                meta["xmp:CreatorTool"] = f"imephu {imephu.__version__}"
            updated_pdf = BytesIO()
            document.save(updated_pdf)
            return updated_pdf.getvalue()
