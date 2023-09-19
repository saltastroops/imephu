from __future__ import annotations

import bisect
import os
import warnings
from datetime import datetime
from io import BytesIO
from typing import (
    Any,
    BinaryIO,
    Callable,
    Dict,
    Generator,
    List,
    Optional,
    Tuple,
    Union,
)

import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.visualization.interval import AsymmetricPercentileInterval
from astropy.visualization.mpl_normalize import simple_norm
from astropy.visualization.wcsaxes.core import WCSAxesSubplot
from astropy.wcs import WCS, FITSFixedWarning
from matplotlib.figure import Figure

from imephu.annotation import Annotation
from imephu.geometry import pixel_scales
from imephu.service.survey import load_fits
from imephu.utils import Ephemeris

"""A finder chart."""


class FinderChart:
    """A finder chart for an astronomical observation.

    The finder chart is generated from a FITS file. By default the finder chart just
    shows the region with an overlaid coordinate grid. The axes show WCS coordinates.
    Use the `add_annotation` method to add more content to the finder chart.

    By default, the size of the finder chart is the same as that of the FITS file. You
    may set the `max_size` property to choose a non-default maximum size.

    You can display the finder chart on the screen or save it as a file.

    This class uses Matplotlib for generating the finder chart.

    .. warning:: Be careful when working with several finder charts as all the finder
       charts use the same Matplotlib figure.

    Parameters
    ----------
    name: `str`, `path-like` or `binary file-like`
        FITS file to display.
    max_size: `~astropy.coordinates.Angle`
        Maximum size of the generated finder chart.
    """

    def __init__(
        self,
        name: Union[str, BinaryIO, os.PathLike[Any]],
        max_size: Optional[Angle] = None,
    ):
        self._hdu = fits.open(name)[0]
        # The FITS data is read in only when it is needed. To avoid trying to read from
        # a closed stream later on, we thus force the data to be read in immediately.
        self._data = self._hdu.data
        self._metadata: Dict[str, Any] = dict()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FITSFixedWarning)
            self._wcs = WCS(self._hdu)
        self._annotations: List[Annotation] = []
        self._max_size = max_size

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

    @staticmethod
    def for_time_interval(
        start: datetime,
        end: datetime,
        ephemerides: List[Ephemeris],
        max_track_length: Angle,
        create_finder_chart: Callable[[List[Ephemeris]], "FinderChart"],
    ) -> Generator[Tuple["FinderChart", Tuple[datetime, datetime]], None, None]:
        """Create finder charts for a time interval.

        This method is intended for non-sidereal targets, which may require multiple
        finder charts to cover their track over a time interval. It creates finder
        charts according to the following rules:

        * Taken together, the finder charts cover the whole time interval.
        * The track length on each finder chart does not exceed `max_track_length`.
        * Each finder chart is created using ``create_finder_chart``.

        The method is completely agnostic as to what a finder chart should look like;
        this decision is left completely to ``create_finder_chart``. In particular, if
        you need any annotations for the non-sidereal nature, it is up to
        ``create_finder_chart`` to provide them.

        The ``create_finder_chart`` function has to accept as its single argument the
        list of ephemerides to include on the finder chart.

        The finder charts are returned along with the time interval they cover. The time
        interval is also added as a tuple of `~datetime.datetime` values with the key
        ``valid_for`` to the finder chart's metadata.

        Parameters
        ----------
        start: `~datetime.datetime`
            The start time of the interval to be covered by the finder charts. This must
            be a timezone-aware datetime.
        end: `~datetime.datetime`
            The end time of the interval to be covered by the finder charts. This must
            be a timezone-aware datetime.
        ephemerides: list of `~imephu.utils.Ephemeris`
            The list of ephemerides. The time interval from ``start`` to ``end`` must be
            fully covered by the ephemerides.
        max_track_length: float
            The maximum length a track may have on a finder chart, as an angle on the
            sky.
        create_finder_chart: function
            The function for creating a finder chart from a list of ephemerides.

        Yields
        ------
        tuple of a `~imephu.finder_chart.FinderChart` and a datetime interval
            The generated finder charts along with the time intervals they cover.
        """
        # The start and end time must be timezone-aware
        if start.tzinfo is None or start.tzinfo.utcoffset(None) is None:
            raise ValueError("The start time must be timezone-aware.")
        if end.tzinfo is None or end.tzinfo.utcoffset(None) is None:
            raise ValueError("The end time must be timezone-aware.")

        if start >= end:
            raise ValueError("The start time must be earlier than the end time.")

        if start < ephemerides[0].epoch:
            raise ValueError("The start time must not be earlier than the first epoch.")
        if end > ephemerides[-1].epoch:
            raise ValueError("The end time must not be later than the last epoch.")

        if max_track_length.to_value(u.arcmin) <= 0:
            raise ValueError("The maximum track length must be positive.")

        # Find the smallest interval covering the time interval
        all_times = [e.epoch for e in ephemerides]
        start_index = bisect.bisect_right(all_times, start) - 1
        end_index = bisect.bisect(all_times, end)
        if end_index > 0 and all_times[end_index - 1] == end:
            end_index -= 1

        # This should never happen, but let's rule out intervals with zero length
        if start_index == end_index:
            raise ValueError("The interval must have a positive length")

        # Split the ephemerides so that the maximum track length isn't exceeded for each
        # group. We assume a linear path, so that the track length is equal to the angle
        # between the first and last position.
        current_group = [ephemerides[start_index]]
        groups = [current_group]
        for i in range(start_index + 1, end_index + 1):
            next_ephemeris = ephemerides[i]
            # Subsequent positions must not be more than max_track_length apart
            if (
                current_group[-1].position.separation(next_ephemeris.position)
                > max_track_length
            ):
                raise ValueError(
                    "The maximum track length on a finder chart is exceeded."
                )
            track_length = current_group[0].position.separation(next_ephemeris.position)
            if track_length <= max_track_length:
                current_group.append(next_ephemeris)
            else:
                current_group = [current_group[-1], next_ephemeris]
                groups.append(current_group)

        # Create the finder charts
        for group in groups:
            valid_for = (group[0].epoch, group[-1].epoch)
            finder_chart = create_finder_chart(group)
            finder_chart.add_metadata("valid_for", valid_for)
            yield finder_chart, valid_for

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

    @property
    def metadata(self) -> Dict[str, Any]:
        """
        Return the metadata for the finder chart.

        Metadata can  be added with the `add_metadata` method. A shallow copy of the
        metadata is returned.

        Returns
        -------
        dict
            The metadata of the finder chart.
        """
        return self._metadata.copy()

    def add_metadata(self, key: str, value: Any) -> None:
        """
        Add a key-value to the finder chart's metadata.

        If the key existrs in the metadata already, the existing value for the key is
        replaced.

        Parameters
        ----------
        key: `str`
            Key.
        value: `~typing.Any`
            Value.
        """
        self._metadata[key] = value

    @property
    def max_size(self) -> Optional[Angle]:
        return self._max_size

    @max_size.setter
    def max_size(self, new_max_size: Optional[Angle]) -> None:
        self._max_size = new_max_size

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
            if hasattr(name, "write"):
                name.write(pdf.getvalue())
            else:
                with open(name, "wb") as f:
                    f.write(pdf.getvalue())
        else:
            plt.savefig(name, format=format, bbox_inches="tight")
        plt.close(figure)

    def _create_plot(self) -> Figure:
        figure = plt.figure(figsize=(10, 9))
        ax = plt.subplot(projection=self._wcs)

        self._add_fits_content(ax)

        # Enforce the maximum finder chart size. This needs to be done after the FITS
        # file has been added, as that sets the plot limits according to the image size.
        self._enforce_max_size(ax)

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

    def _enforce_max_size(self, ax: WCSAxesSubplot) -> None:
        if self.max_size is None:
            return

        # Get the current plot limits
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        # Get the angular size in x and y direction, and assume the finder chart size
        # to be the greater of the two
        x_scale, y_scale = pixel_scales(ax.wcs)
        x_size_pixels = x_max - x_min
        y_size_pixels = y_max - y_min
        x_size = x_size_pixels * x_scale
        y_size = y_size_pixels * y_scale
        size = max(x_size, y_size)

        # If we are within 5% of the maximum size, there is nothing to do
        if size < 1.05 * self.max_size:
            return

        # Update the limits in x direction
        shrink_factor = float(self.max_size / size)
        x_center = (x_min + x_max) / 2
        x_size_data_coords_new = shrink_factor * x_size_pixels
        x_min_new = x_center - x_size_data_coords_new / 2
        x_max_new = x_center + x_size_data_coords_new / 2
        ax.set_xlim(left=x_min_new, right=x_max_new)

        # Update the limits in y direction
        y_center = (y_min + y_max) / 2
        y_size_data_coords_new = shrink_factor * y_size_pixels
        y_min_new = y_center - y_size_data_coords_new / 2
        y_max_new = y_center + y_size_data_coords_new / 2
        ax.set_ylim(bottom=y_min_new, top=y_max_new)

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
