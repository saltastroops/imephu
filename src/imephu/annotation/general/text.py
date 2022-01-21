from copy import deepcopy
from typing import Any, Literal

from astropy.coordinates import Angle, SkyCoord
from astropy.units import Quantity
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import WCS

from imephu.annotation import Annotation
from imephu.geometry import rotate, sky_position_to_pixel, translate


class TextAnnotation(Annotation):
    """An annotation for adding text to a plot.

    .. note::
        The `rotate` method rotates the text position, but not the text itself. If you
        want to rotate the text rather than the text position, you should use the
        ``rotation`` parameter with the angle in degrees.


    Parameters
    ----------
    position: `~astropy.coordinates.SkyCoord`
        The right ascension and declination of the point where to put the text.
    s: str
        The text.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    color: str, default: "black"
        The text color.
    horizontalalignment: {"center", "left", "right"}, default: "center"
        The horizontal alignment of the text.
    verticalalignment: {"baseline", "bottom", "center", "center_baseline", "top"}, default: "center"
        The vertical alignment of the text.
    **kwargs: dict, optional
        Additional keyword arguments, which will be passed to Matplotlib's
        `~matplotlib.axes.Axes.text` method for adding text to plot Axes.
    """  # noqa: E501

    def __init__(
        self,
        position: SkyCoord,
        s: str,
        wcs: WCS,
        color: str = "black",
        horizontalalignment: Literal["center", "left", "right"] = "center",
        verticalalignment: Literal[
            "baseline", "bottom", "center", "center_baseline", "top"
        ] = "center",
        **kwargs: Any,
    ):
        self._position = position
        self._s = s
        self._wcs = wcs
        self._kwargs = deepcopy(kwargs)
        self._kwargs["color"] = color
        self._kwargs["horizontalalignment"] = horizontalalignment
        self._kwargs["verticalalignment"] = verticalalignment

    def add_to(self, ax: WCSAxes) -> None:
        """Add the text to a finder chart.

        ax: WCSAxes
            WCS axes object.
        """
        position_px = sky_position_to_pixel(self._position, self._wcs)
        ax.text(position_px[0], position_px[1], self._s, **self._kwargs)  # noqa

    def rotate(self, center: SkyCoord, angle: Angle) -> Annotation:
        """Rotate this annotation around a center and return the result.

        The rotation angle is an angle on the sky, measured from north to east.

        Parameters
        ----------
        center: `~astropy.coordinates.SkyCoord`
            Point around which to rotate the annotation.
        angle: `~astropy.units.Quantity` ["angle"]
            Angle of rotation, measured from north to east.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the rotation.
        """
        rotated_annotation = deepcopy(self)
        rotated_position = rotate(self._position, center, angle, self._wcs)
        rotated_annotation._position = rotated_position
        return rotated_annotation

    def translate(self, displacement: Quantity) -> Annotation:
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the translation.
        """
        translated_annotation = deepcopy(self)
        translated_position = translate(self._position, displacement)
        translated_annotation._position = translated_position
        return translated_annotation
