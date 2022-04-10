from copy import deepcopy
from typing import Any, Literal, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt
from astropy.coordinates import Angle, SkyCoord
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import WCS

from imephu.annotation import Annotation
from imephu.geometry import rotate, sky_position_to_pixel, translate


class TextAnnotation(Annotation):
    """An annotation for adding text to a plot.

    The text position can be specified as a position on the sky (as an instance of
    `~astropy.coordinates.SkyCoord`) or as a point in axes coordinates (as an array of
    `float` values. Axes coordinates range from (0, 0) (bottom left of the axes) to
    (1, 1) (top right of the axes).

    You can only rotate and translate the annotation if the text position is specified
    as a position on the sky.

    .. note::
        The `rotate` method rotates the text position, but not the text itself. If you
        want to rotate the text rather than the text position, you should use the
        ``rotation`` parameter with the angle in degrees.


    Parameters
    ----------
    position: `~astropy.coordinates.SkyCoord` or
        The right ascension and declination of the point where to put the text.
    s: str
        The text.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    color: color, default: "black"
        The text color.
    horizontalalignment: {"center", "left", "right"}, default: "center"
        The horizontal alignment of the text.
    verticalalignment: {"baseline", "bottom", "center", "center_baseline", "top"}, default: "center"
        The vertical alignment of the text.
    clip_on: bool, default: True
        Whether to clip the text at the plot boundaries.
    **kwargs: dict, optional
        Additional keyword arguments, which will be passed to Matplotlib's
        `~matplotlib.axes.Axes.text` method for adding text to plot Axes.
    """  # noqa: E501

    def __init__(
        self,
        position: Union[SkyCoord, Sequence[float], npt.NDArray[np.float_]],
        s: str,
        wcs: WCS,
        color: Union[str, Tuple[float, float, float]] = "black",
        horizontalalignment: Literal["center", "left", "right"] = "center",
        verticalalignment: Literal[
            "baseline", "bottom", "center", "center_baseline", "top"
        ] = "center",
        clip_on: bool = True,
        **kwargs: Any,
    ):
        self._position = position
        self._s = s
        self._wcs = wcs
        self._kwargs = deepcopy(kwargs)
        self._kwargs["color"] = color
        self._kwargs["horizontalalignment"] = horizontalalignment
        self._kwargs["verticalalignment"] = verticalalignment
        self._kwargs["clip_on"] = clip_on

    def add_to(self, ax: WCSAxes) -> None:
        """Add the text to a finder chart.

        ax:  `~astropy.visualization.wcsaxes.WCSAxes`
            WCS axes object.
        """
        if type(self._position) == SkyCoord:
            position_px = sky_position_to_pixel(self._position, self._wcs)
            ax.text(position_px[0], position_px[1], self._s, **self._kwargs)  # noqa
        else:
            ax.text(
                self._position[0],
                self._position[1],
                self._s,
                transform=ax.transAxes,
                **self._kwargs,
            )

    def rotate(self, pivot: SkyCoord, angle: Angle) -> "TextAnnotation":
        """Rotate this annotation around a pivot and return the result.

        The rotation angle is an angle on the sky, measured from north to east.

        Parameters
        ----------
        pivot: `~astropy.coordinates.SkyCoord`
            Point around which to rotate the annotation.
        angle: `~astropy.coordinates.Angle`
            Angle of rotation, measured from north to east.

        Returns
        -------
        `~imephu.annotation.general.TextAnnotation`
            The annotation resulting from the rotation.
        """
        if type(self._position) is not SkyCoord:
            raise NotImplementedError(
                "The text can only be rotated if its position is an instance of "
                "astropy.coord.SkyCoord."
            )
        rotated_annotation = deepcopy(self)
        rotated_position = rotate(self._position, pivot, angle, self._wcs)
        rotated_annotation._position = rotated_position
        return rotated_annotation

    def translate(self, displacement: Angle) -> "TextAnnotation":
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.general.TextAnnotation`
            The annotation resulting from the translation.
        """
        if type(self._position) is not SkyCoord:
            raise NotImplementedError(
                "The text can only be translated if its position is an instance of "
                "astropy.coord.SkyCoord."
            )
        translated_annotation = deepcopy(self)
        translated_position = translate(self._position, displacement)
        translated_annotation._position = translated_position
        return translated_annotation
