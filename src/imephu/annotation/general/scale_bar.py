import dataclasses
import math
from copy import deepcopy
from typing import Any, Tuple

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import WCS
from matplotlib.patches import PathPatch
from matplotlib.path import Path

from imephu.annotation import Annotation
from imephu.geometry import pixel_scales


@dataclasses.dataclass()
class _ScaleBarParameters:
    pixels: float
    angle: float
    units: str


class ScaleBarLineAnnotation(Annotation):
    """An annotation for drawing the line of a scale bar.

    The position of the scale bar's centre must be specified in axes coordinates. Axes
    coordinates range from (0, 0) (bottom left of the axes) to (1, 1) (top right of the
    axes).

    Parameters
    ----------
    left_edge: pair of float
        The position of the left edge of the scale bar's horizontal line, in pixels.
    minimum_length: float
        The minimum length of the scale bar's horizontal line, in pixels.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    color: color, default: "black"
        The color of the scale bar.
    linewidth: float, default: 0.5
        The line width.
    **kwargs: dict, optional
        Additional keyword arguments, which must be valid arguments for both paths and
        text.
    """

    def __init__(
        self,
        left_edge: Tuple[float, float],
        minimum_length: float,
        wcs: WCS,
        color: str = "black",
        linewidth: float = 0.5,
        **kwargs: Any,
    ):
        self._left_edge = left_edge
        self._minimum_length = minimum_length
        self._wcs = wcs
        self._kwargs = deepcopy(kwargs)
        self._kwargs["edgecolor"] = color
        self._kwargs["facecolor"] = color
        self._kwargs["linewidth"] = linewidth

    def add_to(self, ax: WCSAxes) -> None:
        """Add the scale bar to a finder chart.

        ax: `~astropy.visualization.wcsaxes.WCSAxes`
            WCS axes object.
        """
        # Add the horizontal line
        pixel_scale = pixel_scales(self._wcs)[0]
        parameters = self._parameters(self._minimum_length, pixel_scale)
        x_min = ax.get_xlim()[0] + self._left_edge[0]
        x_max = x_min + parameters.pixels
        y_mid = ax.get_ylim()[0] + self._left_edge[1]
        horizontal_path = Path([(x_min, y_mid), (x_max, y_mid)], closed=False)  # noqa
        horizontal_path_patch = PathPatch(horizontal_path, **self._kwargs)
        horizontal_path_patch.set_clip_on(False)
        ax.add_patch(horizontal_path_patch)

        # Add the vertical lines
        y_min = y_mid - 10
        y_max = y_mid + 10
        left_vertical_path = Path(
            [(x_min, y_min), (x_min, y_max)], closed=False  # noqa
        )  # noqa
        left_vertical_path_patch = PathPatch(left_vertical_path, **self._kwargs)
        left_vertical_path_patch.set_clip_on(False)
        ax.add_patch(left_vertical_path_patch)
        right_vertical_path = Path(
            [(x_max, y_min), (x_max, y_max)], closed=False  # noqa
        )  # noqa
        right_vertical_path_patch = PathPatch(right_vertical_path, **self._kwargs)
        right_vertical_path_patch.set_clip_on(False)
        ax.add_patch(right_vertical_path_patch)

        # Add the text
        text_kwargs = self._kwargs.copy()
        del text_kwargs["edgecolor"]
        del text_kwargs["facecolor"]
        del text_kwargs["linewidth"]
        text_kwargs["color"] = self._kwargs["facecolor"]
        x_mid = (x_min + x_max) / 2
        text = f"{parameters.angle:.1f} {parameters.units}"
        ax.text(
            x_mid,
            y_mid - 20,
            text,
            fontsize=10,
            horizontalalignment="center",
            **text_kwargs,
        )

    def rotate(self, pivot: SkyCoord, angle: Angle) -> "Annotation":
        """Rotate this annotation around a pivot and return the result.

        The rotation angle is an angle on the sky, measured from north to east.

        Parameters
        ----------
        pivot: `~astropy.coordinates.SkyCoord`
            Point around which to rotate the annotation.
        angle: `~astropy.coordinates.Angle`
            Angle of rotation, measured from north towards the east.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the rotation.
        """
        raise NotImplementedError()

    def translate(self, displacement: Angle) -> "Annotation":
        """
        Move this annotation along a displacement vector and return the result.

        Not implemented.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the translation.
        """
        raise NotImplementedError()

    @staticmethod
    def _parameters(minimum_length: float, pixel_scale: Angle) -> _ScaleBarParameters:
        minimum_angle = minimum_length * pixel_scale
        if minimum_angle <= 50 * u.arcsec:
            minimum_angle_arcsec = minimum_angle.to_value(u.arcsec)
            angle_value = ScaleBarLineAnnotation._preferred_length_value(
                minimum_angle_arcsec
            )
            angle = angle_value * u.arcsec
            units = "arcsec"
        else:
            minimum_angle_arcmin = minimum_angle.to_value(u.arcmin)
            angle_value = ScaleBarLineAnnotation._preferred_length_value(
                minimum_angle_arcmin
            )
            angle = angle_value * u.arcmin
            units = "arcmin"
        pixels = float(angle / pixel_scale)
        return _ScaleBarParameters(pixels=pixels, angle=angle_value, units=units)

    @staticmethod
    def _preferred_length_value(minimum_length: float) -> float:
        """
        Return the preferred length given a minimum length value.

        The preferred value p is chosen as follows for a minimum length l_min.

        1. If the l_min is less than or equal to 1, p is equal to 1.
        2. Otherwise, if l_min is less than or equal to 1, p is l_min rounded up to the
           nearest integer.
        3. Otherwise, p is l_min rounded up to the next integer divisible by 10.

        Parameters
        ----------
        minimum_length: float
            Minimum scale bar length.

        Returns
        -------
        float
            The preferred scale bar length.

        """
        if minimum_length <= 1:
            return 1
        if minimum_length <= 10:
            return math.ceil(minimum_length)

        normalizing_factor = math.pow(10, math.floor(math.log10(minimum_length)))
        minimum_length /= normalizing_factor
        minimum_length = math.ceil(minimum_length)
        minimum_length *= normalizing_factor
        return minimum_length
