from copy import deepcopy
from typing import Any, Tuple, Union, cast

from astropy.coordinates import Angle, SkyCoord
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import WCS
from matplotlib.patches import Circle

from imephu.annotation import Annotation
from imephu.geometry import pixel_scales, rotate, sky_position_to_pixel, translate


class CircleAnnotation(Annotation):
    """An annotation for adding a circle to a plot.

    The radius of the circle is given as an angular distance on the sky, in the
    direction of the right ascension. A circle is drawn irrespective of the pixel scales
    for right ascension and declination. Hence the circle in general is not a circle on
    the sky but just on the finder chart.

    Parameters
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The right ascension and declination of the circle's center.
    radius: `~astropy.coordinates.Angle`
        The radius of the circle, as an angular distance on the sky in right ascension
        direction.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    edgecolor: color, default: "black"
        The edge color.
    facecolor: color, default: "none"
        The filling color.
    **kwargs: dict, optional
        Additional keyword arguments, which will be passed to Matplotlib's
        `~matplotlib.patches.Circle` patch constructor.
    """

    def __init__(
        self,
        center: SkyCoord,
        radius: Angle,
        wcs: WCS,
        edgecolor: Union[str, Tuple[float, float, float]] = "black",
        facecolor: Union[str, Tuple[float, float, float]] = "none",
        **kwargs: Any,
    ):
        self._center = center
        self._radius = radius
        self._wcs = wcs
        self._kwargs = deepcopy(kwargs)
        self._kwargs["edgecolor"] = edgecolor
        self._kwargs["facecolor"] = facecolor

    def add_to(self, ax: WCSAxes) -> None:
        """Add the circle to a finder chart.

        ax: `~astropy.visualization.wcsaxes.WCSAxes`
            WCS axes object.
        """
        center_px = cast(
            Tuple[float, float], sky_position_to_pixel(self._center, self._wcs)
        )
        ra_pixel_scale, _ = pixel_scales(self._wcs)
        radius_px = float(self._radius / ra_pixel_scale)
        circle_patch = Circle(center_px, radius_px, **self._kwargs)
        ax.add_patch(circle_patch)

    def rotate(self, pivot: SkyCoord, angle: Angle) -> "CircleAnnotation":
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
        `~imephu.annotation.general.CircleAnnotation`
            The annotation resulting from the rotation.
        """
        rotated_annotation = deepcopy(self)
        rotated_center = rotate(self._center, pivot, angle, self._wcs)
        rotated_annotation._center = rotated_center
        return rotated_annotation

    def translate(self, displacement: Angle) -> "CircleAnnotation":
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.general.CircleAnnotation`
            The annotation resulting from the translation.
        """
        translated_annotation = deepcopy(self)
        translated_center = translate(self._center, displacement)
        translated_annotation._center = translated_center
        return translated_annotation
