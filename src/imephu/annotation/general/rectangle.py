from typing import Any, List, cast

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.units import Quantity
from astropy.wcs import WCS

from imephu.annotation.general.line_path import LinePathAnnotation
from imephu.geometry import translate


class RectangleAnnotation(LinePathAnnotation):
    """An annotation for plotting a rectangle on a finder chart.

    The rectangle is defined by its center, width and height. The center must be a sky
    position in right ascension and declination, the width and height must be an angle
    pm the sky. The width is taken to be the extension in the right ascension direction,
    the height the extension in the declination direction.

    The corners of the rectangle are mapped to their corresponding pixel positions on
    the finder chart. This implies that in general the figure plotted on the finder
    chart will be a quadrilateral, but not exactly a rectangle.

    .. note::
        The center rather than "bottom left" corner is specified, as this is more
        convenient in the context of finder charts.

    Parameters
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The center of the rectangle, in right ascension and declination.
    width: `~astropy.units.Quantity` ["angle"]
        The width of the rectangle as an angle on the sky in the direction of the
        right ascension.
    height: ~astropy.units.Quantity` ["angle"]
        The height of the rectangle as an angle on the sky in the direction of the
        declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    edgecolor: str, default: "black"
        The edge color.
    facecolor: str, default: "none"
        The filling color.
    **kwargs: dict, optional
        Additional keyword arguments, which will be passed to Matplotlib's
        `~matplotlib.patches.PathPatch` patch constructor.
    """

    def __init__(
        self,
        center: SkyCoord,
        width: Quantity,
        height: Quantity,
        wcs: WCS,
        edgecolor: str = "black",
        facecolor: str = "none",
        **kwargs: Any,
    ):
        super().__init__(
            vertices=RectangleAnnotation._corners(center, width, height),
            wcs=wcs,
            closed=True,
            edgecolor=edgecolor,
            facecolor=facecolor,
            **kwargs,
        )

    @staticmethod
    def _corners(center: SkyCoord, width: Quantity, height: Quantity) -> List[SkyCoord]:
        corners: List[SkyCoord] = []
        half_width_arcsec = width.to_value(u.arcsec) / 2
        half_height_arcsec = height.to_value(u.arcsec) / 2
        directions = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        for direction in directions:
            center_to_corner = (
                np.array(
                    (
                        direction[0] * half_width_arcsec,
                        direction[1] * half_height_arcsec,
                    )
                )
                * u.arcsec
            )
            corner = translate(center, center_to_corner)
            corners.append(corner)
        return corners

    def rotate(self, pivot: SkyCoord, angle: Angle) -> "RectangleAnnotation":
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
        `~imephu.annotation.general.RectangleAnnotation`
            The annotation resulting from the rotation.
        """
        return cast(RectangleAnnotation, super().rotate(pivot, angle))

    def translate(self, displacement: Quantity) -> "RectangleAnnotation":
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.general.RectangleAnnotation`
            The annotation resulting from the translation.
        """
        return cast(RectangleAnnotation, super().translate(displacement))
