from typing import Any, Tuple, Union

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import GroupAnnotation, LinePathAnnotation
from imephu.geometry import translate


class CrosshairsAnnotation(GroupAnnotation):
    """An annotation for plotting crosshairs.

    The crosshairs are defined by their center (i.e. where the two lines intersect), the
    size (i.e. the length of the lines) and the angle of rotation around the center.

    The lines bisect each other. The crosshair lines are in the direction of right
    ascension and declination.

    Parameters
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The center of the crosshairs, i.e. the point of intersection of the two lines,
        as a position on the sky (in right ascension and declination).
    size: `~astropy.coordinates.Angle`
        The length of each line, as an angle on the sky.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    color: color, default: "black"
        The line color.
    **kwargs: dict, optional
        Additional keyword arguments, which will be passed to Matplotlib's
        `~matplotlib.patches.PathPatch` patch constructor.
    """

    def __init__(
        self,
        center: SkyCoord,
        size: Angle,
        wcs: WCS,
        color: Union[str, Tuple[float, float, float]] = "black",
        **kwargs: Any,
    ) -> None:
        super().__init__(
            items=CrosshairsAnnotation._crosshair_lines(
                center=center, size=size, wcs=wcs, color=color, **kwargs
            )
        )

    @staticmethod
    def _crosshair_lines(
        center: SkyCoord,
        size: Angle,
        wcs: WCS,
        color: Union[str, Tuple[float, float, float]] = "black",
        **kwargs: Any,
    ) -> Tuple[LinePathAnnotation, LinePathAnnotation]:
        size_arcsec = size.to_value(u.arcsec)
        line_ra = LinePathAnnotation(
            vertices=[
                translate(center, Angle((-size_arcsec / 2, 0) * u.arcsec)),
                translate(center, Angle((size_arcsec / 2, 0) * u.arcsec)),
            ],
            wcs=wcs,
            edgecolor=color,
            **kwargs,
        )
        line_dec = LinePathAnnotation(
            vertices=[
                translate(center, Angle((0, -size_arcsec / 2) * u.arcsec)),
                translate(center, Angle((0, size_arcsec / 2) * u.arcsec)),
            ],
            wcs=wcs,
            edgecolor=color,
            **kwargs,
        )
        return line_ra, line_dec
