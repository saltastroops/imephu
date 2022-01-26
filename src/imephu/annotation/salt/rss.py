from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import RectangleAnnotation


def longslit_annotation(
    fits_center: SkyCoord,
    slit_width: Angle,
    slit_height: Angle,
    position_angle: Angle,
    wcs: WCS,
) -> RectangleAnnotation:
    """Return an annotation showing a longslit.

    Parameters
    ----------
    fits_center: `~astropy.coordinates.SkyCoord`
        The central position of the finder chart, in right ascension and declination
    slit_width: `~astropy.coordinates.Angle`
        The width of the slit, as an angle on the sky.
    slit_height: `~astropy.coordinates.Angle`
        The height of the slit, as an angle on the sky.
    position_angle: `~astropy.coordinates.Angle`
        The position angle, as an angle on the sky measured from north to east.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `imephu.annotation.general.RectangleAnnotation`
        The annotation showing the longslit.
    """
    return RectangleAnnotation(
        fits_center, slit_width, slit_height, wcs=wcs, edgecolor="red", alpha=0.5
    ).rotate(fits_center, position_angle)
