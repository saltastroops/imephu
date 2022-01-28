from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import (
    CircleAnnotation,
    GroupAnnotation,
    RectangleAnnotation,
    TextAnnotation,
)


def field_of_view_annotation(fits_center: SkyCoord, wcs: WCS) -> GroupAnnotation:
    """
    Return an annotation with the Salticam field of view.

    Parameters
    ----------
    fits_center: `~astropy.coordinates.SkyCoord`
        The central position of the finder chart, in right ascension and declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.GroupAnnotation`
        An annotation with the Salticam field of view.
    """
    fov_annotation = CircleAnnotation(
        fits_center, 5 * u.arcmin, wcs=wcs, edgecolor="green"
    )
    name_annotation = TextAnnotation(
        (0.86, 0.86),
        "SCAM",
        wcs=wcs,
        style="italic",
        weight="bold",
        size="large",
        horizontalalignment="left",
        color=(0, 0, 1),
    )
    return GroupAnnotation([fov_annotation, name_annotation])


def slot_annotation(center: SkyCoord, position_angle: Angle, wcs: WCS):
    """Return an annotation showing the slot.

    Parameters
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The center of the slot on the sky, in right ascension and declination.
    position_angle: `~astropy.coordinates.Angle`
        The position angle as an angle on the sky, measured from north to east.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.RectangleAnnotation`
        The slot annotation.
    """
    return RectangleAnnotation(
        center=center,
        width=1 * u.arcmin / 3,
        height=10 * u.arcmin,
        wcs=wcs,
        edgecolor="red",
        alpha=0.5,
        linewidth=2,
    ).rotate(center, position_angle + 90 * u.deg)
