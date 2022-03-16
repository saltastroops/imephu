from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import CircleAnnotation, GroupAnnotation, TextAnnotation
from imephu.geometry import translate


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
    label_position = translate(fits_center, (-3.6, 3.6) * u.arcmin)
    name_annotation = TextAnnotation(
        label_position,
        "SCAM",
        wcs=wcs,
        style="italic",
        weight="bold",
        size="large",
        horizontalalignment="left",
        color=(0, 0, 1),
    )
    return GroupAnnotation([fov_annotation, name_annotation])
