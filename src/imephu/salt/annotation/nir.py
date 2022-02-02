from typing import List

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import LinePathAnnotation, GroupAnnotation, \
    CrosshairsAnnotation, CircleAnnotation
from imephu.geometry import translate


def nir_annotation(center: SkyCoord, wcs: WCS) -> GroupAnnotation:
    return _science_fiber_bundle_annotation(center, wcs)


def _science_fiber_bundle_annotation(center: SkyCoord, wcs: WCS) -> GroupAnnotation:
    return fiber_bundle_annotation(center=center, min_width=22 * u.arcsec, max_width=29 * u.arcsec, max_height=18 * u.arcsec, wcs=wcs)


def fiber_bundle_annotation(center: SkyCoord, min_width: Angle, max_width: Angle, max_height: Angle, wcs: WCS) -> GroupAnnotation:
    """
    Return the annotation for an NIR fiber bundle.

    The bundle is assumed to be a flattened hexagon, which is characterized by a minimum
    width ("at the top and bottom"), a maximum width ("in the middle") and a maximum
    height ("in the middle"). Width and height are measured in the direction of right
    ascension and declination, respectively.

    The outline of this hexagon is shown together with all the fibers. At the center of
    the hexagon crosshairs are added.

    Parameters
    ----------
    center: `astropy.coordinates.SkyCoord`
        The center of the bundle, as a position on the sky, in right ascension and
        declination.
    min_width: `~astropy.coordinates.Angle`
        The minimum bundle width (in the direction of right ascension), as an angle on
        the sky.
    max_width: `~astropy.coordinates.Angle`
        The maximum bundle width (in the direction of right ascension), as an angle on
        the sky.
    max_height: `~astropy.coordinates.Angle`
        The maximum bundle height (in the direction of declination), as an angle on the
        sky.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `imephu.annotation.general.GroupAnnotation`
        The annotation for the fiber bundle.
    """
    bundle_annotation = GroupAnnotation([])
    bundle_annotation.add_item(_fiber_bundle_outline(center=center, min_width=min_width, max_width=max_width, height=max_height, wcs=wcs))
    bundle_annotation.add_item(CrosshairsAnnotation(center=center, size=10 * u.arcsec, wcs=wcs, color="red"))
    return bundle_annotation


def _fiber_bundle_outline(center: SkyCoord, min_width: Angle, max_width: Angle, height: Angle, wcs: WCS) -> LinePathAnnotation:
    min_width_arcsec = min_width.to_value(u.arcsec)
    max_width_arcsec = max_width.to_value(u.arcsec)
    height_arcsec = height.to_value(u.arcsec)
    vertices = [
        translate(center, (-min_width_arcsec / 2, -height_arcsec / 2) * u.arcsec),
        translate(center, (min_width_arcsec / 2, -height_arcsec / 2) * u.arcsec),
        translate(center, (max_width_arcsec / 2, 0) * u.arcsec),
        translate(center, (min_width_arcsec / 2, height_arcsec / 2) * u.arcsec),
        translate(center, (-min_width_arcsec / 2, height_arcsec / 2) * u.arcsec),
        translate(center, (-max_width_arcsec / 2, 0) * u.arcsec)
        ]
    return LinePathAnnotation(vertices=vertices, wcs=wcs, edgecolor="red")
