from typing import List

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import LinePathAnnotation, GroupAnnotation, \
    CrosshairsAnnotation, CircleAnnotation
from imephu.geometry import translate


def nir_annotation(science_bundle_center: SkyCoord, bundle_separation: Angle, position_angle: Angle, wcs: WCS) -> GroupAnnotation:
    """Return the annotation for an NIR setup.

    The annotation contains the outline of the science fiber bundle (with its center
    marked by crosshairs) and the outlines of the sky fiber bundles.

    Parameters
    ----------
    science_bundle_center: `~astropy.coordinates.SkyCoord`
        The center position of the sky fiber bundle, as a position on the sky, in
        right ascension and declination.
    bundle_separation: ~astropy.coordinates.Angle`
        The separation between the science fiber bundle and the sky fiber bundles, as an
        angle on the sky. The separation is measured between the center of the science
        bundle and the midpoint of the line between the centers of the sky bundles.
    position_angle: `~astropy.coordinates.Angle`
        The position angle, as an angle on the sky measured from north to east.
    wcs:`astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `imephu.annotation.general.GroupAnnotation`
        The annotation for an NIR setup.
    """
    science_bundle = _science_fiber_bundle_annotation(science_bundle_center, wcs)

    bundle_separation_arcsec = bundle_separation.to_value(u.arcsec)
    sky_bundle_separation = 20 * u.arcsec
    sky_bundle_separation_arcsec = sky_bundle_separation.to_value(u.arcsec)
    sky_bundle1_center = translate(science_bundle_center, (-sky_bundle_separation_arcsec / 2, bundle_separation_arcsec) * u.arcsec)
    sky_bundle1 = _sky_fiber_bundle_annotation(sky_bundle1_center, wcs)

    sky_bundle2_center = translate(science_bundle_center, (sky_bundle_separation_arcsec / 2, bundle_separation_arcsec) * u.arcsec)
    sky_bundle2 = _sky_fiber_bundle_annotation(sky_bundle2_center, wcs)

    annotation = GroupAnnotation([science_bundle, sky_bundle1, sky_bundle2])
    return annotation.rotate(science_bundle_center, position_angle)


def _science_fiber_bundle_annotation(center: SkyCoord, wcs: WCS) -> GroupAnnotation:
    return fiber_bundle_annotation(center=center, min_width=22 * u.arcsec, max_width=29 * u.arcsec, max_height=18 * u.arcsec, with_crosshairs=True, wcs=wcs)


def _sky_fiber_bundle_annotation(center: SkyCoord, wcs: WCS) -> GroupAnnotation:
    bundle_annotation = fiber_bundle_annotation(center, min_width=10 *u.arcsec, max_width=13 * u.arcsec, max_height=8 * u.arcsec, with_crosshairs=False, wcs=wcs)
    return bundle_annotation.rotate(center, 90 * u.deg)



def fiber_bundle_annotation(center: SkyCoord, min_width: Angle, max_width: Angle, max_height: Angle, with_crosshairs: bool, wcs: WCS) -> GroupAnnotation:
    """
    Return the annotation for an NIR fiber bundle.

    The bundle is assumed to be a flattened hexagon, which is characterized by a minimum
    width ("at the top and bottom"), a maximum width ("in the middle") and a maximum
    height ("in the middle"). Width and height are measured in the direction of right
    ascension and declination, respectively.

    The outline of this hexagon is shown together with all the fibers. At the center of
    the hexagon optionally crosshairs are added.

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
    with_crosshairs: `bool`
        Whether to add crosshairs at the center.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `imephu.annotation.general.GroupAnnotation`
        The annotation for the fiber bundle.
    """
    bundle_annotation = GroupAnnotation([])
    bundle_annotation.add_item(_fiber_bundle_outline(center=center, min_width=min_width, max_width=max_width, height=max_height, wcs=wcs))
    if with_crosshairs:
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
