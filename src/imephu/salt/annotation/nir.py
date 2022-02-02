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
    fiber_diameter = 1.3 * u.arcsec
    rows = 13
    columns = 18
    max_row_width = 29 * u.arcsec
    max_column_height = 18 * u.arcsec
    fiber_distance_in_ra = (max_row_width - fiber_diameter) / (columns - 1)
    fiber_distance_in_dec = (max_column_height - fiber_diameter) / (rows - 1)

    fiber_centers = []
    for row in range(0, rows):
        if row in (0, rows - 1):
            index_range = range(2, columns - 2)
        elif row in (1, rows - 2):
            index_range = range(1, columns - 2)
        elif row in (2, rows - 3):
            index_range = range(1, columns - 1)
        elif row % 2 == 0:
            index_range = range(0, columns)
        else:
            index_range = range(0, columns - 1)

        for column in index_range:
            displacement_ra = (column - columns // 2) * fiber_distance_in_ra
            if row % 2 == 0:
                displacement_ra -= 0.5 * fiber_distance_in_ra
            displacement_dec = (row - rows // 2) * fiber_distance_in_dec
            displacement_ra_arcsec = displacement_ra.to_value(u.arcsec)
            displacement_dec_arcsec = displacement_dec.to_value(u.arcsec)
            fiber_center = translate(center, (displacement_ra_arcsec, displacement_dec_arcsec) * u.arcsec)
            fiber_centers.append(fiber_center)

    return fiber_bundle_annotation(center=center, min_width=7 * u.arcmin, max_width=8 * u.arcmin, max_height=8 * u.arcmin, fiber_centers=fiber_centers, fiber_diameter=fiber_diameter, wcs=wcs)


def fiber_bundle_annotation(center: SkyCoord, min_width: Angle, max_width: Angle, max_height: Angle, fiber_centers: List[SkyCoord], fiber_diameter: Angle, wcs: WCS) -> GroupAnnotation:
    """
    Return the annotation for an NIR fiber bundle.

    The bundle is assumed to be a flattened hexagon, which is characterized by a minimum
    width ("at the top and bottom"), a maximum width ("in the middle") and a maximum
    height ("in the middle"). Width and height are measured in the direction of right
    ascension and declination, respectively.

    The outline of this hexagon is shown together with all the fibers. At the center of
    the hexagon crosshairs are added, the size of which equals a fiber's diameter. The
    outline of all the fibers is shown as well.

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
    fiber_centers: list of `astropy.coordinates.SkyCoord`
        The list of the center positions of all the fibers in the bundle, as positions
        on the sky, in right ascension and declination.
    fiber_diameter: `~astropy.coordinates.Angle`
        The diameter of a fiber, as an angle on the sky.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `imephu.annotation.general.GroupAnnotation`
        The annotation for the fiber bundle.
    """
    bundle_annotation = GroupAnnotation([])
    bundle_annotation.add_item(_fiber_bundle_outline(center=center, min_width=min_width, max_width=max_width, height=max_height, wcs=wcs))
    bundle_annotation.add_item(CrosshairsAnnotation(center=center, size=fiber_diameter, wcs=wcs, color="red"))
    for fiber in fiber_centers:
        bundle_annotation.add_item(CircleAnnotation(center=fiber, radius=fiber_diameter / 2, wcs=wcs, edgecolor="none", facecolor="yellow", alpha=0.2))
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
