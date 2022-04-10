import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import (
    CircleAnnotation,
    CrosshairsAnnotation,
    GroupAnnotation,
    LinePathAnnotation,
)
from imephu.geometry import translate

_OUTLINE_COLOR = "red"


def bundles_annotation(
    science_bundle_center: SkyCoord,
    bundle_separation: Angle,
    position_angle: Angle,
    wcs: WCS,
    include_fibers: bool = False,
) -> GroupAnnotation:
    """Return the annotation for NIR fiber bundles.

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
    include_fibers: `bool`, optional, default False
        Whether to include the fiber outlines in the annotation.

    Returns
    -------
    `imephu.annotation.general.GroupAnnotation`
        The annotation for an NIR setup.
    """
    r = 0.65 * u.arcsec
    d_horizontal = 1.63 * u.arcsec
    d_vertical = 1.42 * u.arcsec
    bundle_separaration_arcsec = bundle_separation.to_value(u.arcsec)
    science_bundle = _science_bundle(
        r=r,
        d_horizontal=d_horizontal,
        d_vertical=d_vertical,
        center=science_bundle_center,
        wcs=wcs,
        include_fibers=include_fibers,
    )
    sky_bundle = _sky_bundle(
        r=r,
        d_horizontal=d_horizontal,
        d_vertical=d_vertical,
        center=science_bundle_center,
        wcs=wcs,
        include_fibers=include_fibers,
    ).translate((0, bundle_separaration_arcsec) * u.arcsec)
    annotation = GroupAnnotation([sky_bundle, science_bundle])
    return annotation.rotate(science_bundle_center, position_angle)


def _science_bundle(
    r: Angle,
    d_horizontal: Angle,
    d_vertical: Angle,
    center: SkyCoord,
    wcs: WCS,
    include_fibers: bool,
) -> GroupAnnotation:
    """
    Return an annotation showing the science bundle.

    The annotation includes the science bundle outline with crosshairs at the center
    and, optionally, the bundle fiber outlines.

    Parameters
    ----------
    r: `~astropy.coord.Angle`
        The fiber radius, as angle on the sky.
    d_horizontal: `~astropy.coord.Angle`
        The distance between adjacent fibers in the direction of right ascension, as an
        angle on the sky.
    d_vertical: `~astropy.coord.Angle`
        The distance between adjacent rows of fiber cdenters in the direction of
        declination, as an angle on the sky.
    center: `~astropy.coord.SkyCoords`
        The position of the bundle center, as right ascension and declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    include_fibers: `bool`
        Whether to include the fiber outlines in the annotation.

    Returns
    -------
    `~imephu.annotation.general.GroupAnnotation`
        The annotation.
    """
    annotation = GroupAnnotation([])
    if include_fibers:
        annotation.add_item(
            _science_fiber_outlines(
                r=r,
                d_horizontal=d_horizontal,
                d_vertical=d_vertical,
                center=center,
                wcs=wcs,
            )
        )
    annotation.add_item(
        _science_bundle_outline(
            r=r,
            d_horizontal=d_horizontal,
            d_vertical=d_vertical,
            center=center,
            wcs=wcs,
        )
    )
    annotation.add_item(
        CrosshairsAnnotation(
            center=center, size=10 * u.arcsec, wcs=wcs, color=_OUTLINE_COLOR
        )
    )
    return annotation


def _science_bundle_outline(
    r: Angle,
    d_horizontal: Angle,
    d_vertical: Angle,
    center: SkyCoord,
    wcs: WCS,
    include_fibers: bool = False,
) -> LinePathAnnotation:
    # For simplicity, we first assume that the bundle center is located at (0, 0).
    # We start by establishing one of the slanting lines. We make use id the fact that
    # the line is parallel to the line through the centers of the fibers to which it is
    # tangential.
    dh = d_horizontal.to_value(u.arcsec)
    dv = d_vertical.to_value(u.arcsec)
    fiber_center = np.array([-8.5 * dh, 2 * dv]) * u.arcsec

    # m is the slope of the line. The vector n from the fiber center to the line is
    # perpendicular to the line.
    m = dv / (dh / 2)
    dr = np.array([-dv, dh / 2]) * u.arcsec
    n = dr / np.linalg.norm(dr)

    # The distance between the fiber center and the line is the fiber radius, and with
    # this in mind we can calculate the y offset c of the line.
    p = fiber_center + r * n
    c = p[1] - m * p[0]

    # The first end point of the line is the intersection of the adjacent horizontal
    # line and y = m x + c. In the following, gap_vertical is the vertical distance
    # between circles in neighbouring rows.
    gap_vertical = d_vertical - 2 * r
    y1 = 6.5 * d_vertical - gap_vertical / 2
    x1 = (y1 - c) / m
    vertices = []
    vertices.append((x1, y1))

    # The second end point of the line is the intersection of the adjacent vertical line
    # and y = m x + c. In the following, gap_horizontal is the horizontal distance
    # between adjacent circles.
    gap_horizontal = d_horizontal - 2 * r
    x2 = -9 * d_horizontal + gap_horizontal / 2
    y2 = m * x2 + c
    vertices.append((x2, y2))

    # The next four vertices are just end points of horizontal and vertical lines.
    x3 = x2
    y3 = 2 * d_vertical - r
    vertices.append((x3, y3))

    x4 = -8 * d_horizontal - r
    y4 = y3
    vertices.append((x4, y4))

    x5 = x4
    y5 = r
    vertices.append((x5, y5))

    x6 = x2
    y6 = y5
    vertices.append((x6, y6))

    # And that's one quarter of the bundle outline covered. The remaining vertices can
    # be obtained by changing signs appropriately.
    quarter = vertices[:]
    quarter_reversed = vertices[::-1]
    vertices.extend((v[0], -v[1]) for v in quarter_reversed)
    vertices.extend((-v[0], -v[1]) for v in quarter)
    vertices.extend((-v[0], v[1]) for v in quarter_reversed)

    # So far we have assumed the bundle center to be at (0, 0), so we have to correct
    # this.
    vertex_positions = [translate(center, v) for v in vertices]

    line = LinePathAnnotation(
        vertex_positions, wcs=wcs, closed=True, edgecolor=_OUTLINE_COLOR
    )
    return line


def _science_fiber_outlines(
    r: Angle, d_horizontal: Angle, d_vertical: Angle, center: SkyCoord, wcs: WCS
) -> GroupAnnotation:
    annotation = GroupAnnotation([])
    for row in range(-6, 7):
        if row == -6 or row == 6:
            values = range(-7, 7)
        elif row == -5 or row == 5:
            values = range(-8, 7)
        elif row == -4 or row == 4:
            values = range(-8, 8)
        elif row == -3 or row == 3:
            values = range(-9, 8)
        elif row == -2 or row == 2:
            values = range(-9, 9)
        elif row == -1 or row == 1:
            values = range(-9, 8)
        else:
            values = range(-9, 9)
        for v in values:
            x = (
                (v + 0.5) * d_horizontal.to_value(u.arcsec)
                if row % 2 == 0
                else (v + 1) * d_horizontal.to_value(u.arcsec)
            )
            y = row * d_vertical.to_value(u.arcsec)
            circle = CircleAnnotation(
                center=center, radius=r, wcs=wcs, edgecolor="darkgray", facecolor="none"
            ).translate((x, y) * u.arcsec)
            annotation.add_item(circle)
    return annotation


def _sky_bundle(
    r: Angle,
    d_horizontal: Angle,
    d_vertical: Angle,
    center: SkyCoord,
    wcs: WCS,
    include_fibers: bool,
) -> GroupAnnotation:
    """
    Return an annotation showing the sky bundle.

    The annotation includes the sky bundle outline and, optionally, the bundle fiber
    outlines.

    Parameters
    ----------
    r: `~astropy.coord.Angle`
        The fiber radius, as angle on the sky.
    d_horizontal: `~astropy.coord.Angle`
        The distance between adjacent fibers in the direction of right ascension, as an
        angle on the sky.
    d_vertical: `~astropy.coord.Angle`
        The distance between adjacent rows of fiber cdenters in the direction of
        declination, as an angle on the sky.
    center: `~astropy.coord.SkyCoords`
        The position of the bundle center, as right ascension and declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    include_fibers: `bool`
        Whether to include the fiber outlines in the annotation.

    Returns
    -------
    `~imephu.annotation.general.GroupAnnotation`
        The annotation.
    """
    annotation = GroupAnnotation([])
    if include_fibers:
        annotation.add_item(
            _sky_fiber_outlines(
                r=r,
                d_horizontal=d_horizontal,
                d_vertical=d_vertical,
                center=center,
                wcs=wcs,
            )
        )
    annotation.add_item(
        _sky_bundle_outline(
            r=r,
            d_horizontal=d_horizontal,
            d_vertical=d_vertical,
            center=center,
            wcs=wcs,
        )
    )
    return annotation


def _sky_bundle_outline(
    r: Angle, d_horizontal: Angle, d_vertical: Angle, center: SkyCoord, wcs: WCS
) -> LinePathAnnotation:
    # For simplicity, we first assume that the bundle center is located at (0, 0).
    vertices = [
        (-8.5 * d_horizontal - r, r),
        (-7 * d_horizontal - r, r),
        (-7 * d_horizontal - r, d_vertical + r),
        (-4 * d_horizontal + r, d_vertical + r),
        (-4 * d_horizontal + r, r),
        (3 * d_horizontal - r, r),
        (3 * d_horizontal - r, d_vertical + r),
        (8 * d_horizontal + r, d_vertical + r),
        (8 * d_horizontal + r, r),
        (8.5 * d_horizontal + r, r),
    ]

    # And that's one half of the bundle outline covered. The remaining vertices can be
    # obtained by changing signs appropriately.
    half = vertices[:]
    vertices.extend((-v[0], -v[1]) for v in half)

    # So far we have assumed the bundle center to be at (0, 0), so we have to correct
    # this.
    vertex_positions = [translate(center, v) for v in vertices]

    line = LinePathAnnotation(
        vertex_positions, wcs=wcs, closed=True, edgecolor=_OUTLINE_COLOR
    )
    return line


def _sky_fiber_outlines(
    r: Angle, d_horizontal: Angle, d_vertical: Angle, center: SkyCoord, wcs: WCS
) -> GroupAnnotation:
    annotation = GroupAnnotation([])

    for v in list(range(-7, -3)) + list(range(3, 9)):
        x = v * d_horizontal.to_value(u.arcsec)
        y = d_vertical.to_value(u.arcsec)
        circle = CircleAnnotation(
            center=center, radius=r, wcs=wcs, edgecolor="darkgray", facecolor="none"
        ).translate((x, y) * u.arcsec)
        annotation.add_item(circle)

    for v in list(range(-9, 9)):
        x = (v + 0.5) * d_horizontal.to_value(u.arcsec)
        y = 0
        circle = CircleAnnotation(
            center=center, radius=r, wcs=wcs, edgecolor="darkgray", facecolor="none"
        ).translate((x, y) * u.arcsec)
        annotation.add_item(circle)

    for v in list(range(-8, -2)) + list(range(4, 8)):
        x = v * d_horizontal.to_value(u.arcsec)
        y = -d_vertical.to_value(u.arcsec)
        circle = CircleAnnotation(
            center=center, radius=r, wcs=wcs, edgecolor="darkgray", facecolor="none"
        ).translate((x, y) * u.arcsec)
        annotation.add_item(circle)

    return annotation
