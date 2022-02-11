import math
from math import atan2, cos, sin
from typing import Any, Tuple, cast

import numpy as np
import numpy.typing as npt
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel


def rotate(v: SkyCoord, pivot: SkyCoord, angle: Angle, wcs: WCS) -> SkyCoord:
    """Rotate a point around a pivot.

    Both the pivot and the point to rotate are assumed to be given as right ascension
    and declination. The rotation angle is taken to be the angle on the sky. If the
    pixel scales for right ascension and declination differ, the angle seen on the
    finder chart will differ from this angle.

    The angle is measured from north to south. Hence it depends on the orientation of
    the coordinate axes whether a positive angle corresponds to a clockwise or
    anti-clockwise orientation.

    Parameters
    ----------
    v: `~astropy.coordinates.SkyCoord`
        The point to rotate.
    pivot: `~astropy.coordinates.SkyCoord`
        The point around which the point ``v`` is rotated.
    angle: `~astropy.coordinates.Angle`
        The angle of rotation, measured from north to east.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        The point after the rotation.
    """
    # Convert from world coordinates to pixels
    v_px = sky_position_to_pixel(v, wcs)
    pivot_px = sky_position_to_pixel(pivot, wcs)

    # Correct the angle for different pixel scales
    sx, sy = pixel_scales(wcs)
    angle_rad = angle.to_value(u.rad)
    x_angle_sky = cos(angle_rad)
    y_angle_sky = sin(angle_rad)
    x_angle_px = x_angle_sky
    y_angle_px = float(sx / sy) * y_angle_sky
    angle_px = atan2(y_angle_px, x_angle_px)

    # If necessary, change the sign of the angle to accommodate the orientation of the
    # axes on the finder chart.
    if not _is_positive_angle_anti_clockwise(wcs):
        angle_px *= -1

    # Perform the rotation
    rotation: Any = np.array(
        [[cos(angle_px), -sin(angle_px)], [sin(angle_px), cos(angle_px)]]
    )
    diff_vector_px: Any = np.array(v_px) - np.array(pivot_px)
    rotated_v_px = np.array(pivot_px) + rotation @ diff_vector_px

    # Convert back from pixels to world coordinates
    return pixel_to_sky_position(rotated_v_px, wcs)


def translate(v: SkyCoord, displacement: Angle) -> SkyCoord:
    """Move a point by a displacement vector.

    The point to move is assumed to be given as right ascension and declination. The
    displacement vector is assumed to be a vector on the sky, with real angles, as
    described for the `~sky_vector_to_pixel` method.

    Parameters
    ----------
    v: `~astropy.coordinates.SkyCoord`
        Point to move.
    displacement: 2D array of angles
        Vector by which to move the point.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        The point after moving.
    """
    # Convert from sky coordinates to real angles
    displacement_delta_ra = displacement[0] / cos(v.dec.to_value(u.rad))
    displacement_delta_dec = displacement[1]

    # Perform the translation
    translated_v_ra = v.ra + displacement_delta_ra
    translated_v_dec = v.dec + displacement_delta_dec
    return SkyCoord(ra=translated_v_ra, dec=translated_v_dec)


def pixel_scales(wcs: WCS) -> Tuple[Angle, Angle]:
    """
    Calculate the pixel scales of a WCS object.

    For a WCS returns pixel scales along each axis of the image pixel at
    the ``CRPIX`` location once it is projected onto the
    "plane of intermediate world coordinates" as defined in
    `Greisen & Calabretta 2002, A&A, 395, 1061
    <https://ui.adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_.

    .. note::
        This function is concerned **only** about the transformation
        "image plane"->"projection plane" and **not** about the
        transformation "celestial sphere"->"projection plane"->"image plane".
        Therefore, this function ignores distortions arising due to
        non-linear nature of most projections.

    .. note::
        In order to compute the scales corresponding to celestial axes only,
        make sure that the input `~astropy.wcs.WCS` object contains
        celestial axes only, e.g., by passing in the
        `~astropy.wcs.WCS.celestial` WCS object.

    The code has been adapted from APLpy.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        A world coordinate system object.

    Returns
    -------
    scale : tuple of `~astropy.coordinates.Angle`
        A vector (`~numpy.ndarray`) of projection plane increments
        corresponding to each pixel side (axis). The units of the returned
        results are the same as the units of `~astropy.wcs.Wcsprm.cdelt`,
        `~astropy.wcs.Wcsprm.crval`, and `~astropy.wcs.Wcsprm.cd` for
        the celestial WCS and can be obtained by inquiring the value
        of `~astropy.wcs.Wcsprm.cunit` property of the input
        `~astropy.wcs.WCS` WCS object.

    See Also
    --------
    astropy.wcs.utils.proj_plane_pixel_area

    """
    scale_factors = np.sqrt((wcs.pixel_scale_matrix**2).sum(axis=0, dtype=float))
    return (
        Angle(scale_factors[0], wcs.world_axis_units[0]),
        Angle(scale_factors[1], wcs.world_axis_units[1]),
    )


def sky_position_to_pixel(position: SkyCoord, wcs: WCS) -> npt.NDArray[np.float_]:
    """Convert a sky position to the corresponding pixel coordinates.

    Parameters
    ----------
    position: `~astropy.coordinates.SkyCoord`
        Sky position as right ascension and declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~numpy.ndarray` of `float`
        The pixel coordinates.
    """
    pixel_coords = skycoord_to_pixel(coords=position, wcs=wcs)
    return np.array([float(pixel_coords[0].item()), float(pixel_coords[1].item())])


def pixel_to_sky_position(position_px: npt.NDArray[np.float_], wcs: WCS) -> SkyCoord:
    """Convert pixel coordinates to the corresponding sky position.

    Parameters
    ----------
    position_px: `~numpy.ndarray` of `float`
        Pixel coordinates.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        The sky coordinates.
    """
    return pixel_to_skycoord(position_px[0], position_px[1], wcs=wcs)  # noqa


def _is_positive_angle_anti_clockwise(wcs: WCS) -> bool:
    """
    Return whether for a given WCS positive angles are anti-clockwise when plotted.

    Angles on the sky are measured from north to east, and depending on the orientation
    of the axes on the finder chart this may correspond to a mathematically positive
    ("anti-clockwise") or negative ("clockwise") angle.

    Parameters
    ----------
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    bool
        ``True`` if angles from north to east correspond to mathematically positive
        ("anti-clockwise") angles on a finder chart, ``False`` otherwise.
    """
    origin_px: npt.NDArray[np.float_] = np.array((0.0, 0.0))
    origin_on_sky = pixel_to_sky_position(origin_px, wcs)
    origin_ra = origin_on_sky.ra
    origin_dec = origin_on_sky.dec

    # A "point to the north" on a finder chart
    point_to_north_ra = origin_ra
    point_to_north_dec = origin_dec + 1 * u.arcmin
    point_to_north = SkyCoord(ra=point_to_north_ra, dec=point_to_north_dec)
    point_to_north_px = sky_position_to_pixel(point_to_north, wcs)

    # A "point to the north-east" on a finder chart
    point_to_north_east_ra = origin_ra + 1 * u.arcmin
    point_to_north_east_dec = origin_dec + 1 * u.arcmin
    point_to_north_east = SkyCoord(
        ra=point_to_north_east_ra, dec=point_to_north_east_dec
    )
    point_to_north_east_px = sky_position_to_pixel(point_to_north_east, wcs)

    # A "point to the north-east" on a finder chart
    point_to_north_west_ra = origin_ra - 1 * u.arcmin
    point_to_north_west_dec = origin_dec + 1 * u.arcmin
    point_to_north_west = SkyCoord(
        ra=point_to_north_west_ra, dec=point_to_north_west_dec
    )
    point_to_north_west_px = sky_position_to_pixel(point_to_north_west, wcs)

    # The difference between the point to the north and the origin gives a vector
    # pointing north - but as the origin is (0, 0) this is just the point to the north.
    # The same argument applies for the vector pointing north-east.
    north_vector_px = point_to_north_px
    north_east_vector_px = point_to_north_east_px
    north_west_vector_px = point_to_north_west_px

    # Get the angle between the vectors to the north and to the north-east. Due to
    # symmetry, this is the same as the angle between the vectors to the north and to
    # north-west
    cos_angle = np.dot(north_vector_px, north_east_vector_px) / (
        np.linalg.norm(north_vector_px) * np.linalg.norm(north_east_vector_px)
    )
    angle = math.acos(cos_angle)

    # Rotate the vector to the north (in anti-clockwise direction) by this angle
    rotation: npt.NDArray[np.float_] = np.array(
        [[cos(angle), -sin(angle)], [sin(angle), cos(angle)]]
    )
    rotated_north_vector_px = np.dot(rotation, north_vector_px)

    # If a rotation from north to east corresponds to an anti-clockwise rotation on a
    # finder chart, the rotated vector must aligned with the vector to the north-east.
    # This implies that their scalar product then must be larger than that between the
    # rotated vector and the vector to the north-west.
    return cast(
        bool,
        np.dot(rotated_north_vector_px, north_east_vector_px)
        > np.dot(rotated_north_vector_px, north_west_vector_px),
    )
