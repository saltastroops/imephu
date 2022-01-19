from math import atan2, cos, sin
from typing import Tuple

import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
from astropy.units import Quantity
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord

AngleLike = Quantity

SkyPosition = Tuple[AngleLike, AngleLike]

SkyVector = Tuple[AngleLike, AngleLike]

PixelPosition = Tuple[float, float]

PixelVector = Tuple[float, float]


def rotate(v: SkyPosition, pivot: SkyPosition, angle: AngleLike, wcs: WCS) -> SkyPosition:
    """Rotate a point around a pivot.

    Both the pivot and the point to rotate are assumed to be given as right ascension
    and declination. The rotation angle is taken to be the angle on the sky. If the
    pixel scales for right ascebsion abd declination differ, the angle seen on the
    finder chart will differ from this angle.

    Parameters
    ----------
    v: tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
        The point to rotate.
    pivot: tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
        The point around which the point ``v`` is rotated.
    angle: Quantity ["angle"]
        The angle of rotation. A positive angle corresponds to an anti-clockwise
        rotation.
    wcs: `~astropy.units.WCS`
        WCS object.

    Returns
    -------
    tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
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

    # Perform the rotation
    rotation = np.array([[cos(angle_px), -sin(angle_px)], [sin(angle_px), cos(angle_px)]])
    diff_vector_px = np.array(v_px) - np.array(pivot_px)
    rotated_v_px = np.array(pivot_px) + rotation @ diff_vector_px

    # Convert back from pixels to world coordinates
    return pixel_to_sky_position(rotated_v_px, wcs)


def translate(v: SkyPosition, displacement: SkyVector, wcs: WCS) -> SkyPosition:
    """Move a point on the by a displacement vector.

    The point to move is assumed to be given as right ascension and declination. The
    displacement vector is assumed to be a vector on the sky, with real angles, as
    described for the `~sky_vector_to_pixel` method.

    Parameters
    ----------
    v: tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
        Point to move.
    displacement: tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
        Vector by which to move the point.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
        The point after moving.
    """
    # Convert from sky coordinates to real angles
    displacement_delta_ra = displacement[0]
    displacement_delta_dec = displacement[1] / cos(v[1].to_value(u.rad))

    # Perform the translation
    w = (v[0] + displacement_delta_ra, v[1] + displacement_delta_dec)
    return w


def pixel_scales(wcs: WCS):
    """
    For a WCS returns pixel scales along each axis of the image pixel at
    the ``CRPIX`` location once it is projected onto the
    "plane of intermediate world coordinates" as defined in
    `Greisen & Calabretta 2002, A&A, 395, 1061 <https://ui.adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_.

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

    Taken from APLpy.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        A world coordinate system object.

    Returns
    -------
    scale : `~numpy.ndarray`
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
    return (Angle(scale_factors[0], wcs.world_axis_units[0]),
            Angle(scale_factors[1], wcs.world_axis_units[1]))


def sky_position_to_pixel(position: SkyPosition, wcs: WCS) -> PixelPosition:
    """Convert a sky position to the corresponding pixel coordinates.

    Parameters
    ----------
    position: tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
        Sky position as right ascension and declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    tuple of (float, float)
        The pixel coordinates.
    """
    print(position)
    sky_coords = SkyCoord(ra=position[0], dec=position[1])
    return tuple(skycoord_to_pixel(coords=sky_coords, wcs=wcs))


def pixel_to_sky_position(position_px: PixelPosition, wcs: WCS):
    """Convert pixel coordinates to the corresponding sky position.

    Parameters
    ----------
    position_px: tuple of (float, float)
        Pixel coordinates.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
        The sky position, as right asension and declination.
    """
    sky_coords = pixel_to_skycoord(position_px[0], position_px[1], wcs)
    return sky_coords.ra, sky_coords.dec


def sky_vector_to_pixel(vector: SkyVector, wcs: WCS) -> PixelVector:
    """Convert a vector on the sky to the corresponding vector in pixels.

    The coordinates of the given vector are assumed to be real angles (rather than
    coordinate differences). For right ascensions the real angle on the sky is the same
    as the coordinate difference, but for all declinations other than 0 the two differ.

    Parameters
    ----------
    vector: tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
        Vector on the sky. The first coordinate is a real angle in the direction of the
        right ascension. The second coordinate is a real angle in the directio of the
        declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    tuple of (float, float)
        The vector in pixels.
    """
    ra_pixel_scale, _ = pixel_scales(wcs)
    vector_px_x = float(vector[0] / ra_pixel_scale)
    vector_px_y = float(vector[1] / ra_pixel_scale)
    return vector_px_x, vector_px_y

