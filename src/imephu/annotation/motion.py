from datetime import timezone
from typing import List, Tuple, Union

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import (
    ArrowAnnotation,
    CircleAnnotation,
    GroupAnnotation,
    LinePathAnnotation,
    TextAnnotation,
)
from imephu.geometry import (
    pixel_scales,
    pixel_to_sky_position,
    sky_position_to_pixel,
    translate,
)
from imephu.utils import Ephemeris

_MAX_SMALL_MOTION_PIXELS = 18.7


def motion_annotation(ephemerides: List[Ephemeris], wcs: WCS) -> GroupAnnotation:
    """Return an annotation showing the motion of an object on the sky.

    If the motion is sufficiently large, it is shown as an arrow from the initial to the
    final position. Labels for the initial and final epoch are added.

    Otherwise an arrow head pointing in the direction of motion and labels for the
    initial and final epoch are shown. A circle is added around the arrowhead for
    highlighting purposes.

    A sufficiently large motion is a motion where the angular distance between the
    initial and final position is larger than 2 / 75 of the finder chart size.

    Parameters
    ----------
    ephemerides: list of `~imephu.utils.Ephemeris`
        The ephemerides for the motion.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.GroupAnnotation`
        The annotation describing the motion.
    """
    sorted_ephemerides = sorted(ephemerides[:], key=lambda e: e.epoch)
    initial_position = sorted_ephemerides[0].position
    final_position = sorted_ephemerides[-1].position
    distance_travelled = initial_position.separation(final_position)
    if distance_travelled == 0 * u.arcsec:
        raise ValueError("The initial and final position must not be identical.")

    color = "blue"
    if distance_travelled > _MAX_SMALL_MOTION_PIXELS * _angle_per_pixel(wcs):
        return _large_motion_annotation(sorted_ephemerides, wcs, color)
    else:
        return _small_motion_annotation(sorted_ephemerides, wcs, color)


def _large_motion_annotation(
    sorted_ephemerides: List[Ephemeris],
    wcs: WCS,
    color: Union[str, Tuple[float, float, float]],
) -> GroupAnnotation:
    # path showing the motion
    path_annotation = LinePathAnnotation(
        vertices=[e.position for e in sorted_ephemerides],
        wcs=wcs,
        closed=False,
        edgecolor=color,
    )

    # arrow head
    angle_per_pixel = _angle_per_pixel(wcs)
    arrow_annotation = ArrowAnnotation(
        start=sorted_ephemerides[-2].position,
        end=sorted_ephemerides[-1].position,
        wcs=wcs,
        head_width=5.5 * angle_per_pixel,
        head_height=8.4 * angle_per_pixel,
        head_filled=False,
        tail_shown=False,
        color=color,
    )

    # epoch labels
    initial_epoch_label = _large_motion_epoch_label(
        sorted_ephemerides=sorted_ephemerides,
        for_initial_epoch=True,
        color=color,
        wcs=wcs,
    )
    final_epoch_label = _large_motion_epoch_label(
        sorted_ephemerides=sorted_ephemerides,
        for_initial_epoch=False,
        color=color,
        wcs=wcs,
    )

    return GroupAnnotation(
        [path_annotation, arrow_annotation, initial_epoch_label, final_epoch_label]
    )


def _small_motion_annotation(
    sorted_ephemerides: List[Ephemeris],
    wcs: WCS,
    color: Union[str, Tuple[float, float, float]],
) -> GroupAnnotation:
    # mid position on the line between initial and final position on the finder chart,
    # which may be slightly different from the actual mid position on a great circle
    # between the points on the sky
    mid_position = _mid_point_on_chart(sorted_ephemerides, wcs)
    mid_position_px = sky_position_to_pixel(mid_position, wcs)

    # calculate an arrow start and end so that the arrow head will be placed at the
    # right direction and will point in the right direction
    initial_position_px = sky_position_to_pixel(sorted_ephemerides[0].position, wcs)
    final_position_px = sky_position_to_pixel(sorted_ephemerides[-1].position, wcs)
    direction = final_position_px - initial_position_px
    direction_unit = direction / np.linalg.norm(direction)
    arrow_start_px = mid_position_px - (_MAX_SMALL_MOTION_PIXELS / 2) * direction_unit
    arrow_start = pixel_to_sky_position(arrow_start_px, wcs)
    arrow_end_px = mid_position_px + (_MAX_SMALL_MOTION_PIXELS / 2) * direction_unit
    arrow_end = pixel_to_sky_position(arrow_end_px, wcs)

    # arrow head
    angle_per_pixel = _angle_per_pixel(wcs)
    arrow_annotation = ArrowAnnotation(
        start=arrow_start,
        end=arrow_end,
        wcs=wcs,
        head_width=15 * angle_per_pixel,
        head_height=15 * angle_per_pixel,
        head_filled=False,
        tail_shown=False,
        color=color,
    )

    # circle for highlighting the arrow
    circle_annotation = CircleAnnotation(
        center=mid_position,
        radius=(_MAX_SMALL_MOTION_PIXELS / 2 + 5) * angle_per_pixel,
        wcs=wcs,
        edgecolor=color,
    )

    # epoch labels
    initial_epoch_label = _small_motion_epoch_label(
        sorted_ephemerides=sorted_ephemerides,
        for_initial_epoch=True,
        color=color,
        wcs=wcs,
    )
    final_epoch_label = _small_motion_epoch_label(
        sorted_ephemerides=sorted_ephemerides,
        for_initial_epoch=False,
        color=color,
        wcs=wcs,
    )

    return GroupAnnotation(
        [arrow_annotation, circle_annotation, initial_epoch_label, final_epoch_label]
    )


def _large_motion_epoch_label(
    sorted_ephemerides: List[Ephemeris],
    for_initial_epoch: bool,
    color: Union[str, Tuple[float, float, float]],
    wcs: WCS,
) -> TextAnnotation:
    # For the initial epoch the label should be put slightly to the south if the object
    # is moving north.
    angle_per_pixel = _angle_per_pixel(wcs)
    initial_dec_displacement = 12 * angle_per_pixel
    if sorted_ephemerides[1].position.dec >= sorted_ephemerides[0].position.dec:
        initial_dec_displacement *= -1

    # The final epoch label should be displaced in the opposite direction as compared to
    # the initial epoch one.
    final_dec_displacement = -initial_dec_displacement

    dec_displacement = (
        initial_dec_displacement if for_initial_epoch else final_dec_displacement
    )
    dec_displacement_arcsec = dec_displacement.to_value(u.arcsec)
    if for_initial_epoch:
        ephemeris = sorted_ephemerides[0]
    else:
        ephemeris = sorted_ephemerides[-1]
    position = translate(ephemeris.position, (0, dec_displacement_arcsec) * u.arcsec)
    text = ephemeris.epoch.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M UT")

    return TextAnnotation(position=position, s=text, wcs=wcs, color=color, fontsize=8)


def _small_motion_epoch_label(
    sorted_ephemerides: List[Ephemeris],
    for_initial_epoch: bool,
    color: Union[str, Tuple[float, float, float]],
    wcs: WCS,
) -> TextAnnotation:
    # For the initial epoch the label should be put slightly to the south if the object
    # is moving north.
    angle_per_pixel = _angle_per_pixel(wcs)
    initial_dec_displacement = (_MAX_SMALL_MOTION_PIXELS / 2 + 13) * angle_per_pixel
    if sorted_ephemerides[1].position.dec >= sorted_ephemerides[0].position.dec:
        initial_dec_displacement *= -1

    # The final epoch label should be displaced in the opposite direction as compared to
    # the initial epoch one.
    final_dec_displacement = -initial_dec_displacement

    dec_displacement = (
        initial_dec_displacement if for_initial_epoch else final_dec_displacement
    )
    dec_displacement_arcsec = dec_displacement.to_value(u.arcsec)
    if for_initial_epoch:
        ephemeris = sorted_ephemerides[0]
    else:
        ephemeris = sorted_ephemerides[-1]
    mid_point = _mid_point_on_chart(sorted_ephemerides, wcs)
    position = translate(mid_point, (0, dec_displacement_arcsec) * u.arcsec)
    text = ephemeris.epoch.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M UT")

    return TextAnnotation(position=position, s=text, wcs=wcs, color=color, fontsize=8)


def _angle_per_pixel(wcs: WCS) -> Angle:
    """Return the angle per pixel for a WCS object.

    The angle per pixel is averaged over right ascension and declination.
    """
    sx, sy = pixel_scales(wcs)
    return (sx + sy) / 2


def _mid_point_on_chart(sorted_ephemerides: List[Ephemeris], wcs: WCS) -> SkyCoord:
    """Return the midpoint between the initial and final position.

    The midpoint is the point on the sky which on the finder chart corresponds to the
    midpoint of the line between the initial and final position. This in general may be
    different from the midpoint on the great circle through the initial and final
    positions on the sky.
    """
    initial_position_px = sky_position_to_pixel(sorted_ephemerides[0].position, wcs)
    final_position_px = sky_position_to_pixel(sorted_ephemerides[-1].position, wcs)
    mid_position_px = (initial_position_px + final_position_px) / 2
    return pixel_to_sky_position(mid_position_px, wcs)
