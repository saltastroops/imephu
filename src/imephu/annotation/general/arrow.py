from typing import Any, Tuple, Union

import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from numpy import typing as npt

from imephu.annotation.general import GroupAnnotation, LinePathAnnotation
from imephu.geometry import pixel_scales, pixel_to_sky_position, sky_position_to_pixel


class ArrowAnnotation(GroupAnnotation):
    """An annotation for plotting an arrow from a start to an end position.

    The start of the arrow's tail and the tip of the arrow's head are given as sky
    positions (in right ascension and declination). The width and height of the arrow
    head are given as an angle on the sky, and this angle is converted to pixels using
    the pixel scale in right ascension direction.

     By default the arrow is not filled, but you can change this behavior with the
     ``head_filled`` parameter.

    Parameters
    ----------
    start: `~astropy.coordinates.SkyCoord`
        Position of the arrow tail's start.
    end: `~astropy.coordinates.SkyCoord`
        Position of the arrow head's tip.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    head_width: `~astropy.coordinates.Angle`
        Width of the arrow head.
    head_height: `~astropy.coordinates.Angle`
        Height of the arrow head.
    head_filled: `bool`, default: False
        Whether the arrow head is filled.
    tail_shown: `bool`, default: True
        Whether the tail should be shown. If ``False``, only the head is displayed.
    color: color, default: "black"
        The arrow color.
    **kwargs: dict, optional
        Additional keyword arguments, which will be passed to the
        `~imephu.annotation.general.LinePathAnnotation` constructor for both tail and
        head.
    """

    def __init__(
        self,
        start: SkyCoord,
        end: SkyCoord,
        wcs: WCS,
        head_width: Angle,
        head_height: Angle,
        head_filled: bool = False,
        tail_shown: bool = True,
        color: Union[str, Tuple[float, float, float]] = "black",
        **kwargs: Any,
    ):
        super().__init__([])
        if tail_shown:
            super().add_item(ArrowAnnotation._tail(start, end, wcs, color, **kwargs))
        super().add_item(
            ArrowAnnotation._head(
                start, end, wcs, head_width, head_height, head_filled, color, **kwargs
            )
        )

    @staticmethod
    def _tail(
        start: SkyCoord,
        end: SkyCoord,
        wcs: WCS,
        color: Union[str, Tuple[float, float, float]],
        **kwargs: Any,
    ) -> LinePathAnnotation:
        return LinePathAnnotation(
            vertices=[start, end], wcs=wcs, closed=False, edgecolor=color, **kwargs
        )

    @staticmethod
    def _head(
        start: SkyCoord,
        end: SkyCoord,
        wcs: WCS,
        head_width: Angle,
        head_height: Angle,
        head_filled: bool,
        color: Union[str, Tuple[float, float, float]],
        **kwargs: Any,
    ) -> LinePathAnnotation:
        # Let S, be the start of the arrow's tail, T be the tip of the arrow and A, B
        # the other two vertices of the head. Then AA' is bisected by the line ST, and
        # AB and ST are perpendicular to each other.
        # Let's get the unit vector n_TS from head to tail in pixels.
        start_px = sky_position_to_pixel(start, wcs)
        end_px = sky_position_to_pixel(end, wcs)
        end_to_start_px = start_px - end_px  # noqa
        n_TS_px = end_to_start_px / np.linalg.norm(end_to_start_px)

        # We can use the vector to get the point of intersection of AA' and ST.
        sx, _ = pixel_scales(wcs)
        head_height_px = float(head_height / sx)
        point_of_intersection = end_px + head_height_px * n_TS_px  # noqa

        # We can also use it to get a unit vector n_AB along AB.
        n_AB_px: npt.NDArray[np.float_] = np.array([n_TS_px[1], -n_TS_px[0]])

        # And hence we can calculate A and B.
        head_width_px = float(head_width / sx) / 2
        v_A_px = point_of_intersection - head_width_px * n_AB_px  # noqa
        v_B_px = point_of_intersection + head_width_px * n_AB_px  # noqa
        v_A = pixel_to_sky_position(v_A_px, wcs)
        v_B = pixel_to_sky_position(v_B_px, wcs)

        return LinePathAnnotation(
            vertices=[v_A, end, v_B],
            wcs=wcs,
            closed=head_filled,
            edgecolor=color,
            facecolor=color if head_filled else "none",
            **kwargs,
        )
