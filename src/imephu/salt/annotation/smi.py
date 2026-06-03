from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterable, List, Protocol, Tuple

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu import finder_chart
from imephu.annotation.general import (
    CircleAnnotation,
    GroupAnnotation,
    LinePathAnnotation,
)
from imephu.geometry import translate

_OUTLINE_COLOR = "red"


def smi_factory(smi_barcode: str) -> SMI:
    match smi_barcode:
        case "PF0200N001":
            return PF0200N001()
        case "PF0300N001":
            return PF0300N001()
        case _:
            raise ValueError(f"Unknown slit mask IFU barcode: {smi_barcode}")


class SMI(ABC):
    def bundle(
        self, center: SkyCoord, wcs: WCS, include_fibers: bool
    ) -> GroupAnnotation:
        return GroupAnnotation(
            [
                self._left_sky_bundle(center, wcs, include_fibers),
                self._target_bundle(center, wcs, include_fibers),
                self._right_sky_bundle(center, wcs, include_fibers),
            ]
        )

    @abstractmethod
    def fibers(self) -> List[Tuple[Angle, Angle]]:
        """Get the fibers for the slit mask IFU."""
        pass

    @abstractmethod
    def fiber_radius(self) -> Angle:
        """Get the fiber radius for the slit mask IFU."""
        pass

    @abstractmethod
    def target_bundle_center_row_extension(self) -> Angle:
        """
        Get the amount by which the center row should be extended in x direction from
        the outermost fiber center for outline purposes.
        """
        pass

    @staticmethod
    def _target_bundle_x_cutoffs() -> Tuple[Angle, Angle]:
        """
        Return a lower and upper cutoff for the target bundle.

        The method returns a tuple x1, x2 so that all left sky bundle fibers have an x
        value less than x1, all target bundle fibers have an x value between x1 and x2,
        and all right sky bundle fibers have an x value greater than x2.

        Returns
        -------
        pair of Angle
            A lower and upper cutoff for the target bundle.
        """
        return -25 * u.arcsec, 25 * u.arcsec

    def _left_sky_bundle_fibers(self) -> List[Tuple[Angle, Angle]]:
        return [c for c in self.fibers() if c[0] < self._target_bundle_x_cutoffs()[0]]

    def _target_bundle_fibers(self) -> List[Tuple[Angle, Angle]]:
        x_range = self._target_bundle_x_cutoffs()
        return [c for c in self.fibers() if x_range[0] < c[0] < x_range[1]]

    def _right_sky_bundle_fibers(self) -> List[Tuple[Angle, Angle]]:
        return [c for c in self.fibers() if c[0] > self._target_bundle_x_cutoffs()[1]]

    @staticmethod
    def _fiber_outlines(
        fibers: Iterable[Tuple[Angle, Angle]], radius: Angle, center: SkyCoord, wcs: WCS
    ) -> GroupAnnotation:
        annotation = GroupAnnotation([])
        for x, y in fibers:
            offset_x = x.to_value(u.arcsec)
            offset_y = y.to_value(u.arcsec)

            fiber_annotation = CircleAnnotation(
                center=center,
                radius=radius,
                wcs=wcs,
                edgecolor="red",
            ).translate((offset_x, offset_y) * u.arcsec)
            annotation.add_item(fiber_annotation)

        return annotation

    def _sky_bundle(
        self,
        fibers: List[Tuple[Angle, Angle]],
        center: SkyCoord,
        wcs: WCS,
        include_fibers: bool,
    ) -> GroupAnnotation:
        annotation = GroupAnnotation([])
        if include_fibers:
            annotation.add_item(
                self._fiber_outlines(fibers, self.fiber_radius(), center, wcs)
            )
        annotation.add_item(self._sky_bundle_outline(fibers, center, wcs))

        return annotation

    def _sky_bundle_outline(
        self, fibers: List[Tuple[Angle, Angle]], center: SkyCoord, wcs: WCS
    ) -> LinePathAnnotation:
        # Approximate the outline by a polygon consisting of horizontal and vertical lines.
        r = self.fiber_radius()
        outline_x_min = min(f[0] for f in fibers) - r
        outline_x_max = max(f[0] for f in fibers) + r
        outline_y_min = min(f[1] for f in fibers) - r
        outline_y_max = max(f[1] for f in fibers) + r

        vertices = [
            (outline_x_min, outline_y_min),
            (outline_x_max, outline_y_min),
            (outline_x_max, outline_y_max),
            (outline_x_min, outline_y_max),
        ]

        # So far we have assumed the bundle center to be at (0, 0), so we have to correct
        # this.
        vertex_positions = [translate(center, v) for v in vertices]

        line = LinePathAnnotation(
            vertex_positions, wcs=wcs, closed=True, edgecolor=_OUTLINE_COLOR
        )
        return line

    def _left_sky_bundle(
        self, center: SkyCoord, wcs: WCS, include_fibers: bool
    ) -> GroupAnnotation:
        return self._sky_bundle(
            self._left_sky_bundle_fibers(), center, wcs, include_fibers
        )

    def _right_sky_bundle(
        self, center: SkyCoord, wcs: WCS, include_fibers: bool
    ) -> GroupAnnotation:
        return self._sky_bundle(
            self._right_sky_bundle_fibers(), center, wcs, include_fibers
        )

    def _target_bundle(
        self, center: SkyCoord, wcs: WCS, include_fibers: bool
    ) -> GroupAnnotation:
        annotation = GroupAnnotation([])
        if include_fibers:
            target_bundle_fibers = self._target_bundle_fibers()
            annotation.add_item(
                PF0200N001._fiber_outlines(
                    target_bundle_fibers, self.fiber_radius(), center, wcs
                )
            )
        annotation.add_item(self._target_bundle_outline(center, wcs))

        return annotation

    @staticmethod
    def _x_range(fibers: list[Tuple[Angle, Angle]]) -> Tuple[Angle, Angle]:
        return min(f[0] for f in fibers), max(f[0] for f in fibers)

    @staticmethod
    def _y_range(fibers: list[Tuple[Angle, Angle]]) -> Tuple[Angle, Angle]:
        return min(f[1] for f in fibers), max(f[1] for f in fibers)

    def _target_bundle_outline(self, center: SkyCoord, wcs: WCS) -> LinePathAnnotation:
        # Approximate the outline by a vertically and horizontally symmetric hexagon.
        r = self.fiber_radius()
        target_bundle_fibers = self._target_bundle_fibers()

        # Get the bottom row, center row and top row.
        x_min, x_max = self._x_range(target_bundle_fibers)
        y_min, y_max = self._y_range(target_bundle_fibers)

        bottom_row = [f for f in target_bundle_fibers if f[1] < y_min + 0.25 * r]
        some_center_row_fiber = [
            f
            for f in target_bundle_fibers
            if x_min - 0.001 * r < f[0] < x_min + 0.001 * r
        ][0]
        center_row = [
            f
            for f in target_bundle_fibers
            if some_center_row_fiber[1] - 0.25 * r
            < f[1]
            < some_center_row_fiber[1] + 0.25 * r
        ]
        top_row = [f for f in target_bundle_fibers if f[1] > y_max - 0.25 * r]

        top_row_outline_x_min = min(c[0] for c in top_row) - r - 0.25 * r
        top_row_outline_x_max = -top_row_outline_x_min
        top_row_outline_y_max = max(c[1] for c in top_row) + r
        bottom_row_outline_x_min = top_row_outline_x_min
        bottom_row_outline_x_max = top_row_outline_x_max
        bottom_row_outline_y_min = min(c[1] for c in bottom_row) - r
        center_row_outline_x_min = (
            min(c[0] for c in center_row) - self.target_bundle_center_row_extension()
        )
        center_row_outline_x_max = -center_row_outline_x_min
        center_row_outline_y = (top_row_outline_y_max + bottom_row_outline_y_min) / 2

        vertices = [
            (top_row_outline_x_min, top_row_outline_y_max),
            (top_row_outline_x_max, top_row_outline_y_max),
            (center_row_outline_x_max, center_row_outline_y),
            (bottom_row_outline_x_max, bottom_row_outline_y_min),
            (bottom_row_outline_x_min, bottom_row_outline_y_min),
            (center_row_outline_x_min, center_row_outline_y),
            (top_row_outline_x_min, top_row_outline_y_max),
        ]

        # So far we have assumed the bundle center to be at (0, 0), so we have to correct
        # this.
        vertex_positions = [translate(center, v) for v in vertices]

        line = LinePathAnnotation(
            vertex_positions, wcs=wcs, closed=True, edgecolor=_OUTLINE_COLOR
        )
        return line


class PF0200N001(SMI):
    def fiber_radius(self) -> Angle:
        return 0.44 * u.arcsec

    def target_bundle_center_row_extension(self) -> Angle:
        return 1.25 * self.fiber_radius()

    def fibers(self) -> List[Tuple[Angle, Angle]]:
        return [
            (-4.345 * u.arcsec, -0.04 * u.arcsec),
            (-3.827 * u.arcsec, 0.841 * u.arcsec),
            (-3.248 * u.arcsec, -3.854 * u.arcsec),
            (-96.903 * u.arcsec, 0.164 * u.arcsec),
            (-1.602 * u.arcsec, -2.836 * u.arcsec),
            (-0.527 * u.arcsec, 2.739 * u.arcsec),
            (-1.619 * u.arcsec, -4.792 * u.arcsec),
            (-3.77 * u.arcsec, -0.969 * u.arcsec),
            (-1.044 * u.arcsec, 3.726 * u.arcsec),
            (-3.279 * u.arcsec, 0.004 * u.arcsec),
            (-1.071 * u.arcsec, -3.854 * u.arcsec),
            (-2.111 * u.arcsec, 3.726 * u.arcsec),
            (-2.142 * u.arcsec, -1.907 * u.arcsec),
            (-1.637 * u.arcsec, -0.969 * u.arcsec),
            (-3.73 * u.arcsec, -2.885 * u.arcsec),
            (-0.54 * u.arcsec, 0.889 * u.arcsec),
            (-3.257 * u.arcsec, -1.903 * u.arcsec),
            (-2.124 * u.arcsec, 1.827 * u.arcsec),
            (-2.146 * u.arcsec, -3.854 * u.arcsec),
            (0.0 * u.arcsec, 0.0 * u.arcsec),
            (-0.496 * u.arcsec, -2.836 * u.arcsec),
            (-96.08 * u.arcsec, 0.181 * u.arcsec),
            (-1.637 * u.arcsec, 2.73 * u.arcsec),
            (-1.075 * u.arcsec, 0.004 * u.arcsec),
            (-3.796 * u.arcsec, 2.739 * u.arcsec),
            (-1.695 * u.arcsec, 0.885 * u.arcsec),
            (-0.504 * u.arcsec, -4.783 * u.arcsec),
            (-1.084 * u.arcsec, -1.912 * u.arcsec),
            (-2.721 * u.arcsec, -2.845 * u.arcsec),
            (-93.973 * u.arcsec, 0.181 * u.arcsec),
            (-95.044 * u.arcsec, 0.186 * u.arcsec),
            (-0.473 * u.arcsec, 4.801 * u.arcsec),
            (-1.655 * u.arcsec, 4.73 * u.arcsec),
            (-4.35 * u.arcsec, -1.85 * u.arcsec),
            (-2.739 * u.arcsec, 2.774 * u.arcsec),
            (-1.062 * u.arcsec, 1.827 * u.arcsec),
            (-2.133 * u.arcsec, 0.004 * u.arcsec),
            (-4.385 * u.arcsec, 1.788 * u.arcsec),
            (-2.752 * u.arcsec, 0.894 * u.arcsec),
            (-97.681 * u.arcsec, 0.279 * u.arcsec),
            (-2.721 * u.arcsec, -0.965 * u.arcsec),
            (-3.332 * u.arcsec, 1.783 * u.arcsec),
            (-0.509 * u.arcsec, -0.965 * u.arcsec),
            (-3.265 * u.arcsec, 3.739 * u.arcsec),
            (-5.416 * u.arcsec, -0.044 * u.arcsec),
            (-2.009 * u.arcsec, 7.597 * u.arcsec),
            (-4.358 * u.arcsec, 3.743 * u.arcsec),
            (-3.252 * u.arcsec, 5.544 * u.arcsec),
            (-2.708 * u.arcsec, 6.704 * u.arcsec),
            (-1.633 * u.arcsec, -6.69 * u.arcsec),
            (-5.425 * u.arcsec, -1.956 * u.arcsec),
            (-89.252 * u.arcsec, -0.739 * u.arcsec),
            (-5.407 * u.arcsec, 1.792 * u.arcsec),
            (-0.903 * u.arcsec, 7.513 * u.arcsec),
            (-4.858 * u.arcsec, 2.748 * u.arcsec),
            (-2.686 * u.arcsec, -6.646 * u.arcsec),
            (-1.031 * u.arcsec, -5.761 * u.arcsec),
            (-2.664 * u.arcsec, -4.788 * u.arcsec),
            (-3.739 * u.arcsec, -4.788 * u.arcsec),
            (-1.447 * u.arcsec, 8.628 * u.arcsec),
            (-4.31 * u.arcsec, -3.858 * u.arcsec),
            (-2.102 * u.arcsec, -5.761 * u.arcsec),
            (-1.049 * u.arcsec, -7.58 * u.arcsec),
            (-94.345 * u.arcsec, -0.735 * u.arcsec),
            (-4.876 * u.arcsec, 0.841 * u.arcsec),
            (-0.522 * u.arcsec, -6.69 * u.arcsec),
            (-5.987 * u.arcsec, 0.845 * u.arcsec),
            (-1.606 * u.arcsec, 6.695 * u.arcsec),
            (0.022 * u.arcsec, -7.58 * u.arcsec),
            (-3.81 * u.arcsec, 4.646 * u.arcsec),
            (-6.482 * u.arcsec, -0.04 * u.arcsec),
            (-5.956 * u.arcsec, -1.013 * u.arcsec),
            (0.031 * u.arcsec, -5.677 * u.arcsec),
            (-2.142 * u.arcsec, 5.549 * u.arcsec),
            (-4.841 * u.arcsec, -2.934 * u.arcsec),
            (-0.681 * u.arcsec, 6.677 * u.arcsec),
            (-91.425 * u.arcsec, -0.712 * u.arcsec),
            (-1.624 * u.arcsec, -8.504 * u.arcsec),
            (-2.73 * u.arcsec, 4.659 * u.arcsec),
            (-3.204 * u.arcsec, -5.765 * u.arcsec),
            (-2.106 * u.arcsec, -7.562 * u.arcsec),
            (-4.841 * u.arcsec, -1.018 * u.arcsec),
            (-1.035 * u.arcsec, 5.673 * u.arcsec),
            (-4.881 * u.arcsec, 4.642 * u.arcsec),
            (-4.319 * u.arcsec, -5.765 * u.arcsec),
            (-4.31 * u.arcsec, 5.527 * u.arcsec),
            (-3.677 * u.arcsec, -6.642 * u.arcsec),
            (-4.805 * u.arcsec, -6.726 * u.arcsec),
            (-8.075 * u.arcsec, -1.009 * u.arcsec),
            (-5.991 * u.arcsec, 4.642 * u.arcsec),
            (-90.801 * u.arcsec, 0.177 * u.arcsec),
            (-4.854 * u.arcsec, -4.801 * u.arcsec),
            (-5.385 * u.arcsec, -3.867 * u.arcsec),
            (-6.531 * u.arcsec, 1.779 * u.arcsec),
            (-5.912 * u.arcsec, -2.805 * u.arcsec),
            (-6.491 * u.arcsec, 3.739 * u.arcsec),
            (-7.071 * u.arcsec, -1.018 * u.arcsec),
            (-5.907 * u.arcsec, -4.801 * u.arcsec),
            (-7.084 * u.arcsec, 2.743 * u.arcsec),
            (-8.106 * u.arcsec, 0.85 * u.arcsec),
            (-4.881 * u.arcsec, 6.58 * u.arcsec),
            (-8.69 * u.arcsec, -0.058 * u.arcsec),
            (-2.686 * u.arcsec, -8.478 * u.arcsec),
            (-88.681 * u.arcsec, 0.173 * u.arcsec),
            (-6.434 * u.arcsec, -3.858 * u.arcsec),
            (-5.327 * u.arcsec, -5.765 * u.arcsec),
            (-5.429 * u.arcsec, 3.708 * u.arcsec),
            (-3.124 * u.arcsec, 7.491 * u.arcsec),
            (-7.013 * u.arcsec, -2.885 * u.arcsec),
            (-7.646 * u.arcsec, 1.774 * u.arcsec),
            (-7.597 * u.arcsec, -1.96 * u.arcsec),
            (-3.204 * u.arcsec, -7.504 * u.arcsec),
            (-4.199 * u.arcsec, 7.544 * u.arcsec),
            (-3.819 * u.arcsec, 6.668 * u.arcsec),
            (-6.5 * u.arcsec, -1.96 * u.arcsec),
            (-7.049 * u.arcsec, 0.845 * u.arcsec),
            (-89.721 * u.arcsec, 0.181 * u.arcsec),
            (-5.978 * u.arcsec, 2.761 * u.arcsec),
            (-3.659 * u.arcsec, 8.456 * u.arcsec),
            (-5.367 * u.arcsec, 5.562 * u.arcsec),
            (-3.695 * u.arcsec, -8.5 * u.arcsec),
            (-2.531 * u.arcsec, 8.593 * u.arcsec),
            (-7.646 * u.arcsec, -0.053 * u.arcsec),
            (-4.274 * u.arcsec, -7.58 * u.arcsec),
            (-9.668 * u.arcsec, -1.973 * u.arcsec),
            (-5.341 * u.arcsec, -7.584 * u.arcsec),
            (-6.429 * u.arcsec, -7.5 * u.arcsec),
            (-6.513 * u.arcsec, 5.558 * u.arcsec),
            (-8.606 * u.arcsec, -3.903 * u.arcsec),
            (-7.527 * u.arcsec, -3.858 * u.arcsec),
            (-7.597 * u.arcsec, 3.739 * u.arcsec),
            (-87.611 * u.arcsec, 0.177 * u.arcsec),
            (-8.155 * u.arcsec, 2.748 * u.arcsec),
            (-9.19 * u.arcsec, -1.031 * u.arcsec),
            (-6.0 * u.arcsec, 6.54 * u.arcsec),
            (-7.013 * u.arcsec, -6.704 * u.arcsec),
            (-4.681 * u.arcsec, 8.54 * u.arcsec),
            (-7.022 * u.arcsec, -4.796 * u.arcsec),
            (-5.336 * u.arcsec, 7.469 * u.arcsec),
            (-9.226 * u.arcsec, 2.735 * u.arcsec),
            (-8.088 * u.arcsec, 4.619 * u.arcsec),
            (-5.85 * u.arcsec, -6.726 * u.arcsec),
            (-9.757 * u.arcsec, 1.783 * u.arcsec),
            (-7.035 * u.arcsec, 4.642 * u.arcsec),
            (-9.757 * u.arcsec, -0.053 * u.arcsec),
            (-4.814 * u.arcsec, -8.513 * u.arcsec),
            (-7.478 * u.arcsec, -5.761 * u.arcsec),
            (-5.73 * u.arcsec, 8.544 * u.arcsec),
            (-7.049 * u.arcsec, 6.553 * u.arcsec),
            (-7.611 * u.arcsec, 5.655 * u.arcsec),
            (-8.022 * u.arcsec, -4.792 * u.arcsec),
            (-8.726 * u.arcsec, 3.717 * u.arcsec),
            (-8.726 * u.arcsec, 1.792 * u.arcsec),
            (-10.765 * u.arcsec, -0.04 * u.arcsec),
            (-88.159 * u.arcsec, -0.752 * u.arcsec),
            (-9.071 * u.arcsec, -2.92 * u.arcsec),
            (-7.969 * u.arcsec, -2.796 * u.arcsec),
            (-10.195 * u.arcsec, -1.022 * u.arcsec),
            (-9.159 * u.arcsec, 0.845 * u.arcsec),
            (-8.646 * u.arcsec, -1.96 * u.arcsec),
            (-6.442 * u.arcsec, -5.77 * u.arcsec),
            (-5.854 * u.arcsec, -8.509 * u.arcsec),
            (-6.487 * u.arcsec, 7.588 * u.arcsec),
            (-10.23 * u.arcsec, 0.85 * u.arcsec),
            (6.513 * u.arcsec, 5.633 * u.arcsec),
            (10.288 * u.arcsec, -0.889 * u.arcsec),
            (90.004 * u.arcsec, -1.159 * u.arcsec),
            (8.142 * u.arcsec, -2.881 * u.arcsec),
            (9.243 * u.arcsec, -0.81 * u.arcsec),
            (7.611 * u.arcsec, 3.889 * u.arcsec),
            (7.0 * u.arcsec, -6.615 * u.arcsec),
            (7.633 * u.arcsec, -3.81 * u.arcsec),
            (5.42 * u.arcsec, -7.535 * u.arcsec),
            (8.146 * u.arcsec, -4.748 * u.arcsec),
            (9.212 * u.arcsec, -2.881 * u.arcsec),
            (8.642 * u.arcsec, -1.867 * u.arcsec),
            (8.549 * u.arcsec, 1.823 * u.arcsec),
            (4.916 * u.arcsec, -8.465 * u.arcsec),
            (5.942 * u.arcsec, -6.606 * u.arcsec),
            (8.066 * u.arcsec, 2.81 * u.arcsec),
            (6.527 * u.arcsec, -7.469 * u.arcsec),
            (91.053 * u.arcsec, -1.212 * u.arcsec),
            (8.726 * u.arcsec, -3.743 * u.arcsec),
            (6.075 * u.arcsec, 8.553 * u.arcsec),
            (7.04 * u.arcsec, -4.743 * u.arcsec),
            (9.173 * u.arcsec, 2.823 * u.arcsec),
            (7.58 * u.arcsec, -5.717 * u.arcsec),
            (5.942 * u.arcsec, 6.58 * u.arcsec),
            (6.513 * u.arcsec, -5.708 * u.arcsec),
            (9.699 * u.arcsec, -1.863 * u.arcsec),
            (7.035 * u.arcsec, 4.735 * u.arcsec),
            (8.177 * u.arcsec, 4.668 * u.arcsec),
            (5.456 * u.arcsec, 7.69 * u.arcsec),
            (90.597 * u.arcsec, -0.155 * u.arcsec),
            (9.177 * u.arcsec, 0.973 * u.arcsec),
            (8.695 * u.arcsec, 3.819 * u.arcsec),
            (10.814 * u.arcsec, 0.088 * u.arcsec),
            (7.566 * u.arcsec, 5.633 * u.arcsec),
            (5.027 * u.arcsec, 8.522 * u.arcsec),
            (9.699 * u.arcsec, 1.881 * u.arcsec),
            (10.283 * u.arcsec, 1.022 * u.arcsec),
            (6.951 * u.arcsec, 6.588 * u.arcsec),
            (3.239 * u.arcsec, -7.584 * u.arcsec),
            (7.08 * u.arcsec, -2.788 * u.arcsec),
            (6.611 * u.arcsec, 3.92 * u.arcsec),
            (5.987 * u.arcsec, -4.739 * u.arcsec),
            (3.867 * u.arcsec, 8.562 * u.arcsec),
            (5.425 * u.arcsec, -5.712 * u.arcsec),
            (8.071 * u.arcsec, 0.938 * u.arcsec),
            (4.296 * u.arcsec, -5.761 * u.arcsec),
            (6.522 * u.arcsec, -3.81 * u.arcsec),
            (7.013 * u.arcsec, 0.934 * u.arcsec),
            (8.655 * u.arcsec, 0.04 * u.arcsec),
            (5.473 * u.arcsec, 5.619 * u.arcsec),
            (6.0 * u.arcsec, -2.836 * u.arcsec),
            (7.588 * u.arcsec, -1.867 * u.arcsec),
            (5.451 * u.arcsec, 3.805 * u.arcsec),
            (6.416 * u.arcsec, 1.827 * u.arcsec),
            (7.058 * u.arcsec, -0.889 * u.arcsec),
            (5.925 * u.arcsec, 2.858 * u.arcsec),
            (4.841 * u.arcsec, 6.575 * u.arcsec),
            (7.487 * u.arcsec, 1.819 * u.arcsec),
            (2.774 * u.arcsec, 8.575 * u.arcsec),
            (4.323 * u.arcsec, 7.575 * u.arcsec),
            (8.173 * u.arcsec, -0.889 * u.arcsec),
            (6.969 * u.arcsec, 2.814 * u.arcsec),
            (6.071 * u.arcsec, 4.85 * u.arcsec),
            (6.522 * u.arcsec, -1.863 * u.arcsec),
            (4.398 * u.arcsec, 5.571 * u.arcsec),
            (3.863 * u.arcsec, -8.513 * u.arcsec),
            (4.929 * u.arcsec, -6.642 * u.arcsec),
            (4.296 * u.arcsec, -7.58 * u.arcsec),
            (7.597 * u.arcsec, 0.044 * u.arcsec),
            (3.279 * u.arcsec, 7.58 * u.arcsec),
            (3.757 * u.arcsec, 6.602 * u.arcsec),
            (2.735 * u.arcsec, -8.473 * u.arcsec),
            (91.637 * u.arcsec, -0.257 * u.arcsec),
            (3.823 * u.arcsec, -6.65 * u.arcsec),
            (92.659 * u.arcsec, -0.252 * u.arcsec),
            (5.416 * u.arcsec, -3.819 * u.arcsec),
            (4.965 * u.arcsec, -4.655 * u.arcsec),
            (93.752 * u.arcsec, -0.257 * u.arcsec),
            (1.659 * u.arcsec, 6.686 * u.arcsec),
            (4.358 * u.arcsec, -3.81 * u.arcsec),
            (0.606 * u.arcsec, -8.509 * u.arcsec),
            (0.677 * u.arcsec, 8.566 * u.arcsec),
            (1.075 * u.arcsec, 7.686 * u.arcsec),
            (0.15 * u.arcsec, 7.611 * u.arcsec),
            (1.633 * u.arcsec, 8.442 * u.arcsec),
            (94.296 * u.arcsec, -1.204 * u.arcsec),
            (4.332 * u.arcsec, 3.765 * u.arcsec),
            (4.85 * u.arcsec, 2.827 * u.arcsec),
            (2.23 * u.arcsec, 5.584 * u.arcsec),
            (2.204 * u.arcsec, -7.575 * u.arcsec),
            (1.146 * u.arcsec, 5.664 * u.arcsec),
            (0.558 * u.arcsec, 6.562 * u.arcsec),
            (0.597 * u.arcsec, -6.646 * u.arcsec),
            (3.894 * u.arcsec, 4.717 * u.arcsec),
            (2.814 * u.arcsec, 4.841 * u.arcsec),
            (2.681 * u.arcsec, -4.792 * u.arcsec),
            (2.137 * u.arcsec, 7.681 * u.arcsec),
            (5.314 * u.arcsec, 1.819 * u.arcsec),
            (2.735 * u.arcsec, 6.677 * u.arcsec),
            (2.235 * u.arcsec, -5.757 * u.arcsec),
            (0.04 * u.arcsec, 5.575 * u.arcsec),
            (3.872 * u.arcsec, -4.792 * u.arcsec),
            (5.434 * u.arcsec, 0.049 * u.arcsec),
            (4.938 * u.arcsec, -2.876 * u.arcsec),
            (1.084 * u.arcsec, -5.765 * u.arcsec),
            (4.903 * u.arcsec, -0.929 * u.arcsec),
            (5.42 * u.arcsec, -1.867 * u.arcsec),
            (6.588 * u.arcsec, -0.013 * u.arcsec),
            (6.009 * u.arcsec, -0.876 * u.arcsec),
            (4.912 * u.arcsec, 0.942 * u.arcsec),
            (93.235 * u.arcsec, -1.159 * u.arcsec),
            (3.301 * u.arcsec, -5.761 * u.arcsec),
            (3.323 * u.arcsec, 5.584 * u.arcsec),
            (5.978 * u.arcsec, 0.934 * u.arcsec),
            (1.668 * u.arcsec, -6.646 * u.arcsec),
            (1.659 * u.arcsec, -8.469 * u.arcsec),
            (1.097 * u.arcsec, -7.571 * u.arcsec),
            (2.743 * u.arcsec, -6.642 * u.arcsec),
            (1.062 * u.arcsec, 3.779 * u.arcsec),
            (3.274 * u.arcsec, 3.752 * u.arcsec),
            (2.19 * u.arcsec, 0.044 * u.arcsec),
            (2.226 * u.arcsec, -1.951 * u.arcsec),
            (0.615 * u.arcsec, -2.929 * u.arcsec),
            (2.743 * u.arcsec, -0.925 * u.arcsec),
            (2.217 * u.arcsec, 3.774 * u.arcsec),
            (1.646 * u.arcsec, 2.823 * u.arcsec),
            (96.956 * u.arcsec, -0.252 * u.arcsec),
            (-0.009 * u.arcsec, -3.85 * u.arcsec),
            (0.606 * u.arcsec, 4.624 * u.arcsec),
            (2.699 * u.arcsec, 0.934 * u.arcsec),
            (0.566 * u.arcsec, -0.969 * u.arcsec),
            (1.659 * u.arcsec, -2.934 * u.arcsec),
            (3.801 * u.arcsec, 2.819 * u.arcsec),
            (3.27 * u.arcsec, 0.049 * u.arcsec),
            (1.615 * u.arcsec, -4.792 * u.arcsec),
            (98.049 * u.arcsec, -0.204 * u.arcsec),
            (3.832 * u.arcsec, -2.925 * u.arcsec),
            (2.199 * u.arcsec, -3.81 * u.arcsec),
            (1.677 * u.arcsec, -0.92 * u.arcsec),
            (3.845 * u.arcsec, -0.929 * u.arcsec),
            (4.327 * u.arcsec, 0.044 * u.arcsec),
            (1.066 * u.arcsec, -3.863 * u.arcsec),
            (1.102 * u.arcsec, 1.823 * u.arcsec),
            (0.571 * u.arcsec, 0.889 * u.arcsec),
            (95.938 * u.arcsec, -0.252 * u.arcsec),
            (0.018 * u.arcsec, 3.73 * u.arcsec),
            (2.704 * u.arcsec, 2.823 * u.arcsec),
            (0.562 * u.arcsec, -4.788 * u.arcsec),
            (2.08 * u.arcsec, 1.77 * u.arcsec),
            (3.261 * u.arcsec, -3.858 * u.arcsec),
            (4.323 * u.arcsec, -1.907 * u.arcsec),
            (2.726 * u.arcsec, -2.929 * u.arcsec),
            (1.712 * u.arcsec, 4.832 * u.arcsec),
            (1.093 * u.arcsec, -1.907 * u.arcsec),
            (4.208 * u.arcsec, 1.774 * u.arcsec),
            (1.137 * u.arcsec, 0.049 * u.arcsec),
            (3.265 * u.arcsec, -1.907 * u.arcsec),
            (0.504 * u.arcsec, 2.708 * u.arcsec),
            (0.0 * u.arcsec, 1.832 * u.arcsec),
            (1.633 * u.arcsec, 0.889 * u.arcsec),
            (3.85 * u.arcsec, 0.938 * u.arcsec),
            (-0.013 * u.arcsec, -1.903 * u.arcsec),
            (3.15 * u.arcsec, 1.77 * u.arcsec),
        ]


class PF0300N001(SMI):
    def fiber_radius(self) -> Angle:
        return 0.665 * u.arcsec

    def target_bundle_center_row_extension(self) -> Angle:
        return 3.8 * self.fiber_radius()

    def fibers(self) -> List[Tuple[Angle, Angle]]:
        return [
            (-8.98 * u.arcsec, -3.23 * u.arcsec),
            (81.18 * u.arcsec, -2.3 * u.arcsec),
            (79.58 * u.arcsec, -2.3 * u.arcsec),
            (-8.98 * u.arcsec, 2.48 * u.arcsec),
            (-6.46 * u.arcsec, 6.76 * u.arcsec),
            (-8.09 * u.arcsec, -7.38 * u.arcsec),
            (80.42 * u.arcsec, -0.97 * u.arcsec),
            (-8.09 * u.arcsec, 3.98 * u.arcsec),
            (-8.93 * u.arcsec, 5.35 * u.arcsec),
            (-7.25 * u.arcsec, -8.93 * u.arcsec),
            (-9.02 * u.arcsec, -0.35 * u.arcsec),
            (-7.3 * u.arcsec, 5.35 * u.arcsec),
            (-9.82 * u.arcsec, -1.81 * u.arcsec),
            (-4.82 * u.arcsec, 6.68 * u.arcsec),
            (-9.82 * u.arcsec, 3.98 * u.arcsec),
            (-8.14 * u.arcsec, -4.64 * u.arcsec),
            (-6.5 * u.arcsec, -7.47 * u.arcsec),
            (-10.66 * u.arcsec, -0.35 * u.arcsec),
            (-9.77 * u.arcsec, -4.64 * u.arcsec),
            (-8.89 * u.arcsec, -5.97 * u.arcsec),
            (-7.25 * u.arcsec, -6.06 * u.arcsec),
            (-4.02 * u.arcsec, 8.14 * u.arcsec),
            (-7.3 * u.arcsec, 8.14 * u.arcsec),
            (-8.05 * u.arcsec, 6.76 * u.arcsec),
            (-10.57 * u.arcsec, -3.23 * u.arcsec),
            (-9.86 * u.arcsec, 1.06 * u.arcsec),
            (-5.66 * u.arcsec, 8.14 * u.arcsec),
            (-10.66 * u.arcsec, 2.48 * u.arcsec),
            (-5.66 * u.arcsec, -8.93 * u.arcsec),
            (82.81 * u.arcsec, 0.49 * u.arcsec),
            (82.02 * u.arcsec, 1.95 * u.arcsec),
            (-13.09 * u.arcsec, -1.72 * u.arcsec),
            (82.86 * u.arcsec, -2.3 * u.arcsec),
            (-13.09 * u.arcsec, 1.06 * u.arcsec),
            (-12.25 * u.arcsec, -0.35 * u.arcsec),
            (83.61 * u.arcsec, -0.97 * u.arcsec),
            (-10.52 * u.arcsec, 8.14 * u.arcsec),
            (-12.25 * u.arcsec, 5.44 * u.arcsec),
            (-13.88 * u.arcsec, 2.56 * u.arcsec),
            (-8.89 * u.arcsec, 8.14 * u.arcsec),
            (-9.73 * u.arcsec, 6.76 * u.arcsec),
            (-11.41 * u.arcsec, 3.98 * u.arcsec),
            (-13.0 * u.arcsec, 3.98 * u.arcsec),
            (-11.5 * u.arcsec, 1.06 * u.arcsec),
            (-11.36 * u.arcsec, 6.9 * u.arcsec),
            (-11.45 * u.arcsec, -1.81 * u.arcsec),
            (-10.52 * u.arcsec, -8.75 * u.arcsec),
            (-8.98 * u.arcsec, -8.75 * u.arcsec),
            (-10.52 * u.arcsec, -5.97 * u.arcsec),
            (-12.2 * u.arcsec, -3.23 * u.arcsec),
            (-12.16 * u.arcsec, -5.97 * u.arcsec),
            (-13.93 * u.arcsec, -0.35 * u.arcsec),
            (-13.88 * u.arcsec, -3.14 * u.arcsec),
            (-10.61 * u.arcsec, 5.44 * u.arcsec),
            (-11.32 * u.arcsec, -7.38 * u.arcsec),
            (-12.25 * u.arcsec, 2.48 * u.arcsec),
            (-11.41 * u.arcsec, -4.64 * u.arcsec),
            (-13.09 * u.arcsec, -4.64 * u.arcsec),
            (-9.73 * u.arcsec, -7.38 * u.arcsec),
            (-2.43 * u.arcsec, -0.53 * u.arcsec),
            (-0.84 * u.arcsec, -3.4 * u.arcsec),
            (-0.71 * u.arcsec, -6.15 * u.arcsec),
            (-0.84 * u.arcsec, 2.34 * u.arcsec),
            (-1.59 * u.arcsec, 0.97 * u.arcsec),
            (-0.84 * u.arcsec, -0.53 * u.arcsec),
            (73.08 * u.arcsec, -2.3 * u.arcsec),
            (-1.72 * u.arcsec, -4.82 * u.arcsec),
            (-3.32 * u.arcsec, 0.97 * u.arcsec),
            (-2.43 * u.arcsec, -3.4 * u.arcsec),
            (-1.77 * u.arcsec, -1.9 * u.arcsec),
            (76.36 * u.arcsec, -2.3 * u.arcsec),
            (75.56 * u.arcsec, -3.67 * u.arcsec),
            (-3.32 * u.arcsec, -1.9 * u.arcsec),
            (-4.91 * u.arcsec, -1.9 * u.arcsec),
            (-0.75 * u.arcsec, 5.22 * u.arcsec),
            (-4.07 * u.arcsec, -3.32 * u.arcsec),
            (-4.86 * u.arcsec, 1.06 * u.arcsec),
            (73.88 * u.arcsec, -0.97 * u.arcsec),
            (71.49 * u.arcsec, -2.3 * u.arcsec),
            (0.0 * u.arcsec, 3.89 * u.arcsec),
            (-2.39 * u.arcsec, -6.15 * u.arcsec),
            (72.33 * u.arcsec, -3.67 * u.arcsec),
            (-4.07 * u.arcsec, -0.53 * u.arcsec),
            (-3.23 * u.arcsec, 3.89 * u.arcsec),
            (-1.59 * u.arcsec, 3.89 * u.arcsec),
            (-2.39 * u.arcsec, 5.22 * u.arcsec),
            (-5.75 * u.arcsec, -0.44 * u.arcsec),
            (74.63 * u.arcsec, 0.53 * u.arcsec),
            (0.0 * u.arcsec, 0.97 * u.arcsec),
            (-4.11 * u.arcsec, 2.39 * u.arcsec),
            (-2.43 * u.arcsec, 2.39 * u.arcsec),
            (-0.71 * u.arcsec, -8.98 * u.arcsec),
            (-4.91 * u.arcsec, -4.73 * u.arcsec),
            (-7.34 * u.arcsec, -0.44 * u.arcsec),
            (-4.02 * u.arcsec, -6.15 * u.arcsec),
            (-5.62 * u.arcsec, -6.06 * u.arcsec),
            (-6.59 * u.arcsec, -1.9 * u.arcsec),
            (-5.7 * u.arcsec, 2.39 * u.arcsec),
            (-3.23 * u.arcsec, -7.47 * u.arcsec),
            (-1.5 * u.arcsec, -7.47 * u.arcsec),
            (-5.7 * u.arcsec, -3.32 * u.arcsec),
            (-6.59 * u.arcsec, 1.06 * u.arcsec),
            (-7.38 * u.arcsec, 2.39 * u.arcsec),
            (-7.3 * u.arcsec, -3.18 * u.arcsec),
            (-8.14 * u.arcsec, -1.9 * u.arcsec),
            (-4.07 * u.arcsec, 5.22 * u.arcsec),
            (78.83 * u.arcsec, -3.67 * u.arcsec),
            (-5.62 * u.arcsec, 5.35 * u.arcsec),
            (0.04 * u.arcsec, 6.68 * u.arcsec),
            (-6.54 * u.arcsec, -4.73 * u.arcsec),
            (-6.54 * u.arcsec, 3.89 * u.arcsec),
            (-2.34 * u.arcsec, -8.98 * u.arcsec),
            (-2.39 * u.arcsec, 8.14 * u.arcsec),
            (-4.91 * u.arcsec, 3.89 * u.arcsec),
            (-4.82 * u.arcsec, -7.47 * u.arcsec),
            (78.74 * u.arcsec, -0.97 * u.arcsec),
            (77.15 * u.arcsec, 1.95 * u.arcsec),
            (-1.64 * u.arcsec, 6.68 * u.arcsec),
            (77.81 * u.arcsec, 0.49 * u.arcsec),
            (-8.18 * u.arcsec, 1.06 * u.arcsec),
            (-3.18 * u.arcsec, 6.68 * u.arcsec),
            (-0.71 * u.arcsec, 8.14 * u.arcsec),
            (12.34 * u.arcsec, -6.06 * u.arcsec),
            (12.29 * u.arcsec, -3.27 * u.arcsec),
            (13.09 * u.arcsec, 0.97 * u.arcsec),
            (13.13 * u.arcsec, -4.69 * u.arcsec),
            (13.97 * u.arcsec, -3.27 * u.arcsec),
            (-82.24 * u.arcsec, -1.59 * u.arcsec),
            (10.74 * u.arcsec, -5.97 * u.arcsec),
            (10.79 * u.arcsec, -8.98 * u.arcsec),
            (-81.35 * u.arcsec, 2.7 * u.arcsec),
            (11.41 * u.arcsec, 3.8 * u.arcsec),
            (11.45 * u.arcsec, 6.68 * u.arcsec),
            (13.04 * u.arcsec, 3.8 * u.arcsec),
            (13.09 * u.arcsec, -1.9 * u.arcsec),
            (13.93 * u.arcsec, -0.49 * u.arcsec),
            (12.29 * u.arcsec, 2.34 * u.arcsec),
            (9.11 * u.arcsec, -8.98 * u.arcsec),
            (12.29 * u.arcsec, -0.49 * u.arcsec),
            (10.66 * u.arcsec, 8.14 * u.arcsec),
            (11.45 * u.arcsec, 0.97 * u.arcsec),
            (11.5 * u.arcsec, -1.9 * u.arcsec),
            (11.45 * u.arcsec, -4.69 * u.arcsec),
            (11.58 * u.arcsec, -7.43 * u.arcsec),
            (9.82 * u.arcsec, 6.68 * u.arcsec),
            (9.9 * u.arcsec, -7.43 * u.arcsec),
            (12.34 * u.arcsec, 5.22 * u.arcsec),
            (10.66 * u.arcsec, 5.22 * u.arcsec),
            (9.02 * u.arcsec, 8.14 * u.arcsec),
            (13.93 * u.arcsec, 2.34 * u.arcsec),
            (8.22 * u.arcsec, 3.8 * u.arcsec),
            (-79.67 * u.arcsec, -0.22 * u.arcsec),
            (-80.64 * u.arcsec, -1.59 * u.arcsec),
            (-78.92 * u.arcsec, 1.24 * u.arcsec),
            (-78.88 * u.arcsec, -1.59 * u.arcsec),
            (8.27 * u.arcsec, -7.43 * u.arcsec),
            (8.14 * u.arcsec, 6.68 * u.arcsec),
            (10.7 * u.arcsec, 2.34 * u.arcsec),
            (10.7 * u.arcsec, -3.27 * u.arcsec),
            (9.02 * u.arcsec, -3.27 * u.arcsec),
            (9.82 * u.arcsec, -4.69 * u.arcsec),
            (9.82 * u.arcsec, 3.8 * u.arcsec),
            (5.0 * u.arcsec, -7.56 * u.arcsec),
            (10.61 * u.arcsec, -0.49 * u.arcsec),
            (9.06 * u.arcsec, 5.22 * u.arcsec),
            (4.2 * u.arcsec, 8.14 * u.arcsec),
            (7.52 * u.arcsec, -8.98 * u.arcsec),
            (4.2 * u.arcsec, -8.98 * u.arcsec),
            (5.84 * u.arcsec, -8.98 * u.arcsec),
            (6.54 * u.arcsec, 6.68 * u.arcsec),
            (6.68 * u.arcsec, -7.52 * u.arcsec),
            (8.14 * u.arcsec, -4.69 * u.arcsec),
            (9.15 * u.arcsec, -6.06 * u.arcsec),
            (9.77 * u.arcsec, -1.9 * u.arcsec),
            (7.47 * u.arcsec, -6.06 * u.arcsec),
            (7.38 * u.arcsec, 5.22 * u.arcsec),
            (9.82 * u.arcsec, 0.97 * u.arcsec),
            (9.02 * u.arcsec, 2.34 * u.arcsec),
            (9.02 * u.arcsec, -0.49 * u.arcsec),
            (7.38 * u.arcsec, 8.14 * u.arcsec),
            (5.79 * u.arcsec, 8.14 * u.arcsec),
            (8.14 * u.arcsec, -1.9 * u.arcsec),
            (7.34 * u.arcsec, -3.27 * u.arcsec),
            (-77.37 * u.arcsec, -1.37 * u.arcsec),
            (2.52 * u.arcsec, -6.15 * u.arcsec),
            (2.52 * u.arcsec, 2.34 * u.arcsec),
            (6.54 * u.arcsec, -4.69 * u.arcsec),
            (8.18 * u.arcsec, 0.97 * u.arcsec),
            (4.91 * u.arcsec, -4.69 * u.arcsec),
            (0.93 * u.arcsec, -8.98 * u.arcsec),
            (3.32 * u.arcsec, -7.43 * u.arcsec),
            (0.04 * u.arcsec, -7.43 * u.arcsec),
            (2.56 * u.arcsec, -8.98 * u.arcsec),
            (1.68 * u.arcsec, -7.43 * u.arcsec),
            (5.79 * u.arcsec, 5.22 * u.arcsec),
            (-76.13 * u.arcsec, 1.02 * u.arcsec),
            (1.68 * u.arcsec, 6.68 * u.arcsec),
            (4.11 * u.arcsec, -6.15 * u.arcsec),
            (5.75 * u.arcsec, -6.06 * u.arcsec),
            (5.66 * u.arcsec, 2.34 * u.arcsec),
            (7.38 * u.arcsec, 2.34 * u.arcsec),
            (6.5 * u.arcsec, 0.97 * u.arcsec),
            (6.5 * u.arcsec, -1.9 * u.arcsec),
            (4.95 * u.arcsec, 3.8 * u.arcsec),
            (5.75 * u.arcsec, -3.27 * u.arcsec),
            (7.34 * u.arcsec, -0.49 * u.arcsec),
            (-76.53 * u.arcsec, 2.48 * u.arcsec),
            (2.56 * u.arcsec, 8.14 * u.arcsec),
            (3.32 * u.arcsec, 6.68 * u.arcsec),
            (6.59 * u.arcsec, 3.8 * u.arcsec),
            (4.91 * u.arcsec, 6.68 * u.arcsec),
            (0.88 * u.arcsec, 8.14 * u.arcsec),
            (-77.99 * u.arcsec, -2.92 * u.arcsec),
            (4.16 * u.arcsec, 5.22 * u.arcsec),
            (1.68 * u.arcsec, 3.89 * u.arcsec),
            (-74.1 * u.arcsec, 1.24 * u.arcsec),
            (-72.51 * u.arcsec, -1.59 * u.arcsec),
            (-71.58 * u.arcsec, -3.09 * u.arcsec),
            (-70.83 * u.arcsec, -1.59 * u.arcsec),
            (-74.85 * u.arcsec, -3.09 * u.arcsec),
            (2.48 * u.arcsec, 5.22 * u.arcsec),
            (0.88 * u.arcsec, 2.34 * u.arcsec),
            (0.93 * u.arcsec, 5.22 * u.arcsec),
            (1.59 * u.arcsec, 0.97 * u.arcsec),
            (-73.26 * u.arcsec, -0.22 * u.arcsec),
            (4.11 * u.arcsec, -0.49 * u.arcsec),
            (1.64 * u.arcsec, -4.82 * u.arcsec),
            (4.16 * u.arcsec, 2.34 * u.arcsec),
            (5.75 * u.arcsec, -0.49 * u.arcsec),
            (2.48 * u.arcsec, -3.4 * u.arcsec),
            (3.23 * u.arcsec, 0.97 * u.arcsec),
            (1.59 * u.arcsec, -1.9 * u.arcsec),
            (4.95 * u.arcsec, 0.97 * u.arcsec),
            (0.84 * u.arcsec, -0.57 * u.arcsec),
            (4.95 * u.arcsec, -1.9 * u.arcsec),
            (4.07 * u.arcsec, -3.27 * u.arcsec),
            (-0.09 * u.arcsec, -1.9 * u.arcsec),
            (3.27 * u.arcsec, -1.9 * u.arcsec),
            (3.27 * u.arcsec, -4.82 * u.arcsec),
            (0.93 * u.arcsec, -6.15 * u.arcsec),
            (2.48 * u.arcsec, -0.49 * u.arcsec),
            (0.84 * u.arcsec, -3.4 * u.arcsec),
            (3.36 * u.arcsec, 3.89 * u.arcsec),
            (-0.04 * u.arcsec, -4.82 * u.arcsec),
        ]
