import pytest
from astropy.coordinates import Angle
import astropy.units as u

from imephu.annotation.general.scale_bar import (
    _ScaleBarParameters,
    ScaleBarLineAnnotation,
)
from imephu.finder_chart import FinderChart


@pytest.mark.parametrize(
    "minimum_length, pixel_scale, expected_angle, expected_units",
    [
        (1, 0.1 * u.arcsec, 1, "arcsec"),
        (1, 1 * u.arcsec, 1, "arcsec"),
        (1, 3.7 * u.arcsec, 4, "arcsec"),
        (1, 9.1 * u.arcsec, 10, "arcsec"),
        (1, 10.1 * u.arcsec, 20, "arcsec"),
        (1, 17 * u.arcsec, 20, "arcsec"),
        (1, 22 * u.arcsec, 30, "arcsec"),
        (1, 44 * u.arcsec, 50, "arcsec"),
        (1, 51 * u.arcsec, 1, "arcmin"),
        (1, 1 * u.arcmin, 1, "arcmin"),
        (1, 3.7 * u.arcmin, 4, "arcmin"),
        (1, 9.1 * u.arcmin, 10, "arcmin"),
        (1, 10.1 * u.arcmin, 20, "arcmin"),
        (1, 17 * u.arcmin, 20, "arcmin"),
        (1, 22 * u.arcmin, 30, "arcmin"),
        (1, 101 * u.arcmin, 200, "arcmin"),
    ],
)
def test_angle_and_unit_parameters(
    minimum_length: float,
    pixel_scale: Angle,
    expected_angle: float,
    expected_units: str,
) -> None:
    parameters = ScaleBarLineAnnotation._parameters(minimum_length, pixel_scale)
    assert parameters.angle == pytest.approx(expected_angle)
    assert parameters.units == expected_units


def test_pixels_parameter():
    parameters = ScaleBarLineAnnotation._parameters(100, 1.5 * u.arcmin)
    assert parameters.pixels == pytest.approx(133.33333)


def test_scale_bar_annotation(fits_file, check_finder):
    """Test rectangle annotations."""
    finder_chart = FinderChart(fits_file)
    rectangle_annotation = ScaleBarLineAnnotation(
        left_edge=(50, 50),
        minimum_length=100,
        wcs=finder_chart.wcs,
        color="black",
    )
    finder_chart.add_annotation(rectangle_annotation)
    check_finder(finder_chart)
