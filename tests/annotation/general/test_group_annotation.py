import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from imephu.annotation.general import (
    CircleAnnotation,
    GroupAnnotation,
    RectangleAnnotation,
)
from imephu.finder_chart import FinderChart


def test_group_annotation(fits_file, check_finder):
    """Test group annotations."""
    finder_chart = FinderChart(fits_file)
    rectangle1_annotation = RectangleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-60d01m00s"),
        width=150 * u.arcsec,
        height=2 * u.arcmin,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="blue",
        alpha=0.2,
    )
    rectangle2_annotation = RectangleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-59d59m00s"),
        width=150 * u.arcsec,
        height=2 * u.arcmin,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="orange",
        alpha=0.2,
    )
    circle_annotation = CircleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-59d58m00s"),
        radius=75 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="green",
        facecolor="none",
        alpha=0.2,
        linewidth=4,
    )
    group_annotation = GroupAnnotation([rectangle1_annotation, rectangle2_annotation])
    group_annotation.add_item(circle_annotation)
    finder_chart.add_annotation(group_annotation)
    check_finder(finder_chart)


@pytest.mark.parametrize("angle", [0 * u.deg, 90 * u.deg])
def test_group_annotation_rotated(angle, fits_file, fits_center, check_finder, legend):
    """Test rotated group annotations."""
    finder_chart = FinderChart(fits_file)
    rectangle1_annotation = RectangleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-60d01m00s"),
        width=150 * u.arcsec,
        height=2 * u.arcmin,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="blue",
        alpha=0.2,
    )
    rectangle2_annotation = RectangleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-59d59m00s"),
        width=150 * u.arcsec,
        height=2 * u.arcmin,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="orange",
        alpha=0.2,
    )
    circle_annotation = CircleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-59d58m00s"),
        radius=75 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="green",
        facecolor="none",
        alpha=0.2,
        linewidth=4,
    )
    group_annotation = GroupAnnotation(
        [rectangle1_annotation, rectangle2_annotation, circle_annotation]
    )
    rotated_group_annotation = group_annotation.rotate(fits_center, angle)
    pivot_marker = CircleAnnotation(
        fits_center,
        12 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="orange",
        alpha=0.7,
    )
    finder_chart.add_annotation(pivot_marker)
    finder_chart.add_annotation(group_annotation)
    finder_chart.add_annotation(rotated_group_annotation)
    finder_chart.add_annotation(
        legend(
            f"Rotated by {angle.to_value(u.deg)} deg",
            wcs=finder_chart.wcs,
        )
    )
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "displacement", [Angle([0, 0] * u.arcmin), Angle([-5, -2] * u.arcmin)]
)
def test_group_annotation_translated(displacement, fits_file, check_finder, legend):
    """Test translated group annotations."""
    finder_chart = FinderChart(fits_file)
    rectangle1_annotation = RectangleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-60d01m00s"),
        width=150 * u.arcsec,
        height=2 * u.arcmin,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="blue",
        alpha=0.2,
    )
    rectangle2_annotation = RectangleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-59d59m00s"),
        width=150 * u.arcsec,
        height=2 * u.arcmin,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="orange",
        alpha=0.2,
    )
    circle_annotation = CircleAnnotation(
        SkyCoord(ra="00h40m30s", dec="-59d58m00s"),
        radius=75 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="green",
        facecolor="none",
        alpha=0.2,
        linewidth=4,
    )
    group_annotation = GroupAnnotation(
        [rectangle1_annotation, rectangle2_annotation, circle_annotation]
    )
    translated_group_annotation = group_annotation.translate(displacement)
    finder_chart.add_annotation(group_annotation)
    finder_chart.add_annotation(translated_group_annotation)
    finder_chart.add_annotation(
        legend(
            f"Translated by {displacement.to_value(u.arcmin)} arcmin",
            wcs=finder_chart.wcs,
        )
    )
    check_finder(finder_chart)
