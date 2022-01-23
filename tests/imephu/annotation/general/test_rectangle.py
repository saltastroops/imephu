import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.general import CircleAnnotation, TextAnnotation
from imephu.annotation.general.rectangle import RectangleAnnotation
from imephu.finder_chart import FinderChart


def test_rectangle_annotation(fits_file, check_finder):
    """Test rectangle annotations."""
    finder_chart = FinderChart(fits_file)
    rectangle_annotation = RectangleAnnotation(
        SkyCoord(ra="00h39m30s", dec="-59d59m00s"),
        width=150 * u.arcsec,
        height=2 * u.arcmin,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="blue",
        alpha=0.2,
    )
    finder_chart.add_annotation(rectangle_annotation)
    check_finder(finder_chart)


@pytest.mark.parametrize("angle", [0 * u.deg, 45 * u.deg])
def test_rectangle_annotation_rotated(angle, fits_file, fits_center, check_finder):
    """Test rotated rectangle annotations."""
    finder_chart = FinderChart(fits_file)
    rectangle_annotation = RectangleAnnotation(
        fits_center,
        width=300 * u.arcsec,
        height=150 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="gray",
        alpha=0.2,
    )
    rotated_rectangle_annotation = rectangle_annotation.rotate(fits_center, angle)
    rotated_rectangle_annotation._kwargs["facecolor"] = "blue"
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"Rotated by {angle.to_value(u.deg)} deg",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    pivot_marker = CircleAnnotation(
        fits_center,
        12 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="orange",
        alpha=0.7,
    )
    finder_chart.add_annotation(pivot_marker)
    finder_chart.add_annotation(rectangle_annotation)
    finder_chart.add_annotation(rotated_rectangle_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize("displacement", [(0, 0) * u.arcmin, (2.5, -2) * u.arcmin])
def test_rectangle_annotation_translated(
    displacement, fits_file, fits_center, check_finder
):
    """Test translated rectangle annotations."""
    finder_chart = FinderChart(fits_file)
    rectangle_annotation = RectangleAnnotation(
        fits_center,
        width=100 * u.arcsec,
        height=200 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="gray",
        alpha=0.2,
    )
    translated_rectangle_annotation = rectangle_annotation.translate(displacement)
    translated_rectangle_annotation._kwargs["facecolor"] = "blue"
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"Translated by {displacement.to_value(u.arcmin)} arcmin",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(rectangle_annotation)
    finder_chart.add_annotation(translated_rectangle_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)
