import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.general import CircleAnnotation, TextAnnotation
from imephu.finder_chart import FinderChart


def test_circle_annotation(fits_file, fits_center, check_finder):
    """Test circle annotations."""
    finder_chart = FinderChart(fits_file)
    circle_annotation = CircleAnnotation(
        fits_center,
        radius=1 * u.arcmin,
        wcs=finder_chart.wcs,
        edgecolor="green",
    )
    finder_chart.add_annotation(circle_annotation)
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "pivot,angle",
    [
        (SkyCoord(ra="00h39m40s", dec=-60 * u.deg), 0 * u.deg),
        (SkyCoord(ra="00h39m40s", dec=-60 * u.deg), -90 * u.deg),
    ],
)
def test_circle_annotation_rotated(pivot, angle, fits_file, check_finder):
    """Test rotated circle annotations."""
    finder_chart = FinderChart(fits_file)
    circle_annotation = CircleAnnotation(
        SkyCoord(ra="00h40m00s", dec="-60d00m00s"),
        15 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="gray",
    )
    rotated_circle_annotation = circle_annotation.rotate(pivot, angle)
    rotated_circle_annotation._kwargs["color"] = "blue"
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"Rotated by {angle.to_value(u.deg)} deg",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    pivot_marker = CircleAnnotation(
        pivot,
        12 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="orange",
        alpha=0.7,
    )
    finder_chart.add_annotation(pivot_marker)
    finder_chart.add_annotation(circle_annotation)
    finder_chart.add_annotation(rotated_circle_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize("displacement", [(0, 0) * u.arcmin, (2.5, -4) * u.arcmin])
def test_circle_annotation_translated(
    displacement, fits_file, fits_center, check_finder
):
    """Test translated circle annotations."""
    finder_chart = FinderChart(fits_file)
    circle_annotation = CircleAnnotation(
        fits_center,
        15 * u.arcsec,
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="gray",
    )
    translated_circle_annotation = circle_annotation.translate(displacement)
    translated_circle_annotation._kwargs["color"] = "blue"
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"Translated by {displacement.to_value(u.arcmin)} arcmin",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(circle_annotation)
    finder_chart.add_annotation(translated_circle_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)
