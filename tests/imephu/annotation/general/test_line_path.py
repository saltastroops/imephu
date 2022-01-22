import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.general import TextAnnotation, CircleAnnotation
from imephu.annotation.general.line_path import LinePathAnnotation
from imephu.finder_chart import FinderChart


def test_line_path_annotation(fits_file, check_finder):
    finder_chart = FinderChart(fits_file)
    closed_path_vertices = [
        SkyCoord(ra="00h40m40s", dec="-59d56m00s"),
        SkyCoord(ra="00h40m20s", dec="-59d55m00s" ),
        SkyCoord(ra="00h40m00s", dec="-60d00m00s"),
    ]
    closed_path_annotation = LinePathAnnotation(closed_path_vertices, wcs=finder_chart.wcs, edgecolor="none", facecolor="green", alpha=0.2)
    open_path_vertices = [
        SkyCoord(ra="00h39m20s", dec="-60d04m00s"),
        SkyCoord(ra="00h39m40s", dec="-60d05m00s" ),
        SkyCoord(ra="00h40m00s", dec="-60d00m00s"),
    ]
    open_path_annotation = LinePathAnnotation(open_path_vertices, wcs=finder_chart.wcs, closed=False, edgecolor="orange")
    finder_chart.add_annotation(closed_path_annotation)
    finder_chart.add_annotation(open_path_annotation)
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "pivot,angle",
    [
        (SkyCoord(ra="00h39m40s", dec=-60 * u.deg), 0 * u.deg),
        (SkyCoord(ra="00h39m40s", dec=-60 * u.deg), -90 * u.deg),
    ],
)
def test_line_path_annotation_rotated(pivot, angle, fits_file, check_finder):
    """Test rotated circle annotations."""
    finder_chart = FinderChart(fits_file)
    line_path_annotation = LinePathAnnotation(
        [SkyCoord(ra="00h39m50s", dec="-60d00m00s"),
         SkyCoord(ra="00h39m50s", dec="-60d01m00s"),
         SkyCoord(ra="00h40m00s", dec="-60d00m00s")],
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="gray",
        alpha=0.2
        )
    rotated_circle_annotation = line_path_annotation.rotate(pivot, angle)
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
    finder_chart.add_annotation(line_path_annotation)
    finder_chart.add_annotation(rotated_circle_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize("displacement", [(0, 0) * u.arcmin, (2.5, -4) * u.arcmin])
def test_line_path_annotation_translated(
        displacement, fits_file, fits_center, check_finder
):
    """Test translated circle annotations."""
    finder_chart = FinderChart(fits_file)
    line_path_annotation = LinePathAnnotation(
        [
            SkyCoord(ra="00h39m40s", dec="-59d58m00s"),
            SkyCoord(ra="00h39m50s", dec="-59d58m00s"),
            SkyCoord(ra="00h39m40s", dec="-59d59m00s")
        ],
        wcs=finder_chart.wcs,
        edgecolor="none",
        facecolor="gray",
        )
    translated_line_path_annotation = line_path_annotation.translate(displacement)
    translated_line_path_annotation._kwargs["color"] = "blue"
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"Translated by {displacement.to_value(u.arcmin)} arcmin",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(line_path_annotation)
    finder_chart.add_annotation(translated_line_path_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)
