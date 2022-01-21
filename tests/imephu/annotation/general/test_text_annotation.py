import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.general import TextAnnotation
from imephu.finder_chart import FinderChart

fits_center = SkyCoord(ra=10 * u.deg, dec=-60 * u.deg)


def test_text_annotation(fits_file, check_finder):
    """Test text annotations."""
    finder_chart = FinderChart(fits_file)
    text_annotation = TextAnnotation(
        fits_center,
        "Green, 18 pt, right-aligned",
        wcs=finder_chart.wcs,
        color="green",
        fontsize=18,
        horizontalalignment="right",
    )
    finder_chart.add_annotation(text_annotation)
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "center,angle",
    [
        (SkyCoord(ra="00h39m40s", dec=-60 * u.deg), 0 * u.deg),
        (SkyCoord(ra="00h39m40s", dec=-60 * u.deg), -90 * u.deg),
    ],
)
def test_text_annotation_rotated(center, angle, fits_file, check_finder):
    """Test rotated text annotations."""
    finder_chart = FinderChart(fits_file)
    non_rotated_text_annotation = TextAnnotation(
        fits_center, "Non-rotated text", wcs=finder_chart.wcs
    )
    s = f"Text center rotated by {angle}"
    text_annotation = TextAnnotation(fits_center, s, wcs=finder_chart.wcs)
    rotated_text_annotation = text_annotation.rotate(center, angle)
    finder_chart.add_annotation(non_rotated_text_annotation)
    finder_chart.add_annotation(rotated_text_annotation)
    check_finder(finder_chart)


@pytest.mark.parametrize("displacement", [(0, 0) * u.arcmin, (2.5, -4) * u.arcmin])
def test_text_annotation_translated(displacement, fits_file, check_finder):
    """Test translated texzt annotations."""
    s = f"Text start displaced by {displacement}"
    finder_chart = FinderChart(fits_file)
    text_annotation = TextAnnotation(
        fits_center, s, wcs=finder_chart.wcs, horizontalalignment="left"
    )
    translated_text_annotation = text_annotation.translate(displacement)
    finder_chart.add_annotation(translated_text_annotation)
    check_finder(finder_chart)
