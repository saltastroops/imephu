from astropy import units as u
from astropy.coordinates import SkyCoord
from imephu.annotation.general import EmptyAnnotation
from imephu.finder_chart import FinderChart


def test_circle_annotation(fits_file, fits_center, check_finder):
    """Test circle annotations."""
    finder_chart = FinderChart(fits_file)
    empty_annotation = EmptyAnnotation()
    finder_chart.add_annotation(empty_annotation)
    check_finder(finder_chart)


def test_circle_annotation_rotated(fits_file, check_finder, legend):
    """Test rotated circle annotations."""
    finder_chart = FinderChart(fits_file)
    empty_annotation = EmptyAnnotation()
    pivot = SkyCoord(ra="00h39m40s", dec=-60 * u.deg)
    angle = 90 * u.deg
    rotated_empty_annotation = empty_annotation.rotate(pivot, angle)
    finder_chart.add_annotation(rotated_empty_annotation)
    check_finder(finder_chart)


def test_empty_annotation_translated(fits_file, fits_center, check_finder, legend):
    """Test translated circle annotations."""
    finder_chart = FinderChart(fits_file)
    empty_annotation = EmptyAnnotation()
    displacement = (2.5, -4) * u.arcmin
    translated_empty_annotation = empty_annotation.translate(displacement)
    finder_chart.add_annotation(translated_empty_annotation)
    check_finder(finder_chart)
