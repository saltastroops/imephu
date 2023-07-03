from astropy import units as u
from imephu.annotation.general import CrosshairsAnnotation
from imephu.finder_chart import FinderChart


def test_crosshairs_annotation(fits_file, fits_center, check_finder):
    """Test the crosshairs annotation."""
    finder_chart = FinderChart(fits_file)
    crosshairs_annotation = CrosshairsAnnotation(
        center=fits_center,
        size=4 * u.arcmin,
        wcs=finder_chart.wcs,
        color="blue",
        alpha=0.2,
        linewidth=2,
    )
    finder_chart.add_annotation(crosshairs_annotation)
    check_finder(finder_chart)
