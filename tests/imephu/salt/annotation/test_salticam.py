from imephu.finder_chart import FinderChart
from imephu.salt.annotation import salticam


def test_salticam_field_of_view_annotation(fits_file, fits_center, check_finder):
    """Test the Salticam field of view annotation."""
    finder_chart = FinderChart(fits_file)
    salticam_fov_annotation = salticam.field_of_view_annotation(
        fits_center, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(salticam_fov_annotation)
    check_finder(finder_chart)
