import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.general import TextAnnotation
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


@pytest.mark.parametrize("position_angle", [0 * u.deg, 60 * u.deg])
def test_salticam_slot_annotation(position_angle, fits_file, fits_center, check_finder):
    """Test the Salticam slot annotation."""
    finder_chart = FinderChart(fits_file)
    slot_annotation = salticam.slot_annotation(
        fits_center, position_angle, finder_chart.wcs
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"position angle: {position_angle}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(slot_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)
