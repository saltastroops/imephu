import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from imephu.annotation.general import TextAnnotation
from imephu.annotation.salt import general
from imephu.finder_chart import FinderChart


def test_title_annotation(fits_file, check_finder):
    """Test the title annotation."""
    finder_chart = FinderChart(fits_file)
    title_annotation = general.title_annotation(
        target="Magrathea",
        proposal_code="2022-1-SCI-042",
        pi_family_name="Adams",
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(title_annotation)
    check_finder(finder_chart)


def test_directions_annotation(fits_file, fits_center, check_finder):
    """Test the directions annotation."""
    finder_chart = FinderChart(fits_file)
    directions_annotation = general.directions_annotation(fits_center, finder_chart.wcs)
    finder_chart.add_annotation(directions_annotation)
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "angle,automated",
    [
        (Angle(1200 * u.arcmin), False),
        (Angle(-67.88 * u.deg), True),
        (Angle(127.14 * u.deg), False),
        (Angle(-77.3 * u.deg), True),
    ],
)
def test_position_angle_annotation(angle, automated, fits_file, check_finder):
    """Test the position angle annotation."""
    finder_chart = FinderChart(fits_file)
    position_angle_annotation = general.position_angle_annotation(
        angle, automated, finder_chart.wcs
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"angle: {angle}, automated: {automated}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(position_angle_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


def test_survey_annotation(fits_file, check_finder):
    finder_chart = FinderChart(fits_file)
    survey_annotation = general.survey_annotation("POSS2/UKSTU Red", finder_chart.wcs)
    finder_chart.add_annotation(survey_annotation)
    check_finder(finder_chart)


def test_salticam_field_of_view_annotation(fits_file, fits_center, check_finder):
    finder_chart = FinderChart(fits_file)
    salticam_fov_annotation = general.salticam_field_of_view_annotation(
        fits_center, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(salticam_fov_annotation)
    check_finder(finder_chart)


def test_rss_field_of_view_annotation(fits_file, fits_center, check_finder):
    finder_chart = FinderChart(fits_file)
    rss_fov_annotation = general.rss_field_of_view_annotation(
        fits_center, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(rss_fov_annotation)
    check_finder(finder_chart)


def test_salt_base_annotations(fits_file, fits_center, check_finder):
    finder_chart = FinderChart(fits_file)
    base_annotations = general.base_annotations(
        target="Magrathea",
        proposal_code="2022-1-SCI-042",
        pi_family_name="Adams",
        position_angle=Angle(60.75 * u.deg),
        automated_position_angle=False,
        survey_name="POSS2/UKSTU Blue",
        fits_center=fits_center,
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(base_annotations)
    check_finder(finder_chart)
