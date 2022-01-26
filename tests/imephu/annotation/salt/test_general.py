import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from imephu.annotation.general import TextAnnotation
from imephu.annotation.salt import general
from imephu.finder_chart import FinderChart


def test_title_annotation(fits_file, check_finder):
    """Test the title annotation."""
    finder_chart = FinderChart(fits_file)
    title_annotation = general.title(
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
    directions_annotation = general.directions(fits_center, finder_chart.wcs)
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
    position_angle_annotation = general.position_angle(
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
    survey_annotation = general.survey("POSS2/UKSTU Red", finder_chart.wcs)
    finder_chart.add_annotation(survey_annotation)
    check_finder(finder_chart)
