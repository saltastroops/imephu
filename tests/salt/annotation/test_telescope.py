import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from imephu.annotation.general import TextAnnotation
from imephu.finder_chart import FinderChart
from imephu.salt.annotation import telescope


def test_title_annotation(fits_file, check_finder):
    """Test the title annotation."""
    finder_chart = FinderChart(fits_file)
    title_annotation = telescope.title_annotation(
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
    directions_annotation = telescope.directions_annotation(
        fits_center, finder_chart.wcs
    )
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
    position_angle_annotation = telescope.position_angle_annotation(
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
    """Test survey annotations."""
    finder_chart = FinderChart(fits_file)
    survey_annotation = telescope.survey_annotation("POSS2/UKSTU Red", finder_chart.wcs)
    finder_chart.add_annotation(survey_annotation)
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "bandpass,min_magnitude,max_magnitude",
    [("V", 16.3, 17), ("R", 18.5, 18.5), ("I", 15.3, 15.4)],
)
def test_magnitude_range_annotation(
    bandpass, min_magnitude, max_magnitude, fits_file, fits_center, check_finder
):
    """Test magnitude range annotations."""
    finder_chart = FinderChart(fits_file)
    magnitude_annotation = telescope.magnitude_range_annotation(
        bandpass=bandpass,
        min_magnitude=min_magnitude,
        max_magnitude=max_magnitude,
        wcs=finder_chart.wcs,
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"bandpass: {bandpass}, min magnitude: {min_magnitude}, "
        f"max magnitude: {max_magnitude}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(magnitude_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize("position_angle", [0 * u.deg, 60 * u.deg])
def test_salticam_slot_annotation(position_angle, fits_file, fits_center, check_finder):
    """Test the Salticam slot annotation."""
    finder_chart = FinderChart(fits_file)
    slot_annotation = telescope.slot_annotation(
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


def test_salt_base_annotations(fits_file, fits_center, check_finder):
    """Test SALT base annotations."""
    finder_chart = FinderChart(fits_file)
    base_annotations = telescope.base_annotations(
        target="Magrathea",
        proposal_code="2022-1-SCI-042",
        pi_family_name="Adams",
        position_angle=Angle(60.75 * u.deg),
        automated_position_angle=False,
        survey="POSS2/UKSTU Blue",
        fits_center=fits_center,
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(base_annotations)
    check_finder(finder_chart)
