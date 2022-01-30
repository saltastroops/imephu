import pytest
from astropy import units as u
from astropy.coordinates import Angle

from imephu.finder_chart import FinderChart
from imephu.salt.annotation import (
    GeneralProperties,
    MagnitudeRange,
    Target,
    rss_fabry_perot_annotation,
    rss_imaging_annotation,
    rss_longslit_annotation,
    salticam_annotation,
)


def _general_properties(target_position, wcs):
    magnitude_range = MagnitudeRange("V", 16.4, 17.1)
    return GeneralProperties(
        target=Target(
            name="Magrathea", position=target_position, magnitude_range=magnitude_range
        ),
        position_angle=Angle(20 * u.deg),
        automated_position_angle=False,
        proposal_code="2022-1-SCI-042",
        pi_family_name="Adams",
        survey="POSS2/UKSTU Red",
        wcs=wcs,
    )


@pytest.mark.parametrize("is_slot_mode", [False, True])
def test_salticam_annotation(is_slot_mode, fits_file, fits_center, check_finder):
    """Test the annotation for Salticam observations."""
    finder_chart = FinderChart(fits_file)
    salticam = salticam_annotation(
        general=_general_properties(fits_center, finder_chart.wcs),
        is_slot_mode=is_slot_mode,
    )
    finder_chart.add_annotation(salticam)
    check_finder(finder_chart)


@pytest.mark.parametrize("is_slot_mode", [False, True])
def test_rss_imaging_annotation(is_slot_mode, fits_file, fits_center, check_finder):
    """Test the annotation for RSS imaging observations."""
    finder_chart = FinderChart(fits_file)
    rss = rss_imaging_annotation(
        general=_general_properties(fits_center, finder_chart.wcs),
        is_slot_mode=is_slot_mode,
    )
    finder_chart.add_annotation(rss)
    check_finder(finder_chart)


def test_rss_longslit_annotation(fits_file, fits_center, check_finder):
    """Test the annotation for RSS longslit observations."""
    finder_chart = FinderChart(fits_file)
    rss = rss_longslit_annotation(
        general=_general_properties(fits_center, finder_chart.wcs),
        slit_width=4 * u.arcsec,
        slit_height=8 * u.arcmin,
    )
    finder_chart.add_annotation(rss)
    check_finder(finder_chart)


def test_rss_fabry_perot_annotation(fits_file, fits_center, check_finder):
    """Test the annotation for RSS Fabry-Pérot observations."""
    finder_chart = FinderChart(fits_file)
    rss = rss_fabry_perot_annotation(
        general=_general_properties(fits_center, finder_chart.wcs),
    )
    finder_chart.add_annotation(rss)
    check_finder(finder_chart)
