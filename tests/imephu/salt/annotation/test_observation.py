import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from imephu.finder_chart import FinderChart
from imephu.salt.annotation import (
    GeneralProperties,
    MagnitudeRange,
    Target,
    hrs_observation_annotation,
    rss_fabry_perot_observation_annotation,
    rss_imaging_observation_annotation,
    rss_longslit_observation_annotation,
    rss_mos_observation_annotation,
    salticam_observation_annotation,
)
from imephu.salt.utils import MosMaskSlit, MosMask

POSITION_ANGLE = Angle(20 * u.deg)


def _general_properties(target_position, wcs):
    magnitude_range = MagnitudeRange("V", 16.4, 17.1)
    return GeneralProperties(
        target=Target(
            name="Magrathea", position=target_position, magnitude_range=magnitude_range
        ),
        position_angle=POSITION_ANGLE,
        automated_position_angle=False,
        proposal_code="2022-1-SCI-042",
        pi_family_name="Adams",
        survey="POSS2/UKSTU Red",
        wcs=wcs,
    )


@pytest.mark.parametrize("is_slot_mode", [False, True])
def test_salticam_observation_annotation(
    is_slot_mode, fits_file, fits_center, check_finder
):
    """Test the annotation for Salticam observations."""
    finder_chart = FinderChart(fits_file)
    salticam = salticam_observation_annotation(
        general=_general_properties(fits_center, finder_chart.wcs),
        is_slot_mode=is_slot_mode,
    )
    finder_chart.add_annotation(salticam)
    check_finder(finder_chart)


@pytest.mark.parametrize("is_slot_mode", [False, True])
def test_rss_imaging_observation_annotation(
    is_slot_mode, fits_file, fits_center, check_finder
):
    """Test the annotation for RSS imaging observations."""
    finder_chart = FinderChart(fits_file)
    rss = rss_imaging_observation_annotation(
        general=_general_properties(fits_center, finder_chart.wcs),
        is_slot_mode=is_slot_mode,
    )
    finder_chart.add_annotation(rss)
    check_finder(finder_chart)


def test_rss_longslit_observation_annotation(fits_file, fits_center, check_finder):
    """Test the annotation for RSS longslit observations."""
    finder_chart = FinderChart(fits_file)
    rss = rss_longslit_observation_annotation(
        general=_general_properties(fits_center, finder_chart.wcs),
        slit_width=4 * u.arcsec,
        slit_height=8 * u.arcmin,
    )
    finder_chart.add_annotation(rss)
    check_finder(finder_chart)


def test_rss_mos_observation_annotation(fits_file, fits_center, mos_mask_xml, check_finder):
    """Test the annotation for RSS MOS observations."""
    finder_chart = FinderChart(fits_file)
    reference_stars = [
        SkyCoord(ra="00h39m30s", dec="-60d02m00s"),
        SkyCoord(ra="00h40m30s", dec="-60d02m00s"),
    ]
    slits = [
        MosMaskSlit(
            center=SkyCoord(ra="00h39m50s", dec="-60d03m00s"),
            width=Angle(2 * u.arcsec),
            height=Angle(10 * u.arcsec),
            tilt=0 * u.deg,
        ),
        MosMaskSlit(
            center=SkyCoord(ra="00h40m10s", dec="-60d03m00s"),
            width=Angle(4 * u.arcsec),
            height=Angle(20 * u.arcsec),
            tilt=-45 * u.deg,
        ),
    ]
    xml = mos_mask_xml(
        center=fits_center,
        position_angle=POSITION_ANGLE,
        reference_stars=reference_stars,
        slits=slits,
    )
    mos_mask = MosMask(xml)
    rss = rss_mos_observation_annotation(general=_general_properties(fits_center, finder_chart.wcs),
                                         mos_mask=mos_mask)
    finder_chart.add_annotation(rss)
    check_finder(finder_chart)


def test_rss_fabry_perot_observation_annotation(fits_file, fits_center, check_finder):
    """Test the annotation for RSS Fabry-Pérot observations."""
    finder_chart = FinderChart(fits_file)
    rss = rss_fabry_perot_observation_annotation(
        general=_general_properties(fits_center, finder_chart.wcs),
    )
    finder_chart.add_annotation(rss)
    check_finder(finder_chart)


def test_hrs_observation_annotation(fits_file, fits_center, check_finder):
    """Test the annotation for an HRS observation."""
    finder_chart = FinderChart(fits_file)
    hrs = hrs_observation_annotation(_general_properties(fits_center, finder_chart.wcs))
    finder_chart.add_annotation(hrs)
    check_finder(finder_chart)
