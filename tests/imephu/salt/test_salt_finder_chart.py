import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from imephu.salt.finder_chart import (
    GeneralProperties,
    MagnitudeRange,
    Target,
    hrs_finder_chart,
    rss_fabry_perot_finder_chart,
    rss_imaging_finder_chart,
    rss_longslit_finder_chart,
    rss_mos_finder_chart,
    salticam_finder_chart,
)
from imephu.salt.utils import MosMask, MosMaskSlit

POSITION_ANGLE = Angle(20 * u.deg)


def _general_properties(target_position):
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
    )


@pytest.mark.parametrize("is_slot_mode", [False, True])
def test_salticam_finder_chart(
    is_slot_mode, fits_file, fits_center, check_finder, mock_salt_load_fits
):
    """Test the finder chart for Salticam observations."""
    finder_chart = salticam_finder_chart(
        general=_general_properties(fits_center),
        is_slot_mode=is_slot_mode,
    )
    check_finder(finder_chart)


@pytest.mark.parametrize("is_slot_mode", [False, True])
def test_rss_imaging_observation_annotation(
    is_slot_mode, fits_file, fits_center, check_finder, mock_salt_load_fits
):
    """Test the finder chart for RSS imaging observations."""
    finder_chart = rss_imaging_finder_chart(
        general=_general_properties(fits_center),
        is_slot_mode=is_slot_mode,
    )
    check_finder(finder_chart)


def test_rss_longslit_finder_chart(
    fits_file, fits_center, check_finder, mock_salt_load_fits
):
    """Test the finder chart for RSS longslit observations."""
    finder_chart = rss_longslit_finder_chart(
        general=_general_properties(fits_center),
        slit_width=4 * u.arcsec,
        slit_height=8 * u.arcmin,
    )
    check_finder(finder_chart)


def test_rss_mos_finder_chart(
    fits_file, fits_center, mos_mask_xml, check_finder, mock_salt_load_fits
):
    """Test the finder chart for RSS MOS observations."""
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
    finder_chart = rss_mos_finder_chart(
        general=_general_properties(fits_center), mos_mask=mos_mask
    )
    check_finder(finder_chart)


def test_rss_fabry_perot_finder_chart(
    fits_file, fits_center, check_finder, mock_salt_load_fits
):
    """Test the finder chart for RSS Fabry-Pérot observations."""
    finder_chart = rss_fabry_perot_finder_chart(
        general=_general_properties(fits_center),
    )
    check_finder(finder_chart)


def test_hrs_finder_chart(fits_file, fits_center, check_finder, mock_salt_load_fits):
    """Test the finder chart for an HRS observation."""
    finder_chart = hrs_finder_chart(_general_properties(fits_center))
    check_finder(finder_chart)
