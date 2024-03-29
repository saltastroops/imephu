from datetime import datetime, timedelta, timezone

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from imephu.geometry import translate
from imephu.salt.finder_chart import (
    GeneralProperties,
    SaltFinderChart,
    Target,
    hrs_finder_chart,
    moving_target_finder_charts,
    nir_finder_chart,
    rss_fabry_perot_finder_chart,
    rss_imaging_finder_chart,
    rss_longslit_finder_chart,
    rss_mos_finder_chart,
    salticam_finder_chart,
)
from imephu.salt.utils import MosMask, MosMaskSlit
from imephu.utils import Ephemeris, MagnitudeRange, SkyCoordRate

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
    is_slot_mode, fits_file, fits_center, check_finder, mock_from_survey
):
    """Test the finder chart for Salticam observations."""
    finder_chart = salticam_finder_chart(
        fits=fits_file,
        general=_general_properties(fits_center),
        is_slot_mode=is_slot_mode,
    )
    check_finder(finder_chart)


@pytest.mark.parametrize("is_slot_mode", [False, True])
def test_rss_imaging_observation_annotation(
    is_slot_mode, fits_file, fits_center, check_finder, mock_from_survey
):
    """Test the finder chart for RSS imaging observations."""
    finder_chart = rss_imaging_finder_chart(
        fits=fits_file,
        general=_general_properties(fits_center),
        is_slot_mode=is_slot_mode,
    )
    check_finder(finder_chart)


@pytest.mark.parametrize("include_reference_star", [True, False])
def test_rss_longslit_finder_chart(
    include_reference_star, fits_file, fits_center, check_finder, mock_from_survey
):
    """Test the finder chart for RSS longslit observations."""
    if include_reference_star:
        reference_star = translate(fits_center, (2, -1) * u.arcmin)
    else:
        reference_star = None
    finder_chart = rss_longslit_finder_chart(
        fits=fits_file,
        general=_general_properties(fits_center),
        slit_width=4 * u.arcsec,
        slit_height=8 * u.arcmin,
        reference_star=reference_star,
    )
    check_finder(finder_chart)


def test_rss_mos_finder_chart(
    fits_file, fits_center, mos_mask_xml, check_finder, mock_from_survey
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
        fits=fits_file, general=_general_properties(fits_center), mos_mask=mos_mask
    )
    check_finder(finder_chart)


def test_rss_fabry_perot_finder_chart(
    fits_file, fits_center, check_finder, mock_from_survey
):
    """Test the finder chart for RSS Fabry-Pérot observations."""
    finder_chart = rss_fabry_perot_finder_chart(
        fits=fits_file,
        general=_general_properties(fits_center),
    )
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "position_angle, target_offset",
    [
        (0 * u.deg, None),
        (30 * u.deg, (0.5, 1) * u.arcmin),
        (-90 * u.deg, (-2, -1) * u.arcmin),
    ],
)
def test_nir_finder_chart(
    position_angle,
    target_offset,
    fits_file,
    fits_center,
    check_finder,
    mock_from_survey,
):
    """Test the finder chart for an NIR observation."""
    general = _general_properties(fits_center)
    if target_offset is not None:
        reference_star = fits_center
        general.target.position = translate(fits_center, target_offset)
    else:
        reference_star = None
    general.position_angle = position_angle
    finder_chart = nir_finder_chart(
        fits=fits_file,
        general=general,
        reference_star=reference_star,
        bundle_separation=1 * u.arcmin,
    )
    check_finder(finder_chart)


def test_hrs_finder_chart(fits_file, fits_center, check_finder, mock_from_survey):
    """Test the finder chart for an HRS observation."""
    finder_chart = hrs_finder_chart(
        fits=fits_file, general=_general_properties(fits_center)
    )
    check_finder(finder_chart)


@pytest.mark.parametrize("which_finder_chart", [0, 1])
def test_moving_target_finder_charts(
    which_finder_chart, fits_center, mock_from_survey, check_finder
):
    """Test that the generated finder charts for a moving target are correct."""
    t = datetime(2022, 2, 18, 0, 0, 0, tzinfo=timezone.utc)
    hour = timedelta(hours=1)
    ephemerides = [
        Ephemeris(
            epoch=t,
            position=SkyCoord(ra="0h40m30s", dec=-60 * u.deg),
            position_rate=SkyCoordRate(
                ra=1 * u.arcsec / u.hour, dec=1 * u.arcsec / u.hour
            ),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + hour,
            position=SkyCoord(ra="0h40m00s", dec=-60 * u.deg),
            position_rate=SkyCoordRate(
                ra=1 * u.arcsec / u.hour, dec=1 * u.arcsec / u.hour
            ),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + 2 * hour,
            position=SkyCoord(ra="0h39m30s", dec=-60 * u.deg),
            position_rate=SkyCoordRate(
                ra=1 * u.arcsec / u.hour, dec=1 * u.arcsec / u.hour
            ),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + 3 * hour,
            position=SkyCoord(ra="0h39m00s", dec=-60 * u.deg),
            position_rate=SkyCoordRate(
                ra=1 * u.arcsec / u.hour, dec=1 * u.arcsec / u.hour
            ),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + 4 * hour,
            position=SkyCoord(ra="0h38m30s", dec=-60 * u.deg),
            position_rate=SkyCoordRate(
                ra=1 * u.arcsec / u.hour, dec=1 * u.arcsec / u.hour
            ),
            magnitude_range=None,
        ),
    ]
    start = t + 0.5 * hour
    end = t + 3.5 * hour
    g = moving_target_finder_charts(
        general=_general_properties(fits_center),
        start=start,
        end=end,
        ephemerides=ephemerides,
        survey="POSS2/UKSTU Red",
    )
    counter = 0
    for finder_chart, _ in g:  # noqa: B007
        # The regression fixture only creates a single file per test run, it seems. So
        # to get all finder charts and ensure a deterministic result, we have to make
        # sure it is used only once.
        if counter == which_finder_chart:
            check_finder(finder_chart)
        counter += 1

    assert counter == 2


def test_target_is_added_to_metadata(fits_file, fits_center):
    """Test that the target coordinates are added to the pdf metadata attachment."""
    np.random.seed(0)
    general = _general_properties(fits_center)
    finder_chart = SaltFinderChart(fits_file, general.target)

    metadata = finder_chart.metadata
    assert metadata["right_ascension"] == "10.0 deg"
    assert metadata["declination"] == "-60.0 deg"
