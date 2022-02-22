from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Optional

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.utils import (
    Ephemeris,
    MagnitudeRange,
    ephemerides_magnitude_range,
    mid_position,
)


@contextmanager
def _does_not_raise():
    yield


@pytest.mark.parametrize(
    "epoch,expectation",
    [
        (datetime(2022, 2, 5), pytest.raises(ValueError)),
        (datetime(2022, 2, 5, tzinfo=timezone.utc), _does_not_raise()),
    ],
)
def test_ephemeris_epoch_must_be_timezone_aware(epoch, expectation):
    """Test that ephemeris epochs must be timezone-aware."""
    with expectation:
        Ephemeris(
            epoch=epoch,
            position=SkyCoord(ra="0h00m00s", dec="0d00m00s"),
            magnitude_range=MagnitudeRange(
                min_magnitude=18, max_magnitude=18, bandpass="V"
            ),
        )


@pytest.mark.parametrize(
    "start,end,expected",
    [
        (
            SkyCoord(ra=10 * u.deg, dec=0 * u.deg),
            SkyCoord(ra=15 * u.deg, dec=0 * u.deg),
            SkyCoord(ra=12.5 * u.deg, dec=0 * u.deg),
        ),
        (
            SkyCoord(ra=130 * u.deg, dec=-20 * u.deg),
            SkyCoord(ra=130 * u.deg, dec=-40 * u.deg),
            SkyCoord(ra=130 * u.deg, dec=-30 * u.deg),
        ),
    ],
)
def test_mid_position_on_sky(start, end, expected):
    """Test that the mid_position function returns the correct position."""
    mp = mid_position(start=start, end=end)
    assert mp.ra.to_value(u.deg) == pytest.approx(expected.ra.to_value(u.deg))
    assert mp.dec.to_value(u.deg) == pytest.approx(expected.dec.to_value(u.deg))


@pytest.mark.parametrize(
    "bandpasses,expectation",
    [
        ((None, None), _does_not_raise()),
        ((None, "V"), _does_not_raise()),
        (("V",), _does_not_raise()),
        (("R", "R"), _does_not_raise()),
        (("V", "R", "V"), pytest.raises(ValueError)),
    ],
)
def test_ephemerides_magnitude_range_requires_consistent_bandpass(
    bandpasses, expectation
):
    """Test that the bandpass must be the same for every ephemeris."""
    ephemerides = []
    for bandpass in bandpasses:
        if bandpass:
            mr: Optional[MagnitudeRange] = MagnitudeRange(
                min_magnitude=17.1, max_magnitude=17.2, bandpass=bandpass
            )
        else:
            mr = None
        ephemerides.append(
            Ephemeris(
                epoch=datetime(2022, 2, 19, 0, 0, 0, 0, tzinfo=timezone.utc),
                position=SkyCoord(ra=0 * u.deg, dec=0 * u.deg),
                magnitude_range=mr,
            )
        )

    with expectation:
        ephemerides_magnitude_range(ephemerides)


@pytest.mark.parametrize(
    "min_magnitudes,max_magnitudes,expected",
    [
        (
            (16.9,),
            (16.9,),
            MagnitudeRange(min_magnitude=16.9, max_magnitude=16.9, bandpass="B"),
        ),
        (
            (16.9,),
            (17.0,),
            MagnitudeRange(min_magnitude=16.9, max_magnitude=17.0, bandpass="B"),
        ),
        (
            (12.5, None, 14.7, 11.3, 16.7, 11.5),
            (12.7, None, 15, 11.3, 16.8, 11.5),
            MagnitudeRange(min_magnitude=11.3, max_magnitude=16.8, bandpass="B"),
        ),
        ((None,), (None,), None),
        ((None, None), (None, None), None),
    ],
)
def test_ephemerides_magnitude_range(min_magnitudes, max_magnitudes, expected):
    """Test that the ephemerides_magnitude_range function returns the correct range."""
    ephemerides = []
    for i, _ in enumerate(min_magnitudes):
        if min_magnitudes[i] is not None:
            mr: Optional[MagnitudeRange] = MagnitudeRange(
                min_magnitude=min_magnitudes[i],
                max_magnitude=max_magnitudes[i],
                bandpass="B",
            )
        else:
            mr = None
        ephemerides.append(
            Ephemeris(
                epoch=datetime(2022, 2, 19, 0, 0, 0, 0, tzinfo=timezone.utc),
                position=SkyCoord(ra=0 * u.deg, dec=0 * u.deg),
                magnitude_range=mr,
            )
        )
    magnitude_range = ephemerides_magnitude_range(ephemerides)
    assert magnitude_range == expected
