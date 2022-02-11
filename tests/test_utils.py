from contextlib import contextmanager
from datetime import datetime, timezone

import pytest
from astropy.coordinates import SkyCoord

from imephu.salt.finder_chart import MagnitudeRange
from imephu.utils import Ephemeris


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