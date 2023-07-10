from datetime import datetime, timedelta, timezone
from unittest import mock

import imephu
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from imephu.service.horizons import HorizonsService
from imephu.utils import Ephemeris

_MOCK_QUERY_RESULT_WITH_MAGNITUDE = Table(
    {
        "datetime_str": ["2022-Feb-08 13:27:14", "2022-Feb-08 13:32:14"],
        "RA": ["127", "127.01"],
        "DEC": ["-35", "-34.9"],
        "V": ["18", "18.1"],
    }
)


_MOCK_QUERY_RESULT_WITHOUT_MAGNITUDE = Table(
    {
        "datetime_str": ["2022-Feb-08 13:27:14", "2022-Feb-08 13:32:14"],
        "RA": ["127", "127.01"],
        "DEC": ["-35", "-34.9"],
    }
)


def test_horizons_constructor_requires_timezone_aware_start_time():
    """Test that the constructor requires a timezone-aware start time."""
    with pytest.raises(ValueError) as excinfo:
        HorizonsService(
            object_id="Ceres",
            location="B31",
            start=datetime(2022, 2, 8, 0, 0, 0, 0),
            end=datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
        )
    assert "start" in str(excinfo.value)


def test_horizons_constructor_requires_timezone_aware_end_time():
    """Test that the constructor requires a timezone-aware end time."""
    with pytest.raises(ValueError) as excinfo:
        HorizonsService(
            object_id="Ceres",
            location="B31",
            start=datetime(2022, 2, 8, 0, 0, 0, 0, tzinfo=timezone.utc),
            end=datetime(2022, 2, 8, 1, 0, 0, 0),
        )
    assert "end" in str(excinfo.value)


@pytest.mark.parametrize(
    "start",
    [
        datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
        datetime(2022, 2, 8, 1, 0, 1, 0, tzinfo=timezone.utc),
    ],
)
def test_horizons_start_and_end_time_must_be_consistent(start):
    """Test that the constructor start time must be earlier than the end time."""
    with pytest.raises(ValueError) as excinfo:
        HorizonsService(
            object_id="Ceres",
            location="B31",
            start=start,
            end=datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
        )
    assert "earlier" in str(excinfo.value)


def test_horizons_stepsize_must_be_at_least_five_minutes():
    """Test that the time between ephemeris epochs must be at least five minutes."""
    with pytest.raises(ValueError) as excinfo:
        HorizonsService(
            object_id="Ceres",
            location="B31",
            start=datetime(2022, 2, 8, 0, 0, 0, 0, tzinfo=timezone.utc),
            end=datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
            stepsize=4.999 * u.minute,
        )
    assert "time between" in str(excinfo.value)


def test_horizons_ephemerides_start_time_must_be_timezone_aware():
    """Test that the ephemerides method requires a timezone-aware start time."""
    horizons = HorizonsService(
        object_id="Ceres",
        location="B31",
        start=datetime(2022, 2, 8, 0, 0, 0, 0, tzinfo=timezone.utc),
        end=datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
    )
    with pytest.raises(ValueError) as excinfo:
        horizons.ephemerides(start=datetime(2022, 2, 8, 0, 0, 0, 0))
    assert "start" in str(excinfo.value)


def test_horizons_ephemerides_end_time_must_be_timezone_aware():
    """Test that the ephemerides method requires a timezone-aware end time."""
    horizons = HorizonsService(
        object_id="Ceres",
        location="B31",
        start=datetime(2022, 2, 8, 0, 0, 0, 0, tzinfo=timezone.utc),
        end=datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
    )
    with pytest.raises(ValueError) as excinfo:
        horizons.ephemerides(end=datetime(2022, 2, 8, 1, 0, 0, 0))
    assert "end" in str(excinfo.value)


@pytest.mark.parametrize(
    "start",
    [
        datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
        datetime(2022, 2, 8, 1, 0, 1, 0, tzinfo=timezone.utc),
    ],
)
def test_horizons_ephemerides_start_and_end_time_must_be_consistent(start):
    """Test that the ephemerides start time must be earlier than the end time."""
    horizons = HorizonsService(
        object_id="Ceres",
        location="B31",
        start=datetime(2022, 2, 8, 0, 0, 0, 0, tzinfo=timezone.utc),
        end=datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
    )
    with pytest.raises(ValueError) as excinfo:
        horizons.ephemerides(
            start=start, end=datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc)
        )
    assert "earlier" in str(excinfo.value)


def test_horizons_ephemerides_interval_must_be_subinterval():
    """Test that the constructor and ephemerides interval are consistent."""
    horizons = HorizonsService(
        object_id="Ceres",
        location="B31",
        start=datetime(2022, 2, 8, 0, 0, 0, 0, tzinfo=timezone.utc),
        end=datetime(2022, 2, 8, 1, 0, 0, 0, tzinfo=timezone.utc),
    )
    with pytest.raises(ValueError) as excinfo:
        horizons.ephemerides(
            start=datetime(2022, 2, 7, 23, 59, 59, 0, tzinfo=timezone.utc)
        )
    assert "start" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        horizons.ephemerides(end=datetime(2022, 2, 8, 1, 0, 1, 0, tzinfo=timezone.utc))
    assert "end" in str(excinfo.value)


def test_horizons_parses_query_results_with_magnitudes():
    """Test that the query result from JPL Horizons is parsed correctly."""
    with mock.patch.object(imephu.service.horizons, "Horizons") as MockHorizons:
        MockHorizons.return_value.ephemerides.return_value = (
            _MOCK_QUERY_RESULT_WITH_MAGNITUDE
        )
        start = datetime(2022, 2, 8, 13, 27, 14, 0, tzinfo=timezone.utc)
        end = datetime(2022, 2, 8, 13, 32, 14, 0, tzinfo=timezone.utc)
        horizons = HorizonsService(
            object_id="SomeAsteroid", location="B31", start=start, end=end
        )
        ephemerides = horizons.ephemerides(start, end)

        assert ephemerides[0].epoch == datetime(
            2022, 2, 8, 13, 27, 14, 0, tzinfo=timezone.utc
        )
        assert ephemerides[1].epoch == datetime(
            2022, 2, 8, 13, 32, 14, 0, tzinfo=timezone.utc
        )

        assert ephemerides[0].position.ra.to_value(u.deg) == pytest.approx(127)
        assert ephemerides[1].position.ra.to_value(u.deg) == pytest.approx(127.01)

        assert ephemerides[0].position.dec.to_value(u.deg) == pytest.approx(-35)
        assert ephemerides[1].position.dec.to_value(u.deg) == pytest.approx(-34.9)

        assert ephemerides[0].magnitude_range.min_magnitude == pytest.approx(18)
        assert ephemerides[0].magnitude_range.max_magnitude == pytest.approx(18)
        assert ephemerides[0].magnitude_range.bandpass == "V"

        assert ephemerides[1].magnitude_range.min_magnitude == pytest.approx(18.1)
        assert ephemerides[1].magnitude_range.max_magnitude == pytest.approx(18.1)
        assert ephemerides[1].magnitude_range.bandpass == "V"


def test_horizons_parses_query_results_without_magnitudes():
    """Test that the query result from JPL Horizons is parsed correctly."""
    with mock.patch.object(imephu.service.horizons, "Horizons") as MockHorizons:
        MockHorizons.return_value.ephemerides.return_value = (
            _MOCK_QUERY_RESULT_WITHOUT_MAGNITUDE
        )
        start = datetime(2022, 2, 8, 13, 27, 14, 0, tzinfo=timezone.utc)
        end = datetime(2022, 2, 8, 13, 32, 14, 0, tzinfo=timezone.utc)
        horizons = HorizonsService(
            object_id="SomeAsteroid", location="B31", start=start, end=end
        )
        ephemerides = horizons.ephemerides(start, end)

        assert ephemerides[0].epoch == datetime(
            2022, 2, 8, 13, 27, 14, 0, tzinfo=timezone.utc
        )
        assert ephemerides[1].epoch == datetime(
            2022, 2, 8, 13, 32, 14, 0, tzinfo=timezone.utc
        )

        assert ephemerides[0].position.ra.to_value(u.deg) == pytest.approx(127)
        assert ephemerides[1].position.ra.to_value(u.deg) == pytest.approx(127.01)

        assert ephemerides[0].position.dec.to_value(u.deg) == pytest.approx(-35)
        assert ephemerides[1].position.dec.to_value(u.deg) == pytest.approx(-34.9)

        assert ephemerides[0].magnitude_range is None
        assert ephemerides[1].magnitude_range is None


def test_jpl_horizons_is_queried_only_once():
    """Test that the JPL Horizons service is queried only once."""
    with mock.patch.object(imephu.service.horizons, "Horizons") as MockHorizons:
        MockHorizons.return_value.ephemerides.return_value = (
            _MOCK_QUERY_RESULT_WITHOUT_MAGNITUDE
        )
        start = datetime(2022, 2, 8, 13, 27, 14, 0, tzinfo=timezone.utc)
        end = datetime(2022, 2, 8, 13, 32, 14, 0, tzinfo=timezone.utc)
        horizons = HorizonsService(
            object_id="Ceres", location="B31", start=start, end=end
        )
        horizons.ephemerides(start, end)
        horizons.ephemerides(start, end)
        MockHorizons.assert_called_once()


def _ephemeris(time_diff):
    return Ephemeris(
        epoch=datetime(2022, 2, 8, 0, 0, 0, 0, tzinfo=timezone.utc)
        + timedelta(minutes=time_diff),
        position=SkyCoord(ra=0 * u.deg, dec=0 * u.deg),
        magnitude_range=None,
    )


@pytest.mark.parametrize(
    "start_time_diff,end_time_diff,expected",
    [
        (0, 10, [_ephemeris(time_diff) for time_diff in range(0, 11)]),
        (0, 5, [_ephemeris(time_diff) for time_diff in range(0, 6)]),
        (8, 10, [_ephemeris(time_diff) for time_diff in range(8, 11)]),
        (2, 6, [_ephemeris(time_diff) for time_diff in range(2, 7)]),
        (0.3, 7.6, [_ephemeris(time_diff) for time_diff in range(0, 9)]),
        (5.3, 9.1, [_ephemeris(time_diff) for time_diff in range(5, 11)]),
        (2.9, 7.7, [_ephemeris(time_diff) for time_diff in range(2, 9)]),
    ],
)
def test_ephemerides_returns_the_correct_ephemerides(
    start_time_diff, end_time_diff, expected
):
    """Test that the ephemerides method returns the correct ephemerides."""
    all_ephemerides = [_ephemeris(time_diff) for time_diff in range(0, 11)]
    horizons = HorizonsService(
        object_id="Ceres",
        location="B31",
        start=all_ephemerides[0].epoch,
        end=all_ephemerides[-1].epoch,
    )
    horizons._ephemerides = all_ephemerides
    start = _ephemeris(start_time_diff).epoch
    end = _ephemeris(end_time_diff).epoch
    ephemerides = horizons.ephemerides(start=start, end=end)
    assert ephemerides == expected
