import io
import pathlib
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
from typing import Any, cast
from unittest import mock

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

import imephu
import imephu.service.survey
from imephu.finder_chart import FinderChart
from imephu.utils import Ephemeris, MagnitudeRange


class _FakeFinderChart:
    def __init__(self, expected_value: Any):
        self.expected_value = expected_value
        self._metadata = {}

    def add_metadata(self, key: str, value: Any) -> None:
        self._metadata[key] = value

    @property
    def metadata(self):
        return self._metadata.copy()


@contextmanager
def _does_not_raise():
    yield


def _ephemeris(epoch, ra=0 * u.deg, dec=0 * u.deg):
    return Ephemeris(
        epoch=epoch,
        position=SkyCoord(ra=ra, dec=dec),
        magnitude_range=MagnitudeRange(
            bandpass="V", min_magnitude=13, max_magnitude=13
        ),
    )


def test_finder_chart_is_generated_from_path(check_finder):
    """Test that finder charts can be generated with a `~pathlib.Path object."""
    path = pathlib.Path(__file__).parent / "data" / "ra10_dec-60.fits"
    finder_chart = FinderChart(path)
    check_finder(finder_chart)


def test_finder_chart_is_generated_from_string(check_finder):
    """Test that finder charts can be generated with a string containing the path."""
    name = str(pathlib.Path(__file__).parent / "data" / "ra10_dec-60.fits")
    finder_chart = FinderChart(name)
    check_finder(finder_chart)


def test_finder_chart_is_generated_from_stream(check_finder):
    """Test that finder charts can be generated with a binary stream."""
    path = pathlib.Path(__file__).parent / "data" / "ra10_dec-60.fits"
    with open(path, "rb") as f:
        finder_chart = FinderChart(f)
    check_finder(finder_chart)


# Formats other than jpg or png may produce different files for different runs,
# so that they cannot be tested with pytest-regressions.
@pytest.mark.parametrize("format", ["jpg", "png"])
def test_finder_chart_export_formats(format, check_image, fits_file):
    """Test that finder charts can be exported to different file formats."""
    np.random.seed(0)
    finder_chart = FinderChart(fits_file)
    contents = io.BytesIO()
    finder_chart.save(contents, format=format)
    check_image(contents)


def test_finder_chart_from_survey_returns_finder_chart(
    fits_file, fits_center, check_finder
):
    """Test that the from_survey method returns a finder chart."""
    with mock.patch.object(
        imephu.finder_chart, "load_fits", autospec=True
    ) as mock_load_fits:
        mock_load_fits.return_value = open(fits_file, "rb")
        finder_chart = FinderChart.from_survey(
            "POSS2/UKSTU Red", fits_center, 10 * u.arcmin
        )
        check_finder(finder_chart)


def test_for_time_interval_start_must_be_timezone_aware():
    """Test that the start time must be timezone-aware."""
    start = datetime(2022, 2, 16, 15, 0, 0, 0)
    end = datetime(2022, 2, 2, 10, 0, 0, 0, tzinfo=timezone.utc)
    ephemerides = [
        _ephemeris(datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc)),
        _ephemeris(datetime(2022, 2, 2, 12, 0, 0, 0, tzinfo=timezone.utc)),
    ]
    with pytest.raises(ValueError) as excinfo:
        g = FinderChart.for_time_interval(
            start, end, ephemerides, 5 * u.arcmin, cast(Any, lambda x: 42)
        )
        next(g)
    assert "start" in str(excinfo.value) and "timezone" in str(excinfo.value)


def test_for_time_interval_end_must_be_timezone_aware():
    """Test that the end time must be timezone-aware."""
    start = datetime(2022, 2, 16, 15, 0, 0, 0, tzinfo=timezone.utc)
    end = datetime(2022, 2, 2, 10, 0, 0, 0)
    ephemerides = [
        _ephemeris(datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc)),
        _ephemeris(datetime(2022, 2, 2, 12, 0, 0, 0, tzinfo=timezone.utc)),
    ]
    with pytest.raises(ValueError) as excinfo:
        g = FinderChart.for_time_interval(
            start, end, ephemerides, 5 * u.arcmin, cast(Any, lambda x: 42)
        )
        next(g)
    assert "end" in str(excinfo.value) and "timezone" in str(excinfo.value)


@pytest.mark.parametrize(
    "start,end,expectation",
    [
        (
            datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc),
            datetime(2022, 2, 16, 12, 0, 1, 0, tzinfo=timezone.utc),
            _does_not_raise(),
        ),
        (
            datetime(2022, 2, 16, 12, 0, 1, 0, tzinfo=timezone.utc),
            datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc),
            pytest.raises(ValueError),
        ),
        (
            datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc),
            datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc),
            pytest.raises(ValueError),
        ),
    ],
)
def test_for_time_interval_start_must_be_earlier_than_end(start, end, expectation):
    """Test that the start time must be earlier than the end time."""
    ephemerides = [
        _ephemeris(datetime(2022, 2, 15, 0, 0, 0, 0, tzinfo=timezone.utc)),
        _ephemeris(datetime(2022, 2, 20, 12, 0, 0, 0, tzinfo=timezone.utc)),
    ]
    with expectation:
        g = FinderChart.for_time_interval(
            start,
            end,
            ephemerides,
            5 * u.arcmin,
            cast(Any, lambda x: _FakeFinderChart(42)),
        )
        next(g)


@pytest.mark.parametrize(
    "start,end,expectation",
    [
        (
            datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc),
            datetime(2022, 2, 20, 12, 0, 0, 0, tzinfo=timezone.utc),
            _does_not_raise(),
        ),
        (
            datetime(2022, 2, 16, 11, 59, 59, 0, tzinfo=timezone.utc),
            datetime(2022, 2, 20, 12, 0, 0, 0, tzinfo=timezone.utc),
            pytest.raises(ValueError),
        ),
        (
            datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc),
            datetime(2022, 2, 20, 12, 0, 1, 0, tzinfo=timezone.utc),
            pytest.raises(ValueError),
        ),
    ],
)
def test_for_time_interval_time_intervals_must_be_covered(start, end, expectation):
    """Test that the time interval is covered by the ephemerides."""
    ephemerides = [
        _ephemeris(datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc)),
        _ephemeris(datetime(2022, 2, 20, 12, 0, 0, 0, tzinfo=timezone.utc)),
    ]
    with expectation:
        g = FinderChart.for_time_interval(
            start,
            end,
            ephemerides,
            5 * u.arcmin,
            cast(Any, lambda x: _FakeFinderChart(42)),
        )
        next(g)


@pytest.mark.parametrize(
    "max_track_length,expectation",
    [
        (-1 * u.arcmin, pytest.raises(ValueError)),
        (0 * u.arcmin, pytest.raises(ValueError)),
        (1 * u.arcmin, _does_not_raise()),
    ],
)
def test_for_time_interval_max_track_length_must_be_positive(
    max_track_length, expectation
):
    """Test that the maximum track length must be positive."""
    start = datetime(2022, 2, 16, 15, 0, 0, 0, tzinfo=timezone.utc)
    end = datetime(2022, 2, 20, 10, 0, 0, 0, tzinfo=timezone.utc)
    ephemerides = [
        _ephemeris(datetime(2022, 2, 16, 12, 0, 0, 0, tzinfo=timezone.utc)),
        _ephemeris(datetime(2022, 2, 20, 12, 0, 0, 0, tzinfo=timezone.utc)),
    ]
    with expectation:
        g = FinderChart.for_time_interval(
            start,
            end,
            ephemerides,
            max_track_length,
            cast(Any, lambda x: _FakeFinderChart(42)),
        )
        next(g)


@pytest.mark.parametrize(
    "ras,decs",
    [
        # For a declination of 60 degrees the actual angle on the sky is half the
        # difference in right ascension.
        (
            [0 * u.deg, 2.1 * u.deg, 4 * u.deg, 6 * u.deg],
            [60 * u.deg, 60 * u.deg, 60 * u.deg, 60 * u.deg],
        ),
        (
            [0 * u.deg, 2 * u.deg, 4.1 * u.deg, 6 * u.deg],
            [60 * u.deg, 60 * u.deg, 60 * u.deg, 60 * u.deg],
        ),
        (
            [0 * u.deg, 2 * u.deg, 4 * u.deg, 6.3 * u.deg],
            [60 * u.deg, 60 * u.deg, 60 * u.deg, 60 * u.deg],
        ),
        (
            [0 * u.deg, 0 * u.deg, 0 * u.deg, 0 * u.deg],
            [0 * u.deg, 1.1 * u.deg, 2 * u.deg, 3 * u.deg],
        ),
        (
            [0 * u.deg, 0 * u.deg, 0 * u.deg, 0 * u.deg],
            [0 * u.deg, 1 * u.deg, 2.1 * u.deg, 3 * u.deg],
        ),
        (
            [0 * u.deg, 0 * u.deg, 0 * u.deg, 0 * u.deg],
            [0 * u.deg, 1 * u.deg, 2 * u.deg, 3.1 * u.deg],
        ),
    ],
)
def test_for_time_interval_requires_sufficiently_close_ephemeris_positions(ras, decs):
    """Test that the maximum track length must not be exceeded on a finder chart."""
    t = datetime(2022, 2, 17, 9, 0, 0, 0, tzinfo=timezone.utc)
    ephemerides = [
        _ephemeris(t + timedelta(hours=i), ra=ras[i], dec=decs[i])
        for i, _ in enumerate(ras)
    ]
    start = ephemerides[0].epoch
    end = ephemerides[-1].epoch
    with pytest.raises(ValueError):
        g = FinderChart.for_time_interval(
            start, end, ephemerides, 1 * u.deg, cast(Any, lambda x: 42)
        )
        next(g)


@pytest.mark.parametrize(
    "positions,expected_indices",
    [
        # For a declination of 60 degrees the actual angle on the sky is half the
        # difference in right ascension.
        (
            (
                (0 * u.deg, 60 * u.deg),
                (0.6 * u.deg, 60 * u.deg),
                (1.8 * u.deg, 60 * u.deg),
                (2.1 * u.deg, 60 * u.deg),
                (3.9 * u.deg, 60 * u.deg),
                (4 * u.deg, 60 * u.deg),
            ),
            ((0, 1, 2), (2, 3), (3, 4, 5)),
        ),
        (
            (
                (0 * u.deg, 0 * u.deg),
                (0 * u.deg, 0.8 * u.deg),
                (0 * u.deg, 1.2 * u.deg),
                (0 * u.deg, 2.1 * u.deg),
            ),
            ((0, 1), (1, 2), (2, 3)),
        ),
    ],
)
def test_for_time_interval_creates_correct_finder_charts(positions, expected_indices):
    """Test that the correct finder charts are created."""
    t = datetime(2022, 2, 17, 0, 0, 0, 0, tzinfo=timezone.utc)
    ephemerides = [
        _ephemeris(
            epoch=t + timedelta(hours=i), ra=positions[i][0], dec=positions[i][1]
        )
        for i, _ in enumerate(positions)
    ]
    expected_ephemerides = []
    expected_valid_for = []
    for index_group in expected_indices:
        expected_ephemerides.append([ephemerides[i] for i in index_group])
        expected_valid_for.append(
            (ephemerides[index_group[0]].epoch, ephemerides[index_group[-1]].epoch)
        )

    start = t
    end = t + timedelta(hours=len(positions) - 1)

    # Fake function for creating a finder chart.
    def fake_create_finder_chart(e):
        return _FakeFinderChart(e)

    g = FinderChart.for_time_interval(
        start, end, ephemerides, 1 * u.deg, cast(Any, fake_create_finder_chart)
    )
    g = list(g)
    assert [cast(Any, f[0]).expected_value for f in g] == expected_ephemerides
    assert [f[0].metadata["valid_for"] for f in g] == expected_valid_for
