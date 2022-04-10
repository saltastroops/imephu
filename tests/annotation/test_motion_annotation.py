from datetime import datetime, timedelta, timezone

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.motion import motion_annotation
from imephu.finder_chart import FinderChart
from imephu.geometry import rotate
from imephu.utils import Ephemeris, MagnitudeRange


@pytest.mark.parametrize(
    "angle",
    [
        0 * u.deg,
        45 * u.deg,
        90 * u.deg,
        135 * u.deg,
        180 * u.deg,
        225 * u.deg,
        270 * u.deg,
        315 * u.deg,
    ],
)
def test_large_motion_annotation(angle, fits_file, fits_center, check_finder, legend):
    """Test the annotation for large motions."""
    finder_chart = FinderChart(fits_file)
    initial_position = rotate(
        v=SkyCoord(ra="00h40m10s", dec="-60d00m00s"),
        pivot=fits_center,
        angle=angle,
        wcs=finder_chart.wcs,
    )
    initial_ephemeris = Ephemeris(
        epoch=datetime(2022, 2, 5, 22, 0, 0, 0, tzinfo=timezone.utc),
        position=initial_position,
        magnitude_range=MagnitudeRange(
            min_magnitude=18, max_magnitude=18, bandpass="V"
        ),
    )
    middle_ephemeris = Ephemeris(
        epoch=datetime(2022, 2, 6, 0, 0, 0, 0, tzinfo=timezone.utc),
        position=fits_center,
        magnitude_range=MagnitudeRange(
            min_magnitude=18, max_magnitude=18, bandpass="V"
        ),
    )
    final_position = rotate(
        v=SkyCoord(ra="00h39m50s", dec="-60d00m00s"),
        pivot=fits_center,
        angle=angle,
        wcs=finder_chart.wcs,
    )
    sast = timezone(timedelta(hours=2), "SAST")
    # note the different timezone for the final ephemeris
    final_ephemeris = Ephemeris(
        epoch=datetime(2022, 2, 6, 4, 0, 0, 0, tzinfo=sast),
        position=final_position,
        magnitude_range=MagnitudeRange(
            min_magnitude=18, max_magnitude=18, bandpass="V"
        ),
    )
    # ephemerides need not be sorted by epoch
    ephemerides = [middle_ephemeris, initial_ephemeris, final_ephemeris]

    motion_annotation_ = motion_annotation(
        ephemerides=ephemerides, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(motion_annotation_)
    finder_chart.add_annotation(legend(f"Angle: {angle}", finder_chart.wcs))
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "angle",
    [
        0 * u.deg,
        45 * u.deg,
        90 * u.deg,
        135 * u.deg,
        180 * u.deg,
        225 * u.deg,
        270 * u.deg,
        315 * u.deg,
    ],
)
def test_small_motion_annotation(angle, fits_file, fits_center, check_finder, legend):
    """Test the annotation for small motions."""
    finder_chart = FinderChart(fits_file)
    initial_position = rotate(
        v=SkyCoord(ra="00h40m00.5s", dec="-60d00m00s"),
        pivot=fits_center,
        angle=angle,
        wcs=finder_chart.wcs,
    )
    initial_ephemeris = Ephemeris(
        epoch=datetime(2022, 2, 5, 22, 0, 0, 0, tzinfo=timezone.utc),
        position=initial_position,
        magnitude_range=MagnitudeRange(
            min_magnitude=18, max_magnitude=18, bandpass="V"
        ),
    )
    middle_ephemeris = Ephemeris(
        epoch=datetime(2022, 2, 6, 0, 0, 0, 0, tzinfo=timezone.utc),
        position=fits_center,
        magnitude_range=MagnitudeRange(
            min_magnitude=18, max_magnitude=18, bandpass="V"
        ),
    )
    final_position = rotate(
        v=SkyCoord(ra="00h39m59.5s", dec="-60d00m00s"),
        pivot=fits_center,
        angle=angle,
        wcs=finder_chart.wcs,
    )
    sast = timezone(timedelta(hours=2), "SAST")
    # note the different timezone for the final ephemeris
    final_ephemeris = Ephemeris(
        epoch=datetime(2022, 2, 6, 4, 0, 0, 0, tzinfo=sast),
        position=final_position,
        magnitude_range=MagnitudeRange(
            min_magnitude=18, max_magnitude=18, bandpass="V"
        ),
    )
    # ephemerides need not be sorted by epoch
    ephemerides = [middle_ephemeris, initial_ephemeris, final_ephemeris]

    motion_annotation_ = motion_annotation(
        ephemerides=ephemerides, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(motion_annotation_)
    finder_chart.add_annotation(legend(f"Angle: {angle}", finder_chart.wcs))
    check_finder(finder_chart)
