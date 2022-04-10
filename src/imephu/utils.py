from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional

from astropy.coordinates import SkyCoord


@dataclass
class MagnitudeRange:
    """A magnitude range.

    Attributes
    ----------
    bandpass: `str`
        The bandpass for which the magnitudes are given.
    min_magnitude: `float`
        The minimum (brightest) magnitude.
    max_magnitude: `float`
        The maximum (faintest) magnitude.
    """

    bandpass: str
    min_magnitude: float
    max_magnitude: float


@dataclass
class Ephemeris:
    """An ephemeris with an epoch, position and magnitude range.

    Parameters
    ----------
    epoch: `~datetime.datetime`
        The epoch, i.e. the datetime for which the position is given. The epoch must be
        a timezone-aware datetime.
    position: `~astropy.coordinates.SkyCoord`
        The position, in right ascension and declination.
    magnitude_range: `~imephu.utils.MagnitudeRange`
        The magnitude range.
    """

    epoch: datetime

    position: SkyCoord

    magnitude_range: Optional[MagnitudeRange]

    def __post_init__(self) -> None:
        """Check that the epoch is timezone-aware."""
        if self.epoch.tzinfo is None or self.epoch.tzinfo.utcoffset(None) is None:
            raise ValueError("The epoch must be a timezone-aware datetime object.")


def mid_position(start: SkyCoord, end: SkyCoord) -> SkyCoord:
    """Return the mid position between a start and end position on the sky.

    The mid position is the mid position on the great arc between the start and end
    position.

    Taken from https://github.com/astropy/astropy/issues/5766.

    Parameters
    ----------
    start: `~astropy.coordinates.SkyCoords`
        Start position on the sky.
    end: `~astropy.coordinates.SkyCoords`
        End position on the sky.
    """
    pa = start.position_angle(end)
    separation = start.separation(end)
    return start.directional_offset_by(pa, separation / 2)


def ephemerides_magnitude_range(
    ephemerides: List[Ephemeris],
) -> Optional[MagnitudeRange]:
    """Return the magnitude range for a list of ephemerides.

    The minimum (maximum) magnitude is the minimum (maximum) magnitude for all
    ephemerides. If none of the ephemerides has a magnitude range, ``None`` is returned.

    Parameters
    ----------
    ephemerides: list of `~imephu.utils.Ephemeris`
        The list of ephemerides.

    Returns
    -------
    `~imephu.utils.MagnitudeRange`, optional
        The magnitude range for the list of ephemerides.
    """
    min_magnitude: Optional[float] = None
    max_magnitude: Optional[float] = None
    bandpass: Optional[str] = None
    for ephemeris in ephemerides:
        if ephemeris.magnitude_range:
            if (
                min_magnitude is None
                or ephemeris.magnitude_range.min_magnitude < min_magnitude
            ):
                min_magnitude = ephemeris.magnitude_range.min_magnitude
            if (
                max_magnitude is None
                or ephemeris.magnitude_range.max_magnitude > max_magnitude
            ):
                max_magnitude = ephemeris.magnitude_range.max_magnitude
            if bandpass is None:
                bandpass = ephemeris.magnitude_range.bandpass
            elif ephemeris.magnitude_range.bandpass != bandpass:
                raise ValueError("The bandpass must be the same for all ephemerides.")

    if min_magnitude is not None and max_magnitude is not None and bandpass is not None:
        return MagnitudeRange(
            min_magnitude=min_magnitude, max_magnitude=max_magnitude, bandpass=bandpass
        )
    else:
        return None
