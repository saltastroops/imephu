from dataclasses import dataclass
from datetime import datetime
from typing import Optional

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
