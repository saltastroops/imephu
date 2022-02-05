from dataclasses import dataclass
from datetime import datetime

from astropy.coordinates import SkyCoord


@dataclass
class Ephemeris:
    """An ephemeris with an epoch and a position.

    Parameters
    ----------
    epoch: `~datetime.datetime`
        The epoch, i.e. the datetime for which the position is given. The epoch must be
        a timezone-aware datetime.
    position: `~astropy.coordinates.SkyCoord`
        The position, in right ascension and declination.
    """

    epoch: datetime

    position: SkyCoord

    def __post_init__(self) -> None:
        """Check that the epoch is timezone-aware."""
        if self.epoch.tzinfo is None or self.epoch.tzinfo.utcoffset(None) is None:
            raise ValueError("The epoch must be a timezone-aware datetime object.")
