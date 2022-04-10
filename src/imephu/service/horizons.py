import bisect
from datetime import datetime, timedelta, timezone
from typing import List, Optional

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astroquery.jplhorizons import Horizons
from dateutil.parser import parse

from imephu.utils import Ephemeris, MagnitudeRange


class HorizonsService:
    """A service for querying ephemerides from JPL Horizons.

    The service is initialised by calling the constructor with a start and end time,
    and ephemerides can be queried with the `ephemerides` method, which takes a start
    and end time as well. The first time this method is called, JPL Horizons is queried
    for the ephemerides for the interval defined by the _constructor's_ start and end
    time, and the returned list of ephemerides is stored. This list is used to find the
    subset of ephemerides covering the interval defined by the start and end time passed
    to the method. Subsequent calls to the `ephemerides` method use this list as well.

    Effectively this means irrespective how often you call the `ephemerides` method,
    only one request is made to the JPL Horizons service. It also means that the time
    interval passed to the `ephemerides` method must be fully contained within the
    interval passed to the constructor.

    The magnitude range for the ephemerides is given for the V band, and the minimum
    and maximum magnitude are always equal. If the Horizons data contain no V band
    magnitude, no magnitude range is included in the ephemeris.

    The time between ephemeris values is rounded to the nearest minute.

    .. note::
       To ensure full coverage of the interval from ``start`` to ``end``, the actual
       start and end time passed when querying JPL Horizons are slightly different from
       ``start`` and ``end``.

    Parameters
    ----------
    object_id: `str`
        The Horizons identifier for this object, such as ``Ceres``. The exact string
        will be used when querying the ephemerides from JPL Horizons.
    location: `str`
        The location as an `MPC observatory
        code<https://minorplanetcenter.net//iau/lists/ObsCodesF.html>`_ For example, the
        location for the Southern African Large Telescope is B31.
    start: `~datetime.datetime`
        The start of the interval for which ephemerides should be queried. This must be
        a timezone-aware `~datetime.datetime` object.
    end: `~datetime.datetime`
        The end of the interval for which ephemerides should be queried. This must be
        a timezone-aware `~datetime.datetime` object.
    stepsize: `~astropy.units.Quantity`, default: 5 minutes
        The time between ephemeris values. This must be at least 5 minutes. The time is
        rounded to the nearest integer when making the query.
    """

    def __init__(
        self,
        object_id: str,
        location: str,
        start: datetime,
        end: datetime,
        stepsize: Quantity = 5 * u.minute,
    ) -> None:
        # check that the start and end time are timezone-aware
        if start.tzinfo is None or start.tzinfo.utcoffset(None) is None:
            raise ValueError("start must be a timezone-aware datetime.")
        if end.tzinfo is None or end.tzinfo.utcoffset(None) is None:
            raise ValueError("end must be a timezone-aware datetime.")

        # check that the start and end time are consistent
        if start >= end:
            raise ValueError("The start time must be earlier than the end time.")

        # avoid overly excessive queries
        if stepsize < 5 * u.minute:
            raise ValueError(
                "The time between ephemeris values must be at least 5 " "minutes"
            )

        self._object_id = object_id
        self._location = location
        self._start = start
        self._end = end
        self._stepsize = stepsize
        self._ephemerides: Optional[List[Ephemeris]] = None

    def ephemerides(
        self, start: Optional[datetime] = None, end: Optional[datetime] = None
    ) -> List[Ephemeris]:
        """
        Return the list of ephemerides covering an interval.

        The start and time of the interval must be timezone-aware, and the interval must
        be a subinterval of that passed to the constructor. The JPL Horizons service is
        only queried the first time you call this method; subsequent calls use cached
        results.

        If the start or end time is ``None``, the start or end time of the constructor
        is used.

        To ensure full coverage, the last ephemeris epoch returned may be slightly later
        than the end time.

        Parameters
        ----------
        start: `~datetime.datetime`, default: None
            The start time of the interval to cover. This must be a timezone-aware
            `~datetime.datetime` object.
        end: `~datetime.datetime`, default: None
            The end time of the interval to cover. This must be a timezone-aware
            `~datetime.datetime` object.

        Returns
        -------
        list of `~imephu.utils.Ephemeris`
            The list of ephemerides covering the interval from ``start`` to ``end``.
        """
        if start is None:
            start = self._start
        if end is None:
            end = self._end

        # check that the start and end time are timezone-aware
        if start.tzinfo is None or start.tzinfo.utcoffset(None) is None:
            raise ValueError("start must be a timezone-aware datetime.")
        if end.tzinfo is None or end.tzinfo.utcoffset(None) is None:
            raise ValueError("end must be a timezone-aware datetime.")

        # check that the start and end time are consistent
        if start >= end:
            raise ValueError("The start time must be earlier than the end time.")

        # check that the passed interval is a subinterval of that passed to the
        # constructor
        if start < self._start:
            raise ValueError(
                "The start time must not be earlier than that of the " "constructor."
            )
        if end > self._end:
            raise ValueError(
                "The end time must not be later than that of the" "constructor."
            )

        if self._ephemerides is None:
            self._query_horizons()

        return self._cover_time_interval(start, end)

    def _query_horizons(self) -> None:
        # Make sure the whole time interval is covered by the queried ephemerides
        start = self._start.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
        end_time_with_margin = self._end + timedelta(
            seconds=self._stepsize.to_value(u.second)
        )
        stop = end_time_with_margin.astimezone(timezone.utc).strftime(
            "%Y-%m-%d %H:%M:%S"
        )

        # Horizons requires an int for the step size. As round() might call NumPy's
        # round method and thus produce a float, we have to round "manually" using
        # the int function.
        step = f"{int(0.5 + self._stepsize.to_value(u.minute))}m"

        obj = Horizons(
            id=self._object_id,
            location=self._location,
            epochs={"start": start, "stop": stop, "step": step},
        )
        ephemerides = obj.ephemerides()

        # store the ephemerides in the format we need
        self._ephemerides = []
        for row in range(len(ephemerides)):
            epoch = parse(ephemerides["datetime_str"][row]).replace(tzinfo=timezone.utc)
            ra = float(ephemerides["RA"][row]) * u.deg
            dec = float(ephemerides["DEC"][row]) * u.deg
            if "V" in ephemerides.keys():
                magnitude = float(ephemerides["V"][row])
                magnitude_range = MagnitudeRange(
                    min_magnitude=magnitude, max_magnitude=magnitude, bandpass="V"
                )
            else:
                magnitude_range = None
            self._ephemerides.append(
                Ephemeris(
                    position=SkyCoord(ra=ra, dec=dec),
                    magnitude_range=magnitude_range,
                    epoch=epoch,
                )
            )

    def _cover_time_interval(self, start: datetime, end: datetime) -> List[Ephemeris]:
        """Find the smallest list of ephemerides covering a time interval."""
        # check that the start and end time are timezone-aware
        if start.tzinfo is None or start.tzinfo.utcoffset(None) is None:
            raise ValueError("start must be a timezone-aware datetime.")
        if end.tzinfo is None or end.tzinfo.utcoffset(None) is None:
            raise ValueError("end must be a timezone-aware datetime.")

        # check that JPL Horizons has been queried already
        if self._ephemerides is None:
            raise ValueError("self._ephemerides must not be None.")

        # check that the ephemerides cover the give interval
        if start < self._ephemerides[0].epoch or self._ephemerides[-1].epoch < end:
            raise ValueError("self._ephemerides does not cover the interval.")

        # find the smallest interval covering the time interval
        all_times = [e.epoch for e in self._ephemerides]
        start_index = bisect.bisect_right(all_times, start)
        end_index = bisect.bisect(all_times, end)
        if end_index > 0 and all_times[end_index - 1] == end:
            end_index -= 1

        return self._ephemerides[start_index - 1 : end_index + 1]
