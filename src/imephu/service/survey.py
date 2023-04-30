import urllib.parse
from io import BytesIO
from typing import BinaryIO, Protocol, NamedTuple, Callable, Dict

import requests
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord


class SurveyError(BaseException):
    """An exception arising when querying a sky survey."""

    pass


class _SurveyDetails(NamedTuple):
    identifier: str
    is_covering_position: Callable[[SkyCoord], bool]


_always_true = lambda x: True


_DSS_DETAILS: Dict[str, _SurveyDetails] = {
    "POSS2/UKSTU Red": _SurveyDetails("poss2ukstu_red", _always_true),
    "POSS2/UKSTU Blue":  _SurveyDetails("poss2ukstu_blue", _always_true),
    "POSS2/UKSTU IR": _SurveyDetails("poss2ukstu_ir", _always_true),
    "POSS1 Red": _SurveyDetails("poss1_red", lambda p: p.dec.to_value(u.deg) > -30),
    "POSS1 Blue": _SurveyDetails("poss1_blue", lambda p: p.dec.to_value(u.deg) > -30),
    "Quick-V": _SurveyDetails("quickv", lambda p: p.dec.to_value(u.deg) > 6),
    "HST Phase2 (GSC2)": _SurveyDetails("phase2_gsc2", _always_true),
    "HST Phase2 (GSC1)": _SurveyDetails("phase2_gsc1", lambda p: p.dec.to_value(u.deg) < -20 or p.dec.to_value(u.deg) > 6)
}


_SKYVIEW_DETAILS = {
    "2MASS-H": _SurveyDetails("2mass-h", _always_true),
    "2MASS-J": _SurveyDetails("2mass-j", _always_true),
    "2MASS-K": _SurveyDetails("2mass-k", _always_true)
}


class SkySurvey(Protocol):
    """A sky survey from which fits files can be requested."""

    def url(self, survey: str, fits_center: SkyCoord, size: Angle) -> str:
        """Return the URL for loading the FITS file for a given position.

        See `load_fits` for the available surveys.

        Parameters
        ----------
        survey: `str`
            The name of the survey to query for the FITS file.
        fits_center: `~astropy.coordinates.SkyCoord`
            The center position of the loaded FITS file.
        size: `~astropy.coordinates.Angle`
            The width and height of the FITS image, as an angle on the sky.

        Returns
        -------
        str
            The URL.
        """
        raise NotImplementedError

    def load_fits(self, survey: str, fits_center: SkyCoord, size: Angle) -> BinaryIO:
        """Load a FITS file for a given position from the sky survey.

        The size of the returned FITS file is determined by the class implementing this
        protocol (or the web service used).

        Parameters
        ----------
        survey: `str`
            The name of the survey to query for the FITS file.
        fits_center: `~astropy.coordinates.SkyCoord`
            The center position of the loaded FITS file.
        size: `~astropy.coordinates.Angle`
            The width and height of the FITS image, as an angle on the sky.

        Returns
        -------
        binary stream
            The FITS file.
        """
        raise NotImplementedError

    def is_covering_position(self, survey: str, position: SkyCoord) -> bool:
        """Return whether a position is covered by the sky survey.

        Parameters
        ----------
        survey: `str`
            The name of the survey.
        position: `~astropy.coordinates.SkyCoord`
            The position.

        Returns
        -------
        bool
            True if the position is covered  by the survey, False otherwise.
        """
        raise NotImplementedError


class DigitizedSkySurvey(SkySurvey):
    """A class for loading FITS files from the Digitized Sky Survey (DSS).

    See the `DSS archive website <https://archive.stsci.edu/dss/index.html>` for more
    details.
    """

    def __init__(self) -> None:
        self._survey_details: Dict[str, _SurveyDetails] = {
            key.lower(): value for key, value in _DSS_DETAILS.items()
        }

    def url(self, survey: str, fits_center: SkyCoord, size: Angle) -> str:
        """Return the URL for loading the FITS file for a given position.

        See `load_fits` for the available surveys.

        Parameters
        ----------
        survey: `str`
            The name of the survey to query for the FITS file.
        fits_center: `~astropy.coordinates.SkyCoord`
            The center position of the loaded FITS file.
        size: `~astropy.coordinates.Angle`
            The width and height of the FITS image, as an angle on the sky.

        Returns
        -------
        str
            The URL.
        """
        url = "https://archive.stsci.edu/cgi-bin/dss_search"
        params = {
            "v": self._survey_identifier(survey),
            "r": fits_center.ra.to_value(u.deg),
            "d": fits_center.dec.to_value(u.deg),
            "e": "J2000",
            "h": size.to_value(u.arcmin),
            "w": size.to_value(u.arcmin),
            "f": "fits",
            "c": "none",
        }
        return url + "?" + urllib.parse.urlencode(params)

    def load_fits(self, survey: str, fits_center: SkyCoord, size: Angle) -> BinaryIO:
        """Load a FITS file for a given position from the sky survey.

        The following surveys can be queried for FITS images:

        - POSS2/UKSTU Red
        - POSS2/UKSTU Blue
        - POSS2/UKSTU IR
        - POSS1 Red
        - POSS1 Blue
        - Quick-V
        - HST Phase2 (GSC2)
        - HST Phase2 (GSC1)

        See the `DSS survey page <http://gsss.stsci.edu/SkySurveys/Surveys.htm>` for
        details about these surveys.

        Parameters
        ----------
        survey: `str`
            The name of the survey to query for the FITS file.
        fits_center: `~astropy.coordinates.SkyCoord`
            The center position of the loaded FITS file.
        size: `~astropy.coordinates.Angle`
            The width and height of the FITS image, as an angle on the sky.

        Returns
        -------
        binary stream
            The FITS file.
        """
        response = requests.get(self.url(survey, fits_center, size))
        if response.status_code != 200:
            raise SurveyError("No FITS file could be loaded.")
        return BytesIO(response.content)

    def is_covering_position(self, survey: str, position: SkyCoord) -> bool:
        """Return whether a position is covered by the sky survey.

        Parameters
        ----------
        survey: `str`
            The name of the survey.
        position: `~astropy.coordinates.SkyCoord`
            The position.

        Returns
        -------
        bool
            True if the position is covered  by the survey, False otherwise.
        """
        return self._survey_details[survey.lower()].is_covering_position(position)


    def _survey_identifier(self, survey: str) -> str:
        if survey.lower() not in self._survey_details:
            raise ValueError(f"Unknown survey: {survey}")

        return self._survey_details[survey.lower()].identifier


class SkyView(SkySurvey):
    """A class for loading FITS files from SkyView.

    See the `SkyView site <https://skyview.gsfc.nasa.gov/current/cgi/titlepage.pl>` for
    more details about SkyView.

    Parameters
    ----------
    pixels: `int`, default: 300
        The image size in pixels.
    """

    def __init__(self, pixels: int = 300) -> None:
        self._pixels = pixels
        self._survey_details: Dict[str, _SurveyDetails] = {
            key.lower(): value for key, value in _SKYVIEW_DETAILS.items()
        }

    def url(self, survey: str, fits_center: SkyCoord, size: Angle) -> str:
        """Return the URL for loading the FITS file for a given position.

        See `load_fits` for the available surveys.

        Parameters
        ----------
        survey: `str`
            The name of the survey to query for the FITS file.
        fits_center: `~astropy.coordinates.SkyCoord`
            The center position of the loaded FITS file.
        size: `~astropy.coordinates.Angle`
            The width and height of the FITS image, as an angle on the sky.

        Returns
        -------
        str
            The URL.
        """
        url = "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
        params = {
            "Position": f"{fits_center.ra.to_value(u.deg)}, "
            f"{fits_center.dec.to_value(u.deg)}",
            "Survey": self._survey_identifier(survey),
            "Coordinates": "ICRS",
            "Return": "FITS",
            "Pixels": self._pixels,
            "Size": size.to_value(u.deg),
        }
        return url + "?" + urllib.parse.urlencode(params)

    def load_fits(self, survey: str, fits_center: SkyCoord, size: Angle) -> BinaryIO:
        """Load a FITS file for a given position from the sky survey.

        The following subset of surveys supported by SkyView can be queried:

        - 2MASS-H
        - 2MASS-J
        - 2MASS-K

        See the `SkyView survey page
        <https://skyview.gsfc.nasa.gov/current/cgi/survey.pl>` for details about these
        surveys.

        Parameters
        ----------
        survey: `str`
            The name of the survey to query for the FITS file.
        fits_center: `~astropy.coordinates.SkyCoord`
            The center position of the loaded FITS file.
        size: `~astropy.coordinates.Angle`
            The width and height of the FITS image, as an angle on the sky.

        Returns
        -------
        binary stream
            The FITS file.
        """
        response = requests.get(self.url(survey, fits_center, size))
        if response.status_code != 200:
            raise SurveyError("No FITS file could be loaded.")
        return BytesIO(response.content)

    def is_covering_position(self, survey: str, position: SkyCoord) -> bool:
        """Return whether a position is covered by the sky survey.

        Parameters
        ----------
        survey: `str`
            The name of the survey.
        position: `~astropy.coordinates.SkyCoord`
            The position.

        Returns
        -------
        bool
            True if the position is covered  by the survey, False otherwise.
        """
        return self._survey_details[survey.lower()].is_covering_position(position)

    def _survey_identifier(self, survey: str) -> str:
        if survey.lower() not in self._survey_details:
            raise ValueError(f"Unknown survey: {survey}")

        return self._survey_details[survey.lower()].identifier


def url(survey: str, fits_center: SkyCoord, size: Angle) -> str:
    """Return the URL for requesting a FITS image from a sky survey.

    In principle, this function is equivalent to the `url` methods of the
    `SkySurvey` implementations. For example,

    .. code:: python

       from astropy import units as u
       from astropy.coordinates import SkyCoord

       url("POSS2/UKSTU Red",
           SkyCoord(ra=120 * u.deg, dec=-30 * u.deg),
           10 * u.arcmin)

    is the same as

    .. code:: python

       survey = DigitizedSkySurvey()
       survey.url("POSS2/UKSTU Red",
                  SkyCoord(ra=120 * u.deg, dec=-30 * u.deg),
                  10 * u.arcmin)

    The advantage of this function is that you don't have to remember which survey
    belongs to which `SkySurvey` implementation.

    See the `SkySurvey` implementations for the supported surveys.

    Parameters
    ----------
    survey: `str`
        The name of the survey to query for the FITS file.
    fits_center: `~astropy.coordinates.SkyCoord`
        The center position of the loaded FITS file.
    size: `~astropy.coordinates.Angle`
        The width and height of the FITS image, as an angle on the sky.

    Returns
    -------
    str
        The URL.
    """
    if survey.lower() in [s.lower() for s in _DSS_DETAILS.keys()]:
        survey_: SkySurvey = DigitizedSkySurvey()
    elif survey.lower() in [s.lower() for s in _SKYVIEW_DETAILS.keys()]:
        survey_ = SkyView(pixels=700)
    else:
        raise ValueError(f"Unknown survey: {survey}")

    return survey_.url(survey, fits_center, size)


def load_fits(survey: str, fits_center: SkyCoord, size: Angle) -> BinaryIO:
    """Request a FITS image from a sky survey.

    In principle, this function is equivalent to the `load_fits` methods of the
    `SkySurvey` implementations. For example,

    .. code:: python

       from astropy import units as u
       from astropy.coordinates import SkyCoord

       load_fits("POSS2/UKSTU Red",
                 SkyCoord(ra=120 * u.deg, dec=-30 * u.deg),
                 10 * u.arcmin)

    is the same as

    .. code:: python

       survey = DigitizedSkySurvey()
       survey.load_fits("POSS2/UKSTU Red",
                        SkyCoord(ra=120 * u.deg, dec=-30 * u.deg),
                        10 * u.arcmin)

    The advantage of this function is that you don't have to remember which survey
    belongs to which `SkySurvey` implementation.

    See the `SkySurvey` implementations for the supported surveys.

    Parameters
    ----------
    survey: `str`
        The name of the survey to query for the FITS file.
    fits_center: `~astropy.coordinates.SkyCoord`
        The center position of the loaded FITS file.
    size: `~astropy.coordinates.Angle`
        The width and height of the FITS image, as an angle on the sky.

    Returns
    -------
    binary stream
        The FITS file.
    """
    response = requests.get(url(survey, fits_center, size))
    if response.status_code != 200:
        raise SurveyError("No FITS file could be loaded.")
    return BytesIO(response.content)


def is_covering_position(survey: str, position: SkyCoord) -> bool:
    """Return whether a position is covered by a sky survey.

    In principle, this function is equivalent to the `is_covering_position` methods of
    the `SkySurvey` implementations. For example,

    .. code:: python

       from astropy import units as u
       from astropy.coordinates import SkyCoord

       is_covering_position("POSS2/UKSTU Red",
                            SkyCoord(ra=120 * u.deg, dec=-30 * u.deg))

    is the same as

    .. code:: python

       survey = DigitizedSkySurvey()
       survey.is_covering_position("POSS2/UKSTU Red",
                                   SkyCoord(ra=120 * u.deg, dec=-30 * u.deg))

    The advantage of this function is that you don't have to remember which survey
    belongs to which `SkySurvey` implementation.

    See the `SkySurvey` implementations for the supported surveys.

    Parameters
    ----------
    survey: `str`
        The name of the survey.
    position: `~astropy.coordinates.SkyCoord`
        The position.

    Returns
    -------
    binary stream
        The FITS file.
    """
    if survey.lower() in [s.lower() for s in _DSS_DETAILS.keys()]:
        survey_: SkySurvey = DigitizedSkySurvey()
    elif survey.lower() in [s.lower() for s in _SKYVIEW_DETAILS.keys()]:
        survey_ = SkyView(pixels=700)
    else:
        raise ValueError(f"Unknown survey: {survey}")
    return survey_.is_covering_position(survey, position)

