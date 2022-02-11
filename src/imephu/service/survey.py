from io import BytesIO
from typing import BinaryIO, Protocol

import requests
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord


class SurveyError(BaseException):
    """An exception arising when querying a sky survey."""

    pass


_DSS_IDENTIFIERS = {
    "POSS2/UKSTU Red": "poss2ukstu_red",
    "POSS2/UKSTU Blue": "poss2ukstu_blue",
    "POSS2/UKSTU IR": "poss2ukstu_ir",
    "POSS1 Red": "poss1_red",
    "POSS1 Blue": "poss1_blue",
    "Quick-V": "quickv",
    "HST Phase2 (GSC2)": "phase2_gsc2",
    "HST Phase2 (GSC1)": "phase2_gsc1",
}


_SKYVIEW_IDENTIFIERS = {
    "2MASS-H": "2mass-h",
    "2MASS-J": "2mass-j",
    "2MASS-K": "2mass-k",
}


class SkySurvey(Protocol):
    """A sky survey from which fits files can be requested."""

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
        """
        raise NotImplementedError


class DigitizedSkySurvey(SkySurvey):
    """A class for loading FITS files from the Digitized Sky Survey (DSS).

    See the `DSS archive website <https://archive.stsci.edu/dss/index.html>` for more
    details.
    """

    def __init__(self) -> None:
        self._survey_identifiers = {
            key.lower(): value for key, value in _DSS_IDENTIFIERS.items()
        }

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
        response = requests.get(url, params=params)
        if response.status_code != 200:
            raise SurveyError("No FITS file could be loaded.")
        return BytesIO(response.content)

    def _survey_identifier(self, survey: str) -> str:
        if survey.lower() not in self._survey_identifiers:
            raise ValueError(f"Unknown survey: {survey}")

        return self._survey_identifiers[survey.lower()]


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
        self._survey_identifiers = {
            key.lower(): value for key, value in _SKYVIEW_IDENTIFIERS.items()
        }

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
        response = requests.get(url, params=params)
        if response.status_code != 200:
            raise SurveyError("No FITS file could be loaded.")
        return BytesIO(response.content)

    def _survey_identifier(self, survey: str) -> str:
        if survey.lower() not in self._survey_identifiers:
            raise ValueError(f"Unknown survey: {survey}")

        return self._survey_identifiers[survey.lower()]


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
    belongs to which `SkySurvey` implementation. The disadvantage is that the
    `SkySurvey` implementations may offer more functionality. (The `SkyView` constructor
    lets yo choose the FITS image size in pixels.)

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
    if survey.lower() in [s.lower() for s in _DSS_IDENTIFIERS.keys()]:
        survey_: SkySurvey = DigitizedSkySurvey()
    elif survey.lower() in [s.lower() for s in _SKYVIEW_IDENTIFIERS.keys()]:
        survey_ = SkyView(pixels=700)
    else:
        raise ValueError(f"Unknown survey: {survey}")

    return survey_.load_fits(survey, fits_center, size)
