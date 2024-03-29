from io import BytesIO

import pytest
import responses
from astropy import units as u
from astropy.coordinates import SkyCoord
from imephu.service.survey import (
    DigitizedSkySurvey,
    SkyView,
    SurveyError,
    is_covering_position,
    load_fits,
)
from responses import matchers

DSS_URL = "https://archive.stsci.edu/cgi-bin/dss_search"

SKYVIEW_URL = "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"


@responses.activate
def test_dss_calls_the_correct_url():
    """Test that the DigitizedSkySurvey class calls the correct URL."""
    expected_params = {
        "v": "poss2ukstu_red",
        "r": "120.0",
        "d": "-30.0",
        "e": "J2000",
        "h": "10.0",
        "w": "10.0",
        "f": "fits",
        "c": "none",
    }
    responses.add(
        responses.GET,
        DSS_URL,
        body=b"I'm a fake FITS file",
        match=[matchers.query_param_matcher(expected_params)],
    )
    dss = DigitizedSkySurvey()
    dss.load_fits(
        "POSS2/UKSTU Red", SkyCoord(ra=120 * u.deg, dec=-30 * u.deg), 10 * u.arcmin
    )


@responses.activate
def test_dss_returns_the_content_sent():
    """Test that the DigitizedSkySurvey class returns the received content."""
    content = b"I'm a fake FITS file"
    responses.add(responses.GET, DSS_URL, body=content)
    dss = DigitizedSkySurvey()
    fits = dss.load_fits(
        "POSS2/UKSTU Red", SkyCoord(ra=0 * u.deg, dec=0 * u.deg), 10 * u.arcmin
    )
    assert fits.read() == BytesIO(b"I'm a fake FITS file").read()


def test_dss_handles_wrong_survey():
    """Test that the DigitizedSkySurvey class handles wrong survey names correctly."""
    dss = DigitizedSkySurvey()
    with pytest.raises(SurveyError):
        dss.load_fits(
            "Invalid Survey", SkyCoord(ra=0 * u.deg, dec=0 * u.deg), 10 * u.arcmin
        )


@responses.activate
@pytest.mark.parametrize("status_code", [201, 304, 400, 401, 403, 404, 500])
def test_dss_handles_http_errors(status_code):
    """Test that the DigitizedSkySurvey class handles HTTP errors correctly."""
    responses.add(responses.GET, DSS_URL, status=status_code, body="Some error")
    with pytest.raises(SurveyError):
        dss = DigitizedSkySurvey()
        dss.load_fits(
            "POSS2/UKSTU Red", SkyCoord(ra=0 * u.deg, dec=0 * u.deg), 10 * u.arcmin
        )


@responses.activate
def test_skyview_calls_the_correct_url():
    """Test that the SkyView class calls the correct URL."""
    expected_params = {
        "Position": "120.0, -30.0",
        "Survey": "2mass-k",
        "Coordinates": "ICRS",
        "Return": "FITS",
        "Pixels": "700",
        "Size": str((10 * u.arcmin).to_value(u.deg)),
    }
    responses.add(
        responses.GET,
        SKYVIEW_URL,
        body=b"I'm a fake FITS file",
        match=[matchers.query_param_matcher(expected_params)],
    )
    dss = SkyView(pixels=700)
    dss.load_fits("2MASS-K", SkyCoord(ra=120 * u.deg, dec=-30 * u.deg), 10 * u.arcmin)


@responses.activate
def test_skyview_returns_the_content_sent():
    """Test that the SkyView class returns the received content."""
    content = b"I'm a fake FITS file"
    responses.add(responses.GET, SKYVIEW_URL, body=content)
    dss = SkyView()
    fits = dss.load_fits(
        "2MASS-K", SkyCoord(ra=0 * u.deg, dec=0 * u.deg), 10 * u.arcmin
    )
    assert fits.read() == BytesIO(b"I'm a fake FITS file").read()


def test_skyview_handles_wrong_survey():
    """Test that the SkyView class handles wrong survey names correctly."""
    dss = SkyView()
    with pytest.raises(BaseException):
        dss.load_fits(
            "Invalid Survey", SkyCoord(ra=0 * u.deg, dec=0 * u.deg), 10 * u.arcmin
        )


@responses.activate
@pytest.mark.parametrize("status_code", [201, 304, 400, 401, 403, 404, 500])
def test_skyview_survey_handles_http_errors(status_code):
    """Test that the SkyView class handles HTTP errors correctly."""
    responses.add(responses.GET, SKYVIEW_URL, status=status_code, body="Some error")
    with pytest.raises(SurveyError):
        dss = SkyView()
        dss.load_fits("2MASS-K", SkyCoord(ra=0 * u.deg, dec=0 * u.deg), 10 * u.arcmin)


def test_load_fits_handles_invalid_survey(fits_center):
    """Test that the load_fits function handles invalid surveys correctly."""
    with pytest.raises(ValueError):
        load_fits(survey="Invalid Survey", fits_center=fits_center, size=10 * u.arcmin)


@pytest.mark.parametrize(
    "survey, ra, dec, expected",
    [
        ("POSS2/UKSTU Red", 1, -60, True),
        ("POSS2/UKSTU Blue", 1, 5, True),
        ("POSS1 Red", 1, -30.001, False),
        ("POSS1 Red", 1, -29.999, True),
        ("POSS1 Blue", 1, -30.01, False),
        ("POSS1 Blue", 1, -29.99, True),
        ("Quick-V", 1, 5.99, False),
        ("Quick-V", 1, 6.01, True),
        ("HST Phase2 (GSC2)", 1, -56, True),
        ("HST Phase2 (GSC1)", 1, -20.01, True),
        ("HST Phase2 (GSC1)", 1, -19.99, False),
        ("HST Phase2 (GSC1)", 1, 5.99, False),
        ("HST Phase2 (GSC1)", 1, 6.01, True),
        ("2MASS-H", 1, 5, True),
        ("2MASS-J", 1, 5, True),
        ("2MASS-K", 1, 5, True),
    ],
)
def test_is_covering_position_returns_correct_value(
    survey: str, ra: float, dec: float, expected: bool
) -> None:
    assert (
        is_covering_position(survey, SkyCoord(ra=ra * u.deg, dec=dec * u.deg))
        == expected
    )
