"""pytest configuration."""

import io
import pathlib
from unittest import mock

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.general import TextAnnotation
from imephu.salt.finder_chart import FinderChart


@pytest.fixture(autouse=True)
def no_http_requests(monkeypatch):
    """Prevent any real HTTP requests.

    Taken (with wording slightly adapted) from
    https://blog.jerrycodes.com/no-http-requests/.
    """

    def urlopen_mock(self, method, url, *args, **kwargs):
        raise RuntimeError(
            f"The test was about to make a {method} request to "
            f"{self.scheme}://{self.host}{url}"
        )

    monkeypatch.setattr(
        "urllib3.connectionpool.HTTPConnectionPool.urlopen", urlopen_mock
    )


@pytest.fixture()
def check_finder(file_regression):
    """
    Return a function for checking finder charts.

    The finder chart is saved as a png, and the png is compared against a previously
    saved version. If no version exists already, the file is saved and the test fails.
    The saved file should be put under version control.

    If the saved png and the previously saved version differ, the test fails.

    In case you need to update the saved files, run ``pytest`` with the
    ``--force-regen`` flag.

    Parameters
    ----------
    file_regression: file regression fixture
        The file regression fixture from the pytest-regressions plugin.

    Returns
    -------
    function
        The function for checking a finder chart.
    """

    def _check_finder(finder_chart):
        np.random.seed(0)
        try:
            contents = io.BytesIO()
            finder_chart.save(contents, format="png")
            file_regression.check(contents.getvalue(), binary=True, extension=".png")
        finally:
            np.random.seed()

    return _check_finder


@pytest.fixture()
def fits_file():
    """
    Return the path of an example FITS file.

    The FITS file whose path is returned shows a 10 arcsecond by 10 arcsecond sky area
    centered on the right ascension 10 degrees and the declination -60 degrees.

    Returns
    -------
    `pathlib.Path`
        The path to the example FITS file.
    """
    return pathlib.Path(__file__).parent / "data" / "ra10_dec-60.fits"


@pytest.fixture()
def fits_center():
    """Return the sky coordinates for the center of the example FITS file."""
    return SkyCoord(ra=10 * u.deg, dec=-60 * u.deg)


@pytest.fixture()
def mos_mask_xml():
    """Return a function for generating XML describing a MOS mask."""

    def _mask_xml(center, position_angle, reference_stars, slits):
        xml = f"""\
<?xml version="1.0" ?>
<slitmask>
<header>
<parameter name="VERSION" value="1.1" />
<parameter name="PROPOSALCODE" value="INDEF" />
<parameter name="MASKNUM" value="0" />
<parameter name="PI" value="INDEF" />
<parameter name="CREATOR" value="Someone" />
<parameter name="ROTANGLE" value="{position_angle.to_value(u.deg)}" />
<parameter name="CENTERRA" value="{center.ra.to_value(u.deg)}" />
<parameter name="CENTERDEC" value="{center.dec.to_value(u.deg)}" />
<parameter name="EQUINOX" value="2000.0" />
<parameter name="NSMODE" value="0" />
<parameter name="COOSYS" value="RADEC" />
<parameter name="VALIDATED" value="FALSE" />
<parameter name="SPECLENGTH" value="12400" />
<parameter name="SPECOFFSET" value="0" />
<parameter name="SPECPOLSPLIT" value="0" />
<parameter name="SPECHEIGHT" value="0" />
</header>
"""
        id = 1
        for star in reference_stars:
            xml += f"""
    <refstar
        id="{id}"
        xce="{star.ra.to_value(u.deg)}"
        yce="{star.dec.to_value(u.deg)}"
        radius="0.5" mag="0.0"
    />"""
            id += 1

        for slit in slits:
            xml += f"""
    <slit
         id="{id}"
         xce="{slit.center.ra.to_value(u.deg)}"
         yce="{slit.center.dec.to_value(u.deg)}"
         width="{slit.width.to_value(u.arcsec)}"
         length="{slit.height.to_value(u.arcsec)}"
         tilt="{slit.tilt.to_value(u.deg)}"
         priority="1.0"
         mag="0.0"
    />"""

        xml += "</slitmask>"

        return xml

    return _mask_xml


@pytest.fixture()
def mock_salt_load_fits(fits_file):
    """Return a fixture for mocking the load_fits function for SALT finder charts.

    This fixture mocks the ``from_survey`` method of the
    `~imephu.finder_chart.FinderChart` class. The mock method always returns the FITS
     image of the `fits_file` fixture.

    .. warning::

       The mock function ignores any arguments - you always get the same FITS file. In
       particular thus implies that you always should use the `fits_center` fixture for
       the center of the FITS image.
    """
    with mock.patch.object(FinderChart, "from_survey", autospec=True) as mock_load_fits:
        mock_load_fits.return_value = FinderChart(open(fits_file, "rb"))
        yield mock_load_fits


@pytest.fixture()
def legend():
    """Return a fixture for adding a legend to a finder chart."""

    def _legend(text, wcs):
        return TextAnnotation(
            SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
            text,
            wcs=wcs,
            color="blue",
            horizontalalignment="left",
        )

    return _legend
