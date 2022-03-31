"""pytest configuration."""
import io
import pathlib
import time
from unittest import mock

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from typer.testing import CliRunner

import imephu
from imephu.annotation.general import TextAnnotation
from imephu.cli import app
from imephu.salt.finder_chart import FinderChart

runner = CliRunner()


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
def check_cli(fits_file, tmp_path_factory, file_regression):
    """
    Return a function for checking the command line interface.

    Parameters
    ----------
    tmp_path_factory: fixture for creating a temporary directory
        Temporary directory.
    file_regression: fixture for regression checking
        Fixture for file regression checking.

    Returns
    -------
    function
        Function for checking the command line interface.
    """

    def _check_cli(
        instrument_yaml,
        fits_source_yaml="fits-source:\n  image-survey: POSS2/UKSTU Red",
    ):
        configuration = f"""\
{fits_source_yaml}
telescope: SALT
pi-family-name: Doe
proposal-code: 2022-1-SCI-042
position-angle: 30d
target:
  name: Magrathea
  ra: 0h 40m 00s
  dec: -60d
  magnitude-range:
    bandpass: V
    minimum: 17
    maximum: 17.3
{instrument_yaml}
"""
        np.random.seed(0)
        try:
            tmp = tmp_path_factory.mktemp(f"finder-chart-{time.time_ns()}")
            config = tmp / "config.yaml"
            config.write_text(configuration)
            output = tmp / "finder_chart.png"
            with mock.patch.object(
                imephu.cli, "load_fits", autospec=True
            ) as mock_load_fits:
                fits = fits_file.read_bytes()
                mock_load_fits.return_value = io.BytesIO(fits)
                runner.invoke(app, ["--config", config, "--out", output])
                finder_chart = output.read_bytes()
                file_regression.check(finder_chart, binary=True, extension=".png")
        finally:
            np.random.seed()

    return _check_cli


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
def fits_file2():
    """
    Return the path of an example FITS file.

    The FITS file whose path is returned shows a 10 arcsecond by 10 arcsecond sky area
    centered on the right ascension 9.75 degrees and the declination -60 degrees.

    Returns
    -------
    `pathlib.Path`
        The path to the example FITS file.
    """
    return pathlib.Path(__file__).parent / "data" / "ra9.75_dec-60.fits"


@pytest.fixture()
def fits_center():
    """Return the sky coordinates for the center of the example FITS file."""
    return SkyCoord(ra=10 * u.deg, dec=-60 * u.deg)


@pytest.fixture()
def fits_center2():
    """Return the sky coordinates for the center of the example FITS file."""
    return SkyCoord(ra=9.75 * u.deg, dec=-60 * u.deg)


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
def mock_from_survey(fits_file, fits_file2):
    """Return a fixture for mocking getting a finder chart from an image survey.

    This fixture mocks the ``from_survey`` method of the
    `~imephu.finder_chart.FinderChart` class. The mock method always returns a finder
    chart with the FITS image of the `fits_file` fixture when called the first time and
    a finder chart with the FITS image of the `fits_file2` fixture when called the
    second time.

    .. warning::

       The mock function ignores any arguments - you always get the same finder chart.
       In particular this implies that you always should use the `fits_center` fixture
       for the center of the FITS image when calling the function for the first time,
       and the fits_center2 fixture when calling it for the second time.
    """
    with mock.patch.object(FinderChart, "from_survey", autospec=True) as mock_load_fits:
        mock_load_fits.side_effect = [
            FinderChart(open(fits_file, "rb")),
            FinderChart(open(fits_file2, "rb")),
        ]
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
