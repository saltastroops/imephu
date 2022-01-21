"""pytest configuration."""

import io
import pathlib

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy import units as u


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

    def _check_fits(finder_chart):
        np.random.seed(0)
        contents = io.BytesIO()
        finder_chart.save(contents, format="png")
        file_regression.check(contents.getvalue(), binary=True, extension=".png")

    return _check_fits


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
