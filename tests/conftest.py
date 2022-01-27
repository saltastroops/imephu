"""pytest configuration."""

import io
import pathlib

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord


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
