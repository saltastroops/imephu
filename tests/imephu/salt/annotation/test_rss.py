import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.general import TextAnnotation
from imephu.finder_chart import FinderChart
from imephu.salt.annotation import rss


@pytest.mark.parametrize(
    "slit_width,slit_height,position_angle",
    [
        (2 * u.arcsec, 8 * u.arcmin, 10 * u.deg),
        (4 * u.arcsec, 4 * u.arcmin, 70 * u.deg),
    ],
)
def test_rss_longslit_annotation(
    slit_width, slit_height, position_angle, fits_file, fits_center, check_finder
):
    """Test RSS longslit annotations."""
    finder_chart = FinderChart(fits_file)
    longslit_annotation = rss.longslit_annotation(
        fits_center=fits_center,
        slit_width=slit_width,
        slit_height=slit_height,
        position_angle=position_angle,
        wcs=finder_chart.wcs,
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"slit width: {slit_width}, slit height: {slit_height}, position angle: "
        f"{position_angle}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(longslit_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)
