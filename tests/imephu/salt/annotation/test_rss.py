import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from imephu.annotation.general import TextAnnotation
from imephu.finder_chart import FinderChart
from imephu.salt.annotation import rss
from imephu.salt.utils import MosMask, MosMaskSlit


def test_rss_field_of_view_annotation(fits_file, fits_center, check_finder):
    """Test the RSS field of view annotation."""
    finder_chart = FinderChart(fits_file)
    rss_fov_annotation = rss.field_of_view_annotation(fits_center, wcs=finder_chart.wcs)
    finder_chart.add_annotation(rss_fov_annotation)
    check_finder(finder_chart)


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


@pytest.mark.parametrize("position_angle", [Angle(0 * u.deg), 70 * u.deg])
def test_mos_mask_annotation(
    position_angle, mos_mask_xml, fits_file, fits_center, check_finder
):
    """Test the MOS mask annotation."""
    finder_chart = FinderChart(fits_file)
    reference_stars = [
        SkyCoord(ra="00h39m30s", dec="-60d02m00s"),
        SkyCoord(ra="00h40m30s", dec="-60d02m00s"),
    ]
    slits = [
        MosMaskSlit(
            center=SkyCoord(ra="00h39m50s", dec="-60d03m00s"),
            width=Angle(2 * u.arcsec),
            height=Angle(10 * u.arcsec),
            tilt=0 * u.deg,
        ),
        MosMaskSlit(
            center=SkyCoord(ra="00h40m10s", dec="-60d03m00s"),
            width=Angle(4 * u.arcsec),
            height=Angle(20 * u.arcsec),
            tilt=-45 * u.deg,
        ),
    ]
    xml = mos_mask_xml(
        center=fits_center,
        position_angle=position_angle,
        reference_stars=reference_stars,
        slits=slits,
    )
    mos_mask = MosMask(xml)
    mask_annotation = rss.mos_mask_annotation(
        mos_mask, finder_chart.wcs, Angle(10 * u.arcsec)
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"position angle: {position_angle}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(mask_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)
