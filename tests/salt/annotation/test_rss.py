from itertools import product
from pathlib import Path

import pytest
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from imephu import geometry
from imephu.annotation.general import TextAnnotation
from imephu.finder_chart import FinderChart
from imephu.geometry import translate
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


def test_rss_reference_star_annotation(fits_file, fits_center, check_finder):
    finder_chart = FinderChart(fits_file)
    reference_star = geometry.translate(fits_center, (1, 3) * u.arcmin)
    reference_star_annotation = rss.reference_star_annotation(
        reference_star, finder_chart.wcs
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        "reference star 1 armin in right ascension and 3 arcmin in declination away "
        "from the center",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(reference_star_annotation)
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


@pytest.mark.parametrize(
    "smi_barcode, include_fibers", product(["PF0200N001", "PF0300N001"], [True, False])
)
def test_smi_annotation_barcode(
    smi_barcode, include_fibers, fits_file, fits_center, check_finder
):
    """Test the SMI annotation with different barcodes, with and without fibers."""
    finder_chart = FinderChart(fits_file)
    smi_annotation = rss.smi_bundles_annotation(
        smi_barcode,
        fits_center,
        fits_center,
        0 * u.deg,
        finder_chart.wcs,
        include_fibers,
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"SMI: {smi_barcode}, {'with' if include_fibers else 'without'} fibers",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(smi_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize("position_angle", [Angle(0 * u.deg), 70 * u.deg, -130 * u.deg])
def test_smi_annotation_position_angle(
    position_angle, fits_file, fits_center, check_finder
):
    """Test the SMI annotation with and without position angles."""
    finder_chart = FinderChart(fits_file)
    smi_annotation = rss.smi_bundles_annotation(
        "PF0200N001", fits_center, fits_center, position_angle, finder_chart.wcs, False
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"position angle: {position_angle}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(smi_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "smi_barcode, position_angle, include_fibers",
    [
        ("PF0200N001", 0 * u.deg, False),
        ("PF0200N001", 70 * u.deg, True),
        ("PF0300N001", 0 * u.deg, False),
        ("PF0300N001", 70 * u.deg, True),
    ],
)
def test_smi_annotation_details(
    smi_barcode, position_angle, include_fibers, check_finder
):
    """Test the SMI annotation with a zoomed in PANSTARRS finder chart."""
    fits_path = Path(__file__).parent.parent.parent / "data" / "ngc6000_1arcmin.fits"
    finder_chart = FinderChart(fits_path)
    center = SkyCoord(ra="15h49m49.5s", dec="-29d23m13s")
    smi_annotation = rss.smi_bundles_annotation(
        smi_barcode, center, center, position_angle, finder_chart.wcs, include_fibers
    )
    legend = TextAnnotation(
        SkyCoord(ra="15h49m51.6s", dec="-29d22m45s"),
        f"{smi_barcode}, {'with' if include_fibers else 'without'} fibers, position angle: {position_angle}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(smi_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize(
    "smi_barcode, dec_offset",
    [
        pytest.param("PF0200N001", -90 * u.arcsec, id="PF0200N001, sky bundle 1"),
        pytest.param("PF0200N001", 0 * u.arcsec, id="PF0200N001, target bundle"),
        pytest.param("PF0200N001", 90 * u.arcsec, id="PF0200N001, sky bundle 2"),
        pytest.param("PF0300N001", -90 * u.arcsec, id="PF0300N001, sky bundle 1"),
        pytest.param("PF0300N001", 0 * u.arcsec, id="PF0300N001, target bundle"),
        pytest.param("PF0300N001", 90 * u.arcsec, id="PF0300N001, sky bundle 2"),
    ],
)
def test_smi_annotation_bundle_part_details(smi_barcode, dec_offset, check_finder):
    """Test the SMI annotation with a zoomed in PANSTARRS finder chart."""
    fits_path = Path(__file__).parent.parent.parent / "data" / "ngc6000_1arcmin.fits"
    finder_chart = FinderChart(fits_path)
    center = SkyCoord(ra="15h49m49.5s", dec="-29d23m13s")
    ifu_center = translate(center, (0 * u.arcsec, dec_offset))
    smi_annotation = rss.smi_bundles_annotation(
        smi_barcode, center, ifu_center, 0 * u.deg, finder_chart.wcs, True
    )
    legend = TextAnnotation(
        SkyCoord(ra="15h49m51.6s", dec="-29d22m45s"),
        f"{smi_barcode}, position angle: 0.0 deg",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(smi_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)
