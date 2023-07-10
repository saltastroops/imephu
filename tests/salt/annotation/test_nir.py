from pathlib import Path

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from imephu.annotation.general import TextAnnotation
from imephu.finder_chart import FinderChart
from imephu.salt.annotation.nir import bundles_annotation  # fiber_bundle_annotation


@pytest.mark.parametrize("position_angle", [0 * u.deg, 30 * u.deg, -90 * u.deg])
def test_nir_bundles_annotation(position_angle, fits_file, fits_center, check_finder):
    """Test the NIR annotation."""
    finder_chart = FinderChart(fits_file)
    nir_bundles_annotation_ = bundles_annotation(
        science_bundle_center=fits_center,
        bundle_separation=2.5 * u.arcmin,
        position_angle=position_angle,
        wcs=finder_chart.wcs,
    )
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"Position angle: {position_angle}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(nir_bundles_annotation_)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize("include_fibers", [False, True])
def test_nir_bundles_annotation_details(include_fibers, check_finder):
    """Test the NIR annotation with a zoomed in PANSTARRS finder chart."""
    fits_path = Path(__file__).parent.parent.parent / "data" / "ngc6000_1arcmin.fits"
    finder_chart = FinderChart(fits_path)
    center = SkyCoord(ra="15h49m49.5s", dec="-29d23m13s")
    nir_bundles_annotation_ = bundles_annotation(
        science_bundle_center=center,
        bundle_separation=20 * u.arcsec,
        position_angle=0 * u.deg,
        wcs=finder_chart.wcs,
        include_fibers=include_fibers,
    )
    finder_chart.add_annotation(nir_bundles_annotation_)
    check_finder(finder_chart)
