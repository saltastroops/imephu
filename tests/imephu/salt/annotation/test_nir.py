import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.annotation.general import TextAnnotation
from imephu.finder_chart import FinderChart
from imephu.salt.annotation.nir import fiber_bundle_annotation, nir_annotation


@pytest.mark.parametrize("with_crosshairs", [True, False])
def test_nir_fiber_bundle_annotation(with_crosshairs, fits_file, fits_center, check_finder):
    """Test the NIR fiber bundle annotation."""
    finder_chart = FinderChart(fits_file)
    bundle_annotation = fiber_bundle_annotation(center=fits_center, min_width=2.5 * u.arcmin, max_width=5 * u.arcmin, max_height=4 * u.arcmin, with_crosshairs=with_crosshairs, wcs=finder_chart.wcs)
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"With crosshairs: {with_crosshairs}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(bundle_annotation)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)


@pytest.mark.parametrize("position_angle", [0 * u.deg, 30 * u.deg, -90 * u.deg])
def test_nir_annotation(position_angle, fits_file, fits_center, check_finder):
    """Test the NIR annotation."""
    finder_chart = FinderChart(fits_file)
    nir_annotation_ = nir_annotation(science_bundle_center=fits_center, bundle_separation=2.5 * u.arcmin, position_angle=position_angle, wcs=finder_chart.wcs)
    legend = TextAnnotation(
        SkyCoord(ra="00h40m36s", dec="-59d55m30s"),
        f"Position angle: {position_angle}",
        wcs=finder_chart.wcs,
        color="blue",
        horizontalalignment="left",
    )
    finder_chart.add_annotation(nir_annotation_)
    finder_chart.add_annotation(legend)
    check_finder(finder_chart)
