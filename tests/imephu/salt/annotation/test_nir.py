from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.finder_chart import FinderChart
from imephu.salt.annotation.nir import fiber_bundle_annotation, nir_annotation


def test_nir_fiber_bundle_annotation(fits_file, fits_center, check_finder):
    """Test the NIR fiber bundle annotation."""
    finder_chart = FinderChart(fits_file)
    fiber_centers = [
        SkyCoord(ra="00h40m00s", dec="-59d59m00s"),
        SkyCoord(ra="00h40m10s", dec="-60d00m00s"),
        SkyCoord(ra="00h40m00s", dec="-60d00m00s"),
        SkyCoord(ra="00h39m50s", dec="-60d00m00s"),
        SkyCoord(ra="00h40m00s", dec="-60d01m00s")
    ]
    bundle_annotation = fiber_bundle_annotation(center=fits_center, min_width=2.5 * u.arcmin, max_width=5 * u.arcmin, max_height=4 * u.arcmin, wcs=finder_chart.wcs, fiber_centers=fiber_centers, fiber_diameter=1 * u.arcmin)
    finder_chart.add_annotation(bundle_annotation)
    check_finder(finder_chart)


def test_nir_annotation(fits_file, fits_center, check_finder):
    """Test the NIR annotation."""
    finder_chart = FinderChart(fits_file)
    nir_annotation_ = nir_annotation(fits_center, finder_chart.wcs)
    finder_chart.add_annotation(nir_annotation_)
    check_finder(finder_chart)
