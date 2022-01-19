import math
import pytest
from astropy import units as u

from imephu.finder_chart import FinderChart
from imephu.geometry import rotate, translate


@pytest.mark.parametrize("v,pivot,angle,expected_rotated_v",
                         [
                             ((10.5 * u.deg, -44.5 * u.deg), (10.7 * u.deg, -44.8 * u.deg), 0 * u.deg, (10.5 * u.deg, -44.5 * u.deg)),
                             ((10.1 * u.deg, -45 * u.deg), (10 * u.deg, -45 * u.deg), 90 * u.deg, (10. * u.deg, -44.9 * u.deg))
                         ])
def test_geometry_rotate(v, pivot, angle, expected_rotated_v, fits_file):
    finder_chart = FinderChart(fits_file)
    wcs = finder_chart._wcs
    v_rotated = rotate(v, pivot, angle, wcs)
    expected_ra_deg = expected_rotated_v[0].to_value(u.deg)
    expected_dec_deg = expected_rotated_v[0].to_value(u.deg)
    ra_deg = v_rotated[0].to_value(u.deg)
    dec_deg = v_rotated[0].to_value(u.deg)

    # Due to non-linearity, the accuracy is limited
    assert (ra_deg, dec_deg) == pytest.approx((expected_ra_deg, expected_dec_deg), rel=1e-4)


@pytest.mark.parametrize("v,displacement,expected_translated_v",
                         [
                             ((9.93 * u.deg, -44.99 * u.deg), (0 * u.deg, 0 * u.deg), (9.93 * u.deg, -44.99 * u.deg)),
                             ((10.1 * u.deg, -45.02 * u.deg), (0.005 * u.deg, -0.02 * u.deg), (10.105 * u.deg, (-45.02 - 0.02 / math.sqrt(0.5)) * u.deg))
                         ])
def test_geometry_translate(v, displacement, expected_translated_v, fits_file):
    finder_chart = FinderChart(fits_file)
    wcs = finder_chart._wcs
    v_translated = translate(v, displacement, wcs)
    expected_ra_deg = expected_translated_v[0].to_value(u.deg)
    expected_dec_deg = expected_translated_v[1].to_value(u.deg)
    ra_deg = v_translated[0].to_value(u.deg)
    dec_deg = v_translated[1].to_value(u.deg)
    assert (ra_deg, dec_deg) == pytest.approx((expected_ra_deg, expected_dec_deg))
