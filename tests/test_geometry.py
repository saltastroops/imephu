import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.finder_chart import FinderChart
from imephu.geometry import rotate, translate


@pytest.mark.parametrize(
    "v,pivot,angle,expected_rotated_v",
    [
        (
            SkyCoord(ra=10.5 * u.deg, dec=-60.5 * u.deg),
            SkyCoord(ra=10.7 * u.deg, dec=-59.8 * u.deg),
            0 * u.deg,
            SkyCoord(ra=10.5 * u.deg, dec=-60.5 * u.deg),
        ),
        (
            SkyCoord(ra=10.2 * u.deg, dec=-60 * u.deg),
            SkyCoord(ra=10 * u.deg, dec=-60 * u.deg),
            90 * u.deg,
            SkyCoord(ra=10.0 * u.deg, dec=-60.1 * u.deg),
        ),
    ],
)
def test_geometry_rotate(v, pivot, angle, expected_rotated_v, fits_file):
    """Test rotations on the sky."""
    finder_chart = FinderChart(fits_file)
    wcs = finder_chart._wcs
    v_rotated = rotate(v, pivot, angle, wcs)
    expected_ra_deg = expected_rotated_v.ra.to_value(u.deg)
    expected_dec_deg = expected_rotated_v.dec.to_value(u.deg)
    ra_deg = v_rotated.ra.to_value(u.deg)
    dec_deg = v_rotated.dec.to_value(u.deg)

    # Due to non-linearity, the accuracy is limited
    assert (ra_deg, dec_deg) == pytest.approx(
        (expected_ra_deg, expected_dec_deg), rel=1e-4
    )


@pytest.mark.parametrize(
    "v,displacement,expected_translated_v",
    [
        (
            SkyCoord(ra=9.93 * u.deg, dec=-59.99 * u.deg),
            (0, 0) * u.deg,
            SkyCoord(ra=9.93 * u.deg, dec=-59.99 * u.deg),
        ),
        (
            SkyCoord(ra=10.1 * u.deg, dec=-60.02 * u.deg),
            (0.03, -0.02) * u.deg,
            SkyCoord(ra=(10.1 + 0.06) * u.deg, dec=(-60.02 - 0.02) * u.deg),
        ),
    ],
)
def test_geometry_translate(v, displacement, expected_translated_v, fits_file):
    """Test translations on the sky."""
    v_translated = translate(v, displacement)
    expected_ra_deg = expected_translated_v.ra.to_value(u.deg)
    expected_dec_deg = expected_translated_v.dec.to_value(u.deg)
    ra_deg = v_translated.ra.to_value(u.deg)
    dec_deg = v_translated.dec.to_value(u.deg)
    assert (ra_deg, dec_deg) == pytest.approx(
        (expected_ra_deg, expected_dec_deg), rel=1e-5
    )
