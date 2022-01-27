import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from imephu.salt.utils import MosMask, MosMaskSlit


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


@pytest.mark.parametrize(
    "center,position_angle,reference_stars,slits",
    [
        (
            SkyCoord(ra="170d06m00s", dec="-55d30m00s"),
            35.67 * u.deg,
            [
                SkyCoord(ra="170d09m00s", dec="-55d33m00"),
                SkyCoord(ra="170d07m30s", dec="-55d36m00s"),
            ],
            [
                (
                    SkyCoord(ra="170d06m30s", dec="-55d30m30s"),
                    2 * u.arcsec,
                    10 * u.arcsec,
                    0 * u.deg,
                ),
                (
                    SkyCoord(ra="170d09m30s", dec="-55d33m30s"),
                    1 * u.arcsec,
                    5 * u.arcsec,
                    45.45 * u.deg,
                ),
            ],
        ),
        (
            SkyCoord(ra="170d06m00s", dec="-55d30m00s"),
            35.67 * u.deg,
            [],
            [
                (
                    SkyCoord(ra="170d06m30s", dec="-55d30m30s"),
                    2 * u.arcsec,
                    10 * u.arcsec,
                    0 * u.deg,
                ),
                (
                    SkyCoord(ra="170d09m30s", dec="-55d33m30s"),
                    1 * u.arcsec,
                    5 * u.arcsec,
                    45.45 * u.deg,
                ),
            ],
        ),
        (
            SkyCoord(ra="170d06m00s", dec="-55d30m00s"),
            35.67 * u.deg,
            [
                SkyCoord(ra="170d09m00s", dec="-55d33m00"),
                SkyCoord(ra="170d07m30s", dec="-55d36m00s"),
            ],
            [],
        ),
    ],
)
def test_init_mos_mask(center, position_angle, reference_stars, slits):
    """Test MOS mask initialisation."""
    mask_slits = []
    for slit_center, width, height, tilt in slits:
        mask_slits.append(
            MosMaskSlit(center=slit_center, width=width, height=height, tilt=tilt)
        )
    xml = _mask_xml(center, position_angle, reference_stars, mask_slits)
    mos_mask = MosMask(xml)

    assert mos_mask.center.ra.to_value(u.deg) == pytest.approx(
        center.ra.to_value(u.deg)
    )
    assert mos_mask.center.dec.to_value(u.deg) == pytest.approx(
        center.dec.to_value(u.deg)
    )
    assert mos_mask.position_angle.to_value(u.deg) == pytest.approx(
        position_angle.to_value(u.deg)
    )

    assert len(mos_mask.reference_stars) == len(reference_stars)
    for i, star in enumerate(mos_mask.reference_stars):
        assert star.ra.to_value(u.deg) == pytest.approx(
            reference_stars[i].ra.to_value(u.deg)
        )
        assert star.dec.to_value(u.deg) == pytest.approx(
            reference_stars[i].dec.to_value(u.deg)
        )

    assert len(mos_mask.slits) == len(slits)
    for i, slit in enumerate(mos_mask.slits):
        slit_center, width, height, tilt = slits[i]
        assert slit.center.ra.to_value(u.deg) == pytest.approx(
            slit_center.ra.to_value(u.deg)
        )
        assert slit.center.dec.to_value(u.deg) == pytest.approx(
            slit_center.dec.to_value(u.deg)
        )
        assert slit.width.to_value(u.arcsec) == pytest.approx(width.to_value(u.arcsec))
        assert slit.height.to_value(u.arcsec) == pytest.approx(
            height.to_value(u.arcsec)
        )
        assert slit.tilt.to_value(u.deg) == pytest.approx(tilt.to_value(u.deg))
