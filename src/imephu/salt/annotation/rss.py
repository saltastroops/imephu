from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import (
    CircleAnnotation,
    GroupAnnotation,
    RectangleAnnotation,
    TextAnnotation,
)
from imephu.geometry import translate
from imephu.salt.annotation.smi import smi_factory
from imephu.salt.utils import MosMask

_OUTLINE_COLOR = "red"


def field_of_view_annotation(fits_center: SkyCoord, wcs: WCS) -> GroupAnnotation:
    """
    Return an annotation with the RSS field of view.

    Parameters
    ----------
    fits_center: `~astropy.coordinates.SkyCoord`
        The center of the finder chart, in right ascension and declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.GroupAnnotation`
        An annotation with the RSS field of view.
    """
    fov_annotation = CircleAnnotation(
        fits_center, 4 * u.arcmin, wcs=wcs, edgecolor="green"
    )
    label_position = translate(fits_center, (-2.9, 2.9) * u.arcmin)
    name_annotation = TextAnnotation(
        label_position,
        "RSS",
        wcs=wcs,
        style="italic",
        weight="bold",
        size="large",
        horizontalalignment="left",
        color=(0, 0, 1),
    )
    return GroupAnnotation([fov_annotation, name_annotation])


def longslit_annotation(
    fits_center: SkyCoord,
    slit_width: Angle,
    slit_height: Angle,
    position_angle: Angle,
    wcs: WCS,
) -> RectangleAnnotation:
    """Return an annotation showing a longslit.

    Parameters
    ----------
    fits_center: `~astropy.coordinates.SkyCoord`
        The central position of the finder chart, in right ascension and declination
    slit_width: `~astropy.coordinates.Angle`
        The width of the slit, as an angle on the sky.
    slit_height: `~astropy.coordinates.Angle`
        The height of the slit, as an angle on the sky.
    position_angle: `~astropy.coordinates.Angle`
        The position angle, as an angle on the sky measured from north to east.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `imephu.annotation.general.RectangleAnnotation`
        The annotation showing the longslit.
    """
    return RectangleAnnotation(
        fits_center, slit_width, slit_height, wcs=wcs, edgecolor="red", alpha=0.5
    ).rotate(fits_center, position_angle)


def reference_star_annotation(reference_star: SkyCoord, wcs: WCS) -> CircleAnnotation:
    """Return a circle highlighting a reference star.

    Parameters
    ----------
    reference_star: `~astropy.coordinates.SkyCoord`
        The central position of the finder chart, in right ascension and declination
    wcs: `~astropy.wcs.WCS`
        WCS object.
    """
    return CircleAnnotation(
        reference_star,
        10 * u.arcsec,
        wcs,
        edgecolor="red",
        alpha=0.5,
    )


def mos_mask_annotation(
    mos_mask: MosMask, wcs: WCS, reference_star_box_width: Angle = 5 * u.arcsec
) -> GroupAnnotation:
    """Return the annotation for a MOS mask.

    The slits and boxes around the reference stars are included in the annotation.

    Parameters
    ----------
    mos_mask: `~imephu.salt.util.MosMask`
        The MOS mask.
    wcs: `~astropy.wcs.WCS`
        WCS object.
    reference_star_box_width: `~astropy.coordinates.Angle`, default: 5 arcseconds
        The width (and height) of the boxes around reference stars, as an angle on the
        sky.

    Returns
    -------
    `~imephu.annotation.general.GroupAnnotation`
        The annotation displaying the MOS mask.
    """
    mask_annotation = GroupAnnotation([])

    position_angle = mos_mask.position_angle
    for star in mos_mask.reference_stars:
        reference_star_annotation = RectangleAnnotation(
            center=star,
            width=reference_star_box_width,
            height=reference_star_box_width,
            wcs=wcs,
            edgecolor=(1, 1, 0),
            linewidth=2,
        ).rotate(star, position_angle)
        mask_annotation.add_item(reference_star_annotation)

    for slit in mos_mask.slits:
        slit_annotation = RectangleAnnotation(
            center=slit.center,
            width=slit.width,
            height=slit.height,
            wcs=wcs,
            edgecolor="red",
        ).rotate(slit.center, position_angle + slit.tilt)
        mask_annotation.add_item(slit_annotation)

    return mask_annotation


def smi_bundles_annotation(
    smi_barcode: str,
    fits_center: SkyCoord,
    science_bundle_center: SkyCoord,
    position_angle: Angle,
    wcs: WCS,
    include_fibers: bool = False,
) -> GroupAnnotation:
    """Return the annotation for NIR fiber bundles.

    The annotation contains the outline of the science fiber bundle (with its center
    marked by crosshairs) and the outlines of the sky fiber bundles.

    Parameters
    ----------
    smi_barcode: str
        The barcode of the slit mask IFU.
    fits_center: `~astropy.coordinates.SkyCoord`
        The finder chart center, as a position on the sky, in right ascension and
        declination.
    science_bundle_center: `~astropy.coordinates.SkyCoord`
        The center position of the science fiber bundle, as a position on the sky, in
        right ascension and declination.
    position_angle: `~astropy.coordinates.Angle`
        The position angle, as an angle on the sky measured from north to east.
    wcs:`astropy.wcs.WCS`
        WCS object.
    include_fibers: `bool`, optional, default False
        Whether to include the fiber outlines in the annotation.

    Returns
    -------
    `imephu.annotation.general.GroupAnnotation`
        The annotation for an NIR setup.
    """
    smi = smi_factory(smi_barcode)
    bundle = smi.bundle(
        center=fits_center,
        wcs=wcs,
        include_fibers=include_fibers,
    ).rotate(fits_center, 90 * u.deg + position_angle)

    target_offsets = fits_center.spherical_offsets_to(science_bundle_center)

    return bundle.translate(target_offsets)
