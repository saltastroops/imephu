from typing import Iterable

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import (
    CircleAnnotation,
    GroupAnnotation,
    RectangleAnnotation,
    TextAnnotation,
    LinePathAnnotation,
)
from imephu.geometry import translate
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
    smi = _SMI()
    bundle = smi.bundle(
        center=fits_center,
        wcs=wcs,
        include_fibers=include_fibers,
    ).rotate(fits_center, 90 * u.deg + position_angle)

    target_offsets = fits_center.spherical_offsets_to(science_bundle_center)

    return bundle.translate(target_offsets)


class _SMI:
    FIBER_RADIUS = 0.44 * u.arcsec

    FIBER_CENTERS = [
        (-4.345 * u.arcsec, -0.04 * u.arcsec),
        (-3.827 * u.arcsec, 0.841 * u.arcsec),
        (-3.248 * u.arcsec, -3.854 * u.arcsec),
        (-96.903 * u.arcsec, 0.164 * u.arcsec),
        (-1.602 * u.arcsec, -2.836 * u.arcsec),
        (-0.527 * u.arcsec, 2.739 * u.arcsec),
        (-1.619 * u.arcsec, -4.792 * u.arcsec),
        (-3.77 * u.arcsec, -0.969 * u.arcsec),
        (-1.044 * u.arcsec, 3.726 * u.arcsec),
        (-3.279 * u.arcsec, 0.004 * u.arcsec),
        (-1.071 * u.arcsec, -3.854 * u.arcsec),
        (-2.111 * u.arcsec, 3.726 * u.arcsec),
        (-2.142 * u.arcsec, -1.907 * u.arcsec),
        (-1.637 * u.arcsec, -0.969 * u.arcsec),
        (-3.73 * u.arcsec, -2.885 * u.arcsec),
        (-0.54 * u.arcsec, 0.889 * u.arcsec),
        (-3.257 * u.arcsec, -1.903 * u.arcsec),
        (-2.124 * u.arcsec, 1.827 * u.arcsec),
        (-2.146 * u.arcsec, -3.854 * u.arcsec),
        (0.0 * u.arcsec, 0.0 * u.arcsec),
        (-0.496 * u.arcsec, -2.836 * u.arcsec),
        (-96.08 * u.arcsec, 0.181 * u.arcsec),
        (-1.637 * u.arcsec, 2.73 * u.arcsec),
        (-1.075 * u.arcsec, 0.004 * u.arcsec),
        (-3.796 * u.arcsec, 2.739 * u.arcsec),
        (-1.695 * u.arcsec, 0.885 * u.arcsec),
        (-0.504 * u.arcsec, -4.783 * u.arcsec),
        (-1.084 * u.arcsec, -1.912 * u.arcsec),
        (-2.721 * u.arcsec, -2.845 * u.arcsec),
        (-93.973 * u.arcsec, 0.181 * u.arcsec),
        (-95.044 * u.arcsec, 0.186 * u.arcsec),
        (-0.473 * u.arcsec, 4.801 * u.arcsec),
        (-1.655 * u.arcsec, 4.73 * u.arcsec),
        (-4.35 * u.arcsec, -1.85 * u.arcsec),
        (-2.739 * u.arcsec, 2.774 * u.arcsec),
        (-1.062 * u.arcsec, 1.827 * u.arcsec),
        (-2.133 * u.arcsec, 0.004 * u.arcsec),
        (-4.385 * u.arcsec, 1.788 * u.arcsec),
        (-2.752 * u.arcsec, 0.894 * u.arcsec),
        (-97.681 * u.arcsec, 0.279 * u.arcsec),
        (-2.721 * u.arcsec, -0.965 * u.arcsec),
        (-3.332 * u.arcsec, 1.783 * u.arcsec),
        (-0.509 * u.arcsec, -0.965 * u.arcsec),
        (-3.265 * u.arcsec, 3.739 * u.arcsec),
        (-5.416 * u.arcsec, -0.044 * u.arcsec),
        (-2.009 * u.arcsec, 7.597 * u.arcsec),
        (-4.358 * u.arcsec, 3.743 * u.arcsec),
        (-3.252 * u.arcsec, 5.544 * u.arcsec),
        (-2.708 * u.arcsec, 6.704 * u.arcsec),
        (-1.633 * u.arcsec, -6.69 * u.arcsec),
        (-5.425 * u.arcsec, -1.956 * u.arcsec),
        (-89.252 * u.arcsec, -0.739 * u.arcsec),
        (-5.407 * u.arcsec, 1.792 * u.arcsec),
        (-0.903 * u.arcsec, 7.513 * u.arcsec),
        (-4.858 * u.arcsec, 2.748 * u.arcsec),
        (-2.686 * u.arcsec, -6.646 * u.arcsec),
        (-1.031 * u.arcsec, -5.761 * u.arcsec),
        (-2.664 * u.arcsec, -4.788 * u.arcsec),
        (-3.739 * u.arcsec, -4.788 * u.arcsec),
        (-1.447 * u.arcsec, 8.628 * u.arcsec),
        (-4.31 * u.arcsec, -3.858 * u.arcsec),
        (-2.102 * u.arcsec, -5.761 * u.arcsec),
        (-1.049 * u.arcsec, -7.58 * u.arcsec),
        (-94.345 * u.arcsec, -0.735 * u.arcsec),
        (-4.876 * u.arcsec, 0.841 * u.arcsec),
        (-0.522 * u.arcsec, -6.69 * u.arcsec),
        (-5.987 * u.arcsec, 0.845 * u.arcsec),
        (-1.606 * u.arcsec, 6.695 * u.arcsec),
        (0.022 * u.arcsec, -7.58 * u.arcsec),
        (-3.81 * u.arcsec, 4.646 * u.arcsec),
        (-6.482 * u.arcsec, -0.04 * u.arcsec),
        (-5.956 * u.arcsec, -1.013 * u.arcsec),
        (0.031 * u.arcsec, -5.677 * u.arcsec),
        (-2.142 * u.arcsec, 5.549 * u.arcsec),
        (-4.841 * u.arcsec, -2.934 * u.arcsec),
        (-0.681 * u.arcsec, 6.677 * u.arcsec),
        (-91.425 * u.arcsec, -0.712 * u.arcsec),
        (-1.624 * u.arcsec, -8.504 * u.arcsec),
        (-2.73 * u.arcsec, 4.659 * u.arcsec),
        (-3.204 * u.arcsec, -5.765 * u.arcsec),
        (-2.106 * u.arcsec, -7.562 * u.arcsec),
        (-4.841 * u.arcsec, -1.018 * u.arcsec),
        (-1.035 * u.arcsec, 5.673 * u.arcsec),
        (-4.881 * u.arcsec, 4.642 * u.arcsec),
        (-4.319 * u.arcsec, -5.765 * u.arcsec),
        (-4.31 * u.arcsec, 5.527 * u.arcsec),
        (-3.677 * u.arcsec, -6.642 * u.arcsec),
        (-4.805 * u.arcsec, -6.726 * u.arcsec),
        (-8.075 * u.arcsec, -1.009 * u.arcsec),
        (-5.991 * u.arcsec, 4.642 * u.arcsec),
        (-90.801 * u.arcsec, 0.177 * u.arcsec),
        (-4.854 * u.arcsec, -4.801 * u.arcsec),
        (-5.385 * u.arcsec, -3.867 * u.arcsec),
        (-6.531 * u.arcsec, 1.779 * u.arcsec),
        (-5.912 * u.arcsec, -2.805 * u.arcsec),
        (-6.491 * u.arcsec, 3.739 * u.arcsec),
        (-7.071 * u.arcsec, -1.018 * u.arcsec),
        (-5.907 * u.arcsec, -4.801 * u.arcsec),
        (-7.084 * u.arcsec, 2.743 * u.arcsec),
        (-8.106 * u.arcsec, 0.85 * u.arcsec),
        (-4.881 * u.arcsec, 6.58 * u.arcsec),
        (-8.69 * u.arcsec, -0.058 * u.arcsec),
        (-2.686 * u.arcsec, -8.478 * u.arcsec),
        (-88.681 * u.arcsec, 0.173 * u.arcsec),
        (-6.434 * u.arcsec, -3.858 * u.arcsec),
        (-5.327 * u.arcsec, -5.765 * u.arcsec),
        (-5.429 * u.arcsec, 3.708 * u.arcsec),
        (-3.124 * u.arcsec, 7.491 * u.arcsec),
        (-7.013 * u.arcsec, -2.885 * u.arcsec),
        (-7.646 * u.arcsec, 1.774 * u.arcsec),
        (-7.597 * u.arcsec, -1.96 * u.arcsec),
        (-3.204 * u.arcsec, -7.504 * u.arcsec),
        (-4.199 * u.arcsec, 7.544 * u.arcsec),
        (-3.819 * u.arcsec, 6.668 * u.arcsec),
        (-6.5 * u.arcsec, -1.96 * u.arcsec),
        (-7.049 * u.arcsec, 0.845 * u.arcsec),
        (-89.721 * u.arcsec, 0.181 * u.arcsec),
        (-5.978 * u.arcsec, 2.761 * u.arcsec),
        (-3.659 * u.arcsec, 8.456 * u.arcsec),
        (-5.367 * u.arcsec, 5.562 * u.arcsec),
        (-3.695 * u.arcsec, -8.5 * u.arcsec),
        (-2.531 * u.arcsec, 8.593 * u.arcsec),
        (-7.646 * u.arcsec, -0.053 * u.arcsec),
        (-4.274 * u.arcsec, -7.58 * u.arcsec),
        (-9.668 * u.arcsec, -1.973 * u.arcsec),
        (-5.341 * u.arcsec, -7.584 * u.arcsec),
        (-6.429 * u.arcsec, -7.5 * u.arcsec),
        (-6.513 * u.arcsec, 5.558 * u.arcsec),
        (-8.606 * u.arcsec, -3.903 * u.arcsec),
        (-7.527 * u.arcsec, -3.858 * u.arcsec),
        (-7.597 * u.arcsec, 3.739 * u.arcsec),
        (-87.611 * u.arcsec, 0.177 * u.arcsec),
        (-8.155 * u.arcsec, 2.748 * u.arcsec),
        (-9.19 * u.arcsec, -1.031 * u.arcsec),
        (-6.0 * u.arcsec, 6.54 * u.arcsec),
        (-7.013 * u.arcsec, -6.704 * u.arcsec),
        (-4.681 * u.arcsec, 8.54 * u.arcsec),
        (-7.022 * u.arcsec, -4.796 * u.arcsec),
        (-5.336 * u.arcsec, 7.469 * u.arcsec),
        (-9.226 * u.arcsec, 2.735 * u.arcsec),
        (-8.088 * u.arcsec, 4.619 * u.arcsec),
        (-5.85 * u.arcsec, -6.726 * u.arcsec),
        (-9.757 * u.arcsec, 1.783 * u.arcsec),
        (-7.035 * u.arcsec, 4.642 * u.arcsec),
        (-9.757 * u.arcsec, -0.053 * u.arcsec),
        (-4.814 * u.arcsec, -8.513 * u.arcsec),
        (-7.478 * u.arcsec, -5.761 * u.arcsec),
        (-5.73 * u.arcsec, 8.544 * u.arcsec),
        (-7.049 * u.arcsec, 6.553 * u.arcsec),
        (-7.611 * u.arcsec, 5.655 * u.arcsec),
        (-8.022 * u.arcsec, -4.792 * u.arcsec),
        (-8.726 * u.arcsec, 3.717 * u.arcsec),
        (-8.726 * u.arcsec, 1.792 * u.arcsec),
        (-10.765 * u.arcsec, -0.04 * u.arcsec),
        (-88.159 * u.arcsec, -0.752 * u.arcsec),
        (-9.071 * u.arcsec, -2.92 * u.arcsec),
        (-7.969 * u.arcsec, -2.796 * u.arcsec),
        (-10.195 * u.arcsec, -1.022 * u.arcsec),
        (-9.159 * u.arcsec, 0.845 * u.arcsec),
        (-8.646 * u.arcsec, -1.96 * u.arcsec),
        (-6.442 * u.arcsec, -5.77 * u.arcsec),
        (-5.854 * u.arcsec, -8.509 * u.arcsec),
        (-6.487 * u.arcsec, 7.588 * u.arcsec),
        (-10.23 * u.arcsec, 0.85 * u.arcsec),
        (6.513 * u.arcsec, 5.633 * u.arcsec),
        (10.288 * u.arcsec, -0.889 * u.arcsec),
        (90.004 * u.arcsec, -1.159 * u.arcsec),
        (8.142 * u.arcsec, -2.881 * u.arcsec),
        (9.243 * u.arcsec, -0.81 * u.arcsec),
        (7.611 * u.arcsec, 3.889 * u.arcsec),
        (7.0 * u.arcsec, -6.615 * u.arcsec),
        (7.633 * u.arcsec, -3.81 * u.arcsec),
        (5.42 * u.arcsec, -7.535 * u.arcsec),
        (8.146 * u.arcsec, -4.748 * u.arcsec),
        (9.212 * u.arcsec, -2.881 * u.arcsec),
        (8.642 * u.arcsec, -1.867 * u.arcsec),
        (8.549 * u.arcsec, 1.823 * u.arcsec),
        (4.916 * u.arcsec, -8.465 * u.arcsec),
        (5.942 * u.arcsec, -6.606 * u.arcsec),
        (8.066 * u.arcsec, 2.81 * u.arcsec),
        (6.527 * u.arcsec, -7.469 * u.arcsec),
        (91.053 * u.arcsec, -1.212 * u.arcsec),
        (8.726 * u.arcsec, -3.743 * u.arcsec),
        (6.075 * u.arcsec, 8.553 * u.arcsec),
        (7.04 * u.arcsec, -4.743 * u.arcsec),
        (9.173 * u.arcsec, 2.823 * u.arcsec),
        (7.58 * u.arcsec, -5.717 * u.arcsec),
        (5.942 * u.arcsec, 6.58 * u.arcsec),
        (6.513 * u.arcsec, -5.708 * u.arcsec),
        (9.699 * u.arcsec, -1.863 * u.arcsec),
        (7.035 * u.arcsec, 4.735 * u.arcsec),
        (8.177 * u.arcsec, 4.668 * u.arcsec),
        (5.456 * u.arcsec, 7.69 * u.arcsec),
        (90.597 * u.arcsec, -0.155 * u.arcsec),
        (9.177 * u.arcsec, 0.973 * u.arcsec),
        (8.695 * u.arcsec, 3.819 * u.arcsec),
        (10.814 * u.arcsec, 0.088 * u.arcsec),
        (7.566 * u.arcsec, 5.633 * u.arcsec),
        (5.027 * u.arcsec, 8.522 * u.arcsec),
        (9.699 * u.arcsec, 1.881 * u.arcsec),
        (10.283 * u.arcsec, 1.022 * u.arcsec),
        (6.951 * u.arcsec, 6.588 * u.arcsec),
        (3.239 * u.arcsec, -7.584 * u.arcsec),
        (7.08 * u.arcsec, -2.788 * u.arcsec),
        (6.611 * u.arcsec, 3.92 * u.arcsec),
        (5.987 * u.arcsec, -4.739 * u.arcsec),
        (3.867 * u.arcsec, 8.562 * u.arcsec),
        (5.425 * u.arcsec, -5.712 * u.arcsec),
        (8.071 * u.arcsec, 0.938 * u.arcsec),
        (4.296 * u.arcsec, -5.761 * u.arcsec),
        (6.522 * u.arcsec, -3.81 * u.arcsec),
        (7.013 * u.arcsec, 0.934 * u.arcsec),
        (8.655 * u.arcsec, 0.04 * u.arcsec),
        (5.473 * u.arcsec, 5.619 * u.arcsec),
        (6.0 * u.arcsec, -2.836 * u.arcsec),
        (7.588 * u.arcsec, -1.867 * u.arcsec),
        (5.451 * u.arcsec, 3.805 * u.arcsec),
        (6.416 * u.arcsec, 1.827 * u.arcsec),
        (7.058 * u.arcsec, -0.889 * u.arcsec),
        (5.925 * u.arcsec, 2.858 * u.arcsec),
        (4.841 * u.arcsec, 6.575 * u.arcsec),
        (7.487 * u.arcsec, 1.819 * u.arcsec),
        (2.774 * u.arcsec, 8.575 * u.arcsec),
        (4.323 * u.arcsec, 7.575 * u.arcsec),
        (8.173 * u.arcsec, -0.889 * u.arcsec),
        (6.969 * u.arcsec, 2.814 * u.arcsec),
        (6.071 * u.arcsec, 4.85 * u.arcsec),
        (6.522 * u.arcsec, -1.863 * u.arcsec),
        (4.398 * u.arcsec, 5.571 * u.arcsec),
        (3.863 * u.arcsec, -8.513 * u.arcsec),
        (4.929 * u.arcsec, -6.642 * u.arcsec),
        (4.296 * u.arcsec, -7.58 * u.arcsec),
        (7.597 * u.arcsec, 0.044 * u.arcsec),
        (3.279 * u.arcsec, 7.58 * u.arcsec),
        (3.757 * u.arcsec, 6.602 * u.arcsec),
        (2.735 * u.arcsec, -8.473 * u.arcsec),
        (91.637 * u.arcsec, -0.257 * u.arcsec),
        (3.823 * u.arcsec, -6.65 * u.arcsec),
        (92.659 * u.arcsec, -0.252 * u.arcsec),
        (5.416 * u.arcsec, -3.819 * u.arcsec),
        (4.965 * u.arcsec, -4.655 * u.arcsec),
        (93.752 * u.arcsec, -0.257 * u.arcsec),
        (1.659 * u.arcsec, 6.686 * u.arcsec),
        (4.358 * u.arcsec, -3.81 * u.arcsec),
        (0.606 * u.arcsec, -8.509 * u.arcsec),
        (0.677 * u.arcsec, 8.566 * u.arcsec),
        (1.075 * u.arcsec, 7.686 * u.arcsec),
        (0.15 * u.arcsec, 7.611 * u.arcsec),
        (1.633 * u.arcsec, 8.442 * u.arcsec),
        (94.296 * u.arcsec, -1.204 * u.arcsec),
        (4.332 * u.arcsec, 3.765 * u.arcsec),
        (4.85 * u.arcsec, 2.827 * u.arcsec),
        (2.23 * u.arcsec, 5.584 * u.arcsec),
        (2.204 * u.arcsec, -7.575 * u.arcsec),
        (1.146 * u.arcsec, 5.664 * u.arcsec),
        (0.558 * u.arcsec, 6.562 * u.arcsec),
        (0.597 * u.arcsec, -6.646 * u.arcsec),
        (3.894 * u.arcsec, 4.717 * u.arcsec),
        (2.814 * u.arcsec, 4.841 * u.arcsec),
        (2.681 * u.arcsec, -4.792 * u.arcsec),
        (2.137 * u.arcsec, 7.681 * u.arcsec),
        (5.314 * u.arcsec, 1.819 * u.arcsec),
        (2.735 * u.arcsec, 6.677 * u.arcsec),
        (2.235 * u.arcsec, -5.757 * u.arcsec),
        (0.04 * u.arcsec, 5.575 * u.arcsec),
        (3.872 * u.arcsec, -4.792 * u.arcsec),
        (5.434 * u.arcsec, 0.049 * u.arcsec),
        (4.938 * u.arcsec, -2.876 * u.arcsec),
        (1.084 * u.arcsec, -5.765 * u.arcsec),
        (4.903 * u.arcsec, -0.929 * u.arcsec),
        (5.42 * u.arcsec, -1.867 * u.arcsec),
        (6.588 * u.arcsec, -0.013 * u.arcsec),
        (6.009 * u.arcsec, -0.876 * u.arcsec),
        (4.912 * u.arcsec, 0.942 * u.arcsec),
        (93.235 * u.arcsec, -1.159 * u.arcsec),
        (3.301 * u.arcsec, -5.761 * u.arcsec),
        (3.323 * u.arcsec, 5.584 * u.arcsec),
        (5.978 * u.arcsec, 0.934 * u.arcsec),
        (1.668 * u.arcsec, -6.646 * u.arcsec),
        (1.659 * u.arcsec, -8.469 * u.arcsec),
        (1.097 * u.arcsec, -7.571 * u.arcsec),
        (2.743 * u.arcsec, -6.642 * u.arcsec),
        (1.062 * u.arcsec, 3.779 * u.arcsec),
        (3.274 * u.arcsec, 3.752 * u.arcsec),
        (2.19 * u.arcsec, 0.044 * u.arcsec),
        (2.226 * u.arcsec, -1.951 * u.arcsec),
        (0.615 * u.arcsec, -2.929 * u.arcsec),
        (2.743 * u.arcsec, -0.925 * u.arcsec),
        (2.217 * u.arcsec, 3.774 * u.arcsec),
        (1.646 * u.arcsec, 2.823 * u.arcsec),
        (96.956 * u.arcsec, -0.252 * u.arcsec),
        (-0.009 * u.arcsec, -3.85 * u.arcsec),
        (0.606 * u.arcsec, 4.624 * u.arcsec),
        (2.699 * u.arcsec, 0.934 * u.arcsec),
        (0.566 * u.arcsec, -0.969 * u.arcsec),
        (1.659 * u.arcsec, -2.934 * u.arcsec),
        (3.801 * u.arcsec, 2.819 * u.arcsec),
        (3.27 * u.arcsec, 0.049 * u.arcsec),
        (1.615 * u.arcsec, -4.792 * u.arcsec),
        (98.049 * u.arcsec, -0.204 * u.arcsec),
        (3.832 * u.arcsec, -2.925 * u.arcsec),
        (2.199 * u.arcsec, -3.81 * u.arcsec),
        (1.677 * u.arcsec, -0.92 * u.arcsec),
        (3.845 * u.arcsec, -0.929 * u.arcsec),
        (4.327 * u.arcsec, 0.044 * u.arcsec),
        (1.066 * u.arcsec, -3.863 * u.arcsec),
        (1.102 * u.arcsec, 1.823 * u.arcsec),
        (0.571 * u.arcsec, 0.889 * u.arcsec),
        (95.938 * u.arcsec, -0.252 * u.arcsec),
        (0.018 * u.arcsec, 3.73 * u.arcsec),
        (2.704 * u.arcsec, 2.823 * u.arcsec),
        (0.562 * u.arcsec, -4.788 * u.arcsec),
        (2.08 * u.arcsec, 1.77 * u.arcsec),
        (3.261 * u.arcsec, -3.858 * u.arcsec),
        (4.323 * u.arcsec, -1.907 * u.arcsec),
        (2.726 * u.arcsec, -2.929 * u.arcsec),
        (1.712 * u.arcsec, 4.832 * u.arcsec),
        (1.093 * u.arcsec, -1.907 * u.arcsec),
        (4.208 * u.arcsec, 1.774 * u.arcsec),
        (1.137 * u.arcsec, 0.049 * u.arcsec),
        (3.265 * u.arcsec, -1.907 * u.arcsec),
        (0.504 * u.arcsec, 2.708 * u.arcsec),
        (0.0 * u.arcsec, 1.832 * u.arcsec),
        (1.633 * u.arcsec, 0.889 * u.arcsec),
        (3.85 * u.arcsec, 0.938 * u.arcsec),
        (-0.013 * u.arcsec, -1.903 * u.arcsec),
        (3.15 * u.arcsec, 1.77 * u.arcsec),
    ]

    def __init__(self):
        self.fibers = _SMI.FIBER_CENTERS
        self.fiber_radius = _SMI.FIBER_RADIUS

    def bundle(
        self, center: SkyCoord, wcs: WCS, include_fibers: bool
    ) -> GroupAnnotation:
        return GroupAnnotation(
            [
                self._left_sky_bundle(center, wcs, include_fibers),
                self._target_bundle(center, wcs, include_fibers),
                self._right_sky_bundle(center, wcs, include_fibers),
            ]
        )

    def _left_sky_bundle_fibers(self):
        return [c for c in self.fibers if c[0] < -50 * u.arcsec]

    def _target_bundle_fibers(self):
        return [
            c
            for c in _SMI.FIBER_CENTERS
            if -12 * u.arcsec < c[0] < 12 * u.arcsec
            and -12 * u.arcsec < c[1] < 12 * u.arcsec
        ]

    def _right_sky_bundle_fibers(self):
        return [c for c in _SMI.FIBER_CENTERS if c[0] > 50 * u.arcsec]

    @staticmethod
    def _fiber_outlines(
        fibers: Iterable, radius: Angle, center: SkyCoord, wcs: WCS
    ) -> GroupAnnotation:
        annotation = GroupAnnotation([])
        for x, y in fibers:
            offset_x = x.to_value(u.arcsec)
            offset_y = y.to_value(u.arcsec)

            fiber_annotation = CircleAnnotation(
                center=center,
                radius=radius,
                wcs=wcs,
                edgecolor="red",
            ).translate((offset_x, offset_y) * u.arcsec)
            annotation.add_item(fiber_annotation)

        return annotation

    def _left_sky_bundle(self, center: SkyCoord, wcs: WCS, include_fibers: bool):
        annotation = GroupAnnotation([])
        if include_fibers:
            sky_bundle_fibers = self._left_sky_bundle_fibers()
            annotation.add_item(
                _SMI._fiber_outlines(sky_bundle_fibers, self.fiber_radius, center, wcs)
            )
        annotation.add_item(self._left_sky_bundle_outline(center, wcs))

        return annotation

    def _left_sky_bundle_outline(self, center: SkyCoord, wcs: WCS):
        # Approximate the outline by a polygon consisting of horizontal and vertical lines.
        r = self.fiber_radius
        left_sky_bundle_fibers = self._left_sky_bundle_fibers()
        upper_row_fibers = [c for c in left_sky_bundle_fibers if c[1] > 0 * u.arcsec]
        lower_row_fibers = [c for c in left_sky_bundle_fibers if c[1] < 0 * u.arcsec]
        upper_row_xmin = min(c[0] for c in upper_row_fibers) - r
        upper_row_xmax = max(c[0] for c in upper_row_fibers) + r
        upper_row_ymin = min(c[1] for c in upper_row_fibers) - r
        upper_row_ymax = max(c[1] for c in upper_row_fibers) + r
        lower_row_xmin = min(c[0] for c in lower_row_fibers) - r
        lower_row_ymin = min(c[1] for c in lower_row_fibers) - r
        vertices = [
            (upper_row_xmin, upper_row_ymax),
            (upper_row_xmax, upper_row_ymax),
            (upper_row_xmax, lower_row_ymin),
            (lower_row_xmin, lower_row_ymin),
            (lower_row_xmin, upper_row_ymin),
            (upper_row_xmin, upper_row_ymin),
            (upper_row_xmin, upper_row_ymax),
        ]

        # So far we have assumed the bundle center to be at (0, 0), so we have to correct
        # this.
        vertex_positions = [translate(center, v) for v in vertices]

        line = LinePathAnnotation(
            vertex_positions, wcs=wcs, closed=True, edgecolor=_OUTLINE_COLOR
        )
        return line

    def _right_sky_bundle(self, center: SkyCoord, wcs: WCS, include_fibers: bool):
        annotation = GroupAnnotation([])
        if include_fibers:
            sky_bundle_fibers = self._right_sky_bundle_fibers()
            annotation.add_item(
                _SMI._fiber_outlines(sky_bundle_fibers, self.fiber_radius, center, wcs)
            )
        annotation.add_item(self._right_sky_bundle_outline(center, wcs))

        return annotation

    def _right_sky_bundle_outline(self, center: SkyCoord, wcs: WCS):
        # Approximate the outline by a polygon consisting of horizontal and vertical lines.
        r = self.fiber_radius
        right_sky_bundle_fibers = self._right_sky_bundle_fibers()
        upper_row_fibers = [
            c for c in right_sky_bundle_fibers if c[1] > -0.5 * u.arcsec
        ]
        lower_row_fibers = [
            c for c in right_sky_bundle_fibers if c[1] < -0.5 * u.arcsec
        ]
        upper_row_xmin = min(c[0] for c in upper_row_fibers) - r
        upper_row_xmax = max(c[0] for c in upper_row_fibers) + r
        upper_row_ymin = min(c[1] for c in upper_row_fibers) - r
        upper_row_ymax = max(c[1] for c in upper_row_fibers) + r
        lower_row_xmin = min(c[0] for c in lower_row_fibers) - r
        lower_row_ymin = min(c[1] for c in lower_row_fibers) - r
        vertices = [
            (upper_row_xmin, upper_row_ymax),
            (upper_row_xmax, upper_row_ymax),
            (upper_row_xmax, lower_row_ymin),
            (lower_row_xmin, lower_row_ymin),
            (lower_row_xmin, upper_row_ymin),
            (upper_row_xmin, upper_row_ymin),
            (upper_row_xmin, upper_row_ymax),
        ]

        # So far we have assumed the bundle center to be at (0, 0), so we have to correct
        # this.
        vertex_positions = [translate(center, v) for v in vertices]

        line = LinePathAnnotation(
            vertex_positions, wcs=wcs, closed=True, edgecolor=_OUTLINE_COLOR
        )
        return line

    def _target_bundle(self, center: SkyCoord, wcs: WCS, include_fibers: bool):
        annotation = GroupAnnotation([])
        if include_fibers:
            target_bundle_fibers = self._target_bundle_fibers()
            annotation.add_item(
                _SMI._fiber_outlines(
                    target_bundle_fibers, self.fiber_radius, center, wcs
                )
            )
        annotation.add_item(self._target_bundle_outline(center, wcs))

        return annotation

    def _target_bundle_outline(self, center: SkyCoord, wcs: WCS):
        # Approximate the outline by a vertically and horizontally symmetric hexagon.
        r = self.fiber_radius
        target_bundle_fibers = self._target_bundle_fibers()
        upper_row_xmin = (
            min(c[0] for c in target_bundle_fibers if c[1] > 8.5 * u.arcsec)
            - r
            - 0.25 * r
        )
        upper_row_xmax = -upper_row_xmin
        upper_row_ymax = max(c[1] for c in target_bundle_fibers) + r
        lower_row_xmin = upper_row_xmin
        lower_row_xmax = upper_row_xmax
        lower_row_ymin = min(c[1] for c in target_bundle_fibers) - r
        center_row_xmin = min(c[0] for c in target_bundle_fibers) - r - 0.25 * r
        center_row_xmax = -center_row_xmin
        center_row_y = (upper_row_ymax + lower_row_ymin) / 2

        vertices = [
            (upper_row_xmin, upper_row_ymax),
            (upper_row_xmax, upper_row_ymax),
            (center_row_xmax, center_row_y),
            (lower_row_xmax, lower_row_ymin),
            (lower_row_xmin, lower_row_ymin),
            (center_row_xmin, center_row_y),
            (upper_row_xmin, upper_row_ymax),
        ]

        # So far we have assumed the bundle center to be at (0, 0), so we have to correct
        # this.
        vertex_positions = [translate(center, v) for v in vertices]

        line = LinePathAnnotation(
            vertex_positions, wcs=wcs, closed=True, edgecolor=_OUTLINE_COLOR
        )
        return line
