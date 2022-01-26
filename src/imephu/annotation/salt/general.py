from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import GroupAnnotation, TextAnnotation

BLUE_COLOR = (0, 0.5, 1)


def title(
    target: str, proposal_code: str, pi_family_name: str, wcs: WCS
) -> TextAnnotation:
    """Return a text annotation with the title to the finder chart.

    The title is of the form "target name (proposal code; surname of the Principal
    Investigator)", such as "Magrathea (2022-1-SCI-042; Adams)".

    Parameters
    ----------
    target: str
        The name of the target for which the finder chart is intended.
    proposal_code: str
        The proposal code for which the finder chart is intended.
    pi_family_name: str
        The family name of the Principal Investigator.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.TextAnnotation`
        The titke annotation.
    """
    title_text = f"{target} ({proposal_code}; {pi_family_name})"
    return TextAnnotation(
        (0.5, 1.03),
        title_text,
        wcs=wcs,
        color="black",
        horizontalalignment="center",
        style="italic",
        weight="bold",
        size="large",
    )


def directions(fits_center: SkyCoord, wcs: WCS) -> GroupAnnotation:
    """Return an annotation for the east and north direction.

    Parameters
    ----------
    fits_center: `~astropy.coordinates.SkyCoord`
        The central position of the finder chart, in right ascension and declination
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~astropy.annotation.general.GroupAnnotation`
        The annotation for the east and north direction.
    """
    east_annotation = TextAnnotation(
        fits_center,
        "E",
        wcs=wcs,
        color=BLUE_COLOR,
        horizontalalignment="right",
        style="italic",
        weight="bold",
        size="large",
    ).translate(Angle((4.8, 0) * u.arcmin))
    north_annotation = TextAnnotation(
        fits_center,
        "N",
        wcs=wcs,
        color=BLUE_COLOR,
        style="italic",
        weight="bold",
        size="large",
    ).translate(Angle((0, 4.8) * u.arcmin))
    return GroupAnnotation([east_annotation, north_annotation])


def position_angle(angle: Angle, automated: bool, wcs: WCS) -> TextAnnotation:
    """Return a text annotation with the position angle.

    The text of the annotation is "PA: angle (auto)", where angle is a value in degrees,
    given with one fractional digit, and where "(auto)" is only included if the
    ``automated`` parameter is ``True``.

    Parameters
    ----------
    angle: `~astropy.coordinates.Angle`
        The position angle, as an angle on the sky from north to east.
    automated: bool
        Whether the position angle has been calculated automatically.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~immephu.annotation.general.TextAnnotation`
        The annotation for the position angle.
    """
    text = f"PA = {angle.to_value(u.degree):.1f}"
    if automated:
        text += " (auto)"
    return TextAnnotation(
        (1, -0.06),
        text,
        wcs=wcs,
        color="black",
        horizontalalignment="right",
        style="italic",
        weight="bold",
        size="large",
    )
