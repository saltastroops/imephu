from typing import Optional

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import (
    CrosshairsAnnotation,
    EmptyAnnotation,
    GroupAnnotation,
    RectangleAnnotation,
    TextAnnotation,
)
from imephu.salt.annotation import rss, salticam


def title_annotation(
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
        The title annotation.
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
        clip_on=False,
    )


def directions_annotation(fits_center: SkyCoord, wcs: WCS) -> GroupAnnotation:
    """Return an annotation for the east and north direction.

    Parameters
    ----------
    fits_center: `~astropy.coordinates.SkyCoord`
        The central position of the finder chart, in right ascension and declination.
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
        color=(0, 0.5, 1),
        horizontalalignment="right",
        style="italic",
        weight="bold",
        size="large",
    ).translate(Angle((4.8, 0) * u.arcmin))
    north_annotation = TextAnnotation(
        fits_center,
        "N",
        wcs=wcs,
        color=(0, 0.5, 1),
        style="italic",
        weight="bold",
        size="large",
    ).translate(Angle((0, 4.8) * u.arcmin))
    return GroupAnnotation([east_annotation, north_annotation])


def position_angle_annotation(
    position_angle: Angle, automated_position_angle: bool, wcs: WCS
) -> TextAnnotation:
    """Return a text annotation with the position angle.

    The text of the annotation is "PA: angle (auto)", where angle is a value in degrees,
    given with one fractional digit, and where "(auto)" is only included if the
    ``automated`` parameter is ``True``.

    Parameters
    ----------
    position_angle: `~astropy.coordinates.Angle`
        The position angle, as an angle on the sky from north to east.
    automated_position_angle: bool
        Whether the position angle has been calculated automatically.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.TextAnnotation`
        The annotation for the position angle.
    """
    text = f"PA = {position_angle.to_value(u.degree):.1f}"
    if automated_position_angle:
        text += " (auto)"
    return TextAnnotation(
        (1, -0.068),
        text,
        wcs=wcs,
        color="black",
        horizontalalignment="right",
        verticalalignment="baseline",
        style="italic",
        weight="bold",
        size="large",
    )


def survey_annotation(survey: str, wcs: WCS) -> TextAnnotation:
    """Return a text annotation with the survey name.

    Parameters
    ----------
    survey: str
        The survey name.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.TextAnnotation`
        The annotation with the survey name.
    """
    return TextAnnotation(
        (0, -0.068),
        survey,
        wcs=wcs,
        color="black",
        horizontalalignment="left",
        verticalalignment="baseline",
        style="italic",
        weight="bold",
        size="large",
    )


def magnitude_range_annotation(
    bandpass: str,
    min_magnitude: float,
    max_magnitude: float,
    fits_center: SkyCoord,
    wcs: WCS,
) -> TextAnnotation:
    """Return a text annotation with the magnitude range.

    The magnitude range is given as a string of the form "bandpass = min magnitude -
    max magnitude", where the magnitudes have one fractional digit. If the difference
    between the minimum and maximum magnitude is less than 0.1, only the maximum
    magnitude is included.

    Parameters
    ----------
    bandpass: `str`
        Bandpass for which the magnitude range is given.
    min_magnitude: `float`
        Minimum (brightest) magnitude.
    max_magnitude: `float`
        Maximum (faintest) magnitude.
    fits_center: `~astropy.coordinates.SkyCoord`
        Center position of the finder chart, as a right ascension and declination.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.TextAnnotation`
        The magnitude range annotation.
    """
    if abs(max_magnitude - min_magnitude) > 0.09999:
        text = f"{bandpass} = {min_magnitude:.1f} - {max_magnitude:.1f}"
    else:
        text = f"{bandpass} = {max_magnitude:.1f}"
    return TextAnnotation(
        fits_center,
        text,
        wcs,
        color=(0, 0.5, 1),
        style="italic",
        weight="bold",
        size="large",
        horizontalalignment="center",
        verticalalignment="baseline",
    ).translate(Angle((0, -4.8) * u.arcmin))


def slot_annotation(
    center: SkyCoord, position_angle: Angle, wcs: WCS
) -> RectangleAnnotation:
    """Return an annotation showing the slot.

    Parameters
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The center of the slot on the sky, in right ascension and declination.
    position_angle: `~astropy.coordinates.Angle`
        The position angle as an angle on the sky, measured from north to east.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.RectangleAnnotation`
        The slot annotation.
    """
    return RectangleAnnotation(
        center=center,
        width=1 * u.arcmin / 3,
        height=10 * u.arcmin,
        wcs=wcs,
        edgecolor="red",
        alpha=0.5,
        linewidth=2,
    ).rotate(center, position_angle + 90 * u.deg)


def base_annotations(
    target: str,
    proposal_code: str,
    pi_family_name: str,
    position_angle: Angle,
    automated_position_angle: bool,
    survey: Optional[str],
    fits_center: SkyCoord,
    wcs: WCS,
) -> GroupAnnotation:
    """Return a group containing all the base annotations for SALT finder charts.

    The base annotations include the title, directions, position angle, survey name and
    the RSS and Salticam fields of view.

    Parameters
    ----------
    target: str
        The name of the target for which the finder chart is intended.
    proposal_code: str
        The proposal code for which the finder chart is intended.
    pi_family_name: str
        The family name of the Principal Investigator.
    position_angle: `~astropy.coordinates.Angle`
        The position angle, as an angle on the sky from north to east.
    automated_position_angle: bool
        Whether the position angle has been calculated automatically.
    survey: str, optional
        The image survey name.
    fits_center: `~astropy.coordinates.SkyCoord`
        The central position of the finder chart, in right ascension and declination
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.GroupAnnotation`
        An annotation with all the SALT base annotations.
    """
    return GroupAnnotation(
        [
            title_annotation(
                target=target,
                proposal_code=proposal_code,
                pi_family_name=pi_family_name,
                wcs=wcs,
            ),
            survey_annotation(survey=survey, wcs=wcs) if survey else EmptyAnnotation(),
            position_angle_annotation(
                position_angle=position_angle,
                automated_position_angle=automated_position_angle,
                wcs=wcs,
            ),
            directions_annotation(fits_center=fits_center, wcs=wcs),
            _crosshairs_annotation(fits_center=fits_center, wcs=wcs),
            salticam.field_of_view_annotation(fits_center=fits_center, wcs=wcs),
            rss.field_of_view_annotation(fits_center=fits_center, wcs=wcs),
        ]
    )


def _crosshairs_annotation(fits_center: SkyCoord, wcs: WCS) -> CrosshairsAnnotation:
    return CrosshairsAnnotation(
        center=fits_center, size=8 * u.arcmin, wcs=wcs, color="green"
    )
