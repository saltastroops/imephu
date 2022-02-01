from dataclasses import dataclass

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import GroupAnnotation
from imephu.finder_chart import FinderChart
from imephu.salt.annotation import rss, telescope
from imephu.salt.utils import MosMask

_FINDER_CHART_SIZE = 10 * u.arcmin


@dataclass
class MagnitudeRange:
    """A magnitude range.

    Attributes
    ----------
    bandpass: `str`
        The bandpass for which the magnitudes are given.
    min_magnitude: `float`
        The minimum (brightest) magnitude.
    max_magnitude: `float`
        The maximum (faintest) magnitude.
    """

    bandpass: str
    min_magnitude: float
    max_magnitude: float


@dataclass
class Target:
    """Target properties.

    Attributes
    ----------
    name: `str`
        The target name.
    position: `~astropy.coordinates.SkyCoord`
        The target position, as a right ascension and declination. This is taken to be
        the center of the finder chart.
    magnitude_range: `MagnitudeRange`
        The magnitude range of the target.
    """

    name: str
    position: SkyCoord
    magnitude_range: MagnitudeRange


@dataclass
class GeneralProperties:
    """Properties which are not specific to a particular instrument.

    Attributes
    ----------
    target: `Target`
        The target to find using this finder chart. Its position is taken to be the
        center of the finder chart.
    position_angle: `~astropy.coordinates.Angle`
        The position angle, as angle on the sky measured from north to east.
    automated_position_angle: `bool`
        Whether the position angle has been calculated automatically.
    proposal_code: `str`
        The proposal code.
    pi_family_name: `str`
        The family name of the Principal Investigator.
    survey: `str`
        The name of the survey from which the finder chart image was taken.
    """

    target: Target
    position_angle: Angle
    automated_position_angle: bool
    proposal_code: str
    pi_family_name: str
    survey: str


def salticam_finder_chart(general: GeneralProperties, is_slot_mode: bool = False) -> FinderChart:
    """Return the finder chart for a Salticam observation.

    Parameters
    ----------
    general: `GeneralProperties`
        Properties which are not specific to the instrument.
    is_slot_mode: `bool`, default: False
        Whether the observation is a slot mode one.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for a Salticam observation.
    """
    finder_chart = FinderChart.from_survey(
        general.survey, general.target.position, _FINDER_CHART_SIZE
    )
    annotation = _salticam_observation_annotation(
        general=general, is_slot_mode=is_slot_mode, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def rss_imaging_finder_chart(
    general: GeneralProperties, is_slot_mode: bool = False
) -> FinderChart:
    """Return the finder chart for an RSS imaging observation.

    Parameters
    ----------
    general: `GeneralProperties`
        Properties which are not specific to the instrument.
    is_slot_mode: `bool`, default: False
        Whether the observation is a slot mode one.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for an RSS imaging observation.
    """
    finder_chart = FinderChart.from_survey(
        general.survey, general.target.position, _FINDER_CHART_SIZE
    )
    annotation = _rss_imaging_observation_annotation(
        general=general, is_slot_mode=is_slot_mode, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def rss_longslit_finder_chart(
    general: GeneralProperties, slit_width: Angle, slit_height: Angle
) -> FinderChart:
    """Return the finder chart for an RSS longslit observation.

    Parameters
    ----------
    general: `GeneralProperties`
        Properties which are not specific to the instrument.
    slit_width: `~astropy.coordinates.Angle`
        The slit width, as an angle on the sky.
    slit_height: `~astropy.coordinates.Angle`
        The slit height, as an angle on the sky.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for an RSS longslit observation.
    """
    finder_chart = FinderChart.from_survey(
        general.survey, general.target.position, _FINDER_CHART_SIZE
    )
    annotation = _rss_longslit_observation_annotation(
        general=general,
        slit_width=slit_width,
        slit_height=slit_height,
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def rss_mos_finder_chart(
    general: GeneralProperties,
    mos_mask: MosMask,
    reference_star_box_width: Angle = Angle(5 * u.arcsec),
) -> FinderChart:
    """Return the finder chart for an RSS MOS observation.

    Parameters
    ----------
    general: `GeneralProperties`
        Properties which are not specific to the instrument.
    mos_mask: `~imephu.salt.utils.MosMask`
        The MOS mask.
    reference_star_box_width: `~astropy.coordinates.Angle`, default: 5 arcseconds
        The width (and height) of the boxes around reference stars, as an angle on the
        sky.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for an RSS MOS observation.
    """
    finder_chart = FinderChart.from_survey(
        general.survey, general.target.position, _FINDER_CHART_SIZE
    )
    annotation = _rss_mos_observation_annotation(
        general=general,
        mos_mask=mos_mask,
        reference_star_box_width=reference_star_box_width,
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def rss_fabry_perot_finder_chart(
    general: GeneralProperties,
) -> FinderChart:
    """Return the finder chart for an RSS Fabry-Pérot observation.

    Parameters
    ----------
    general: `GeneralProperties`
        Properties which are not specific to the instrument.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The annotation for an RSS Fabry-Pérot observation.
    """
    finder_chart = FinderChart.from_survey(
        general.survey, general.target.position, _FINDER_CHART_SIZE
    )
    annotation = _rss_fabry_perot_observation_annotation(
        general=general, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def hrs_finder_chart(general: GeneralProperties) -> FinderChart:
    """Return the annotation for an HRS observation.

    Parameters
    ----------
    general: `GeneralProperties`
        Properties which are not specific to the instrument.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for an HRS observation.
    """
    finder_chart = FinderChart.from_survey(
        general.survey, general.target.position, _FINDER_CHART_SIZE
    )
    annotation = _hrs_observation_annotation(general=general, wcs=finder_chart.wcs)
    finder_chart.add_annotation(annotation)
    return finder_chart


def _salticam_observation_annotation(
    general: GeneralProperties, is_slot_mode: bool, wcs: WCS
) -> GroupAnnotation:
    return _imaging_annotation(general, is_slot_mode, wcs)


def _rss_imaging_observation_annotation(
    general: GeneralProperties, is_slot_mode: bool, wcs: WCS
) -> GroupAnnotation:
    return _imaging_annotation(general, is_slot_mode, wcs)


def _rss_longslit_observation_annotation(
    general: GeneralProperties, slit_width: Angle, slit_height: Angle, wcs: WCS
) -> GroupAnnotation:
    observation_annotation = _base_annotations(general, wcs)
    longslit_annotation = rss.longslit_annotation(
        fits_center=general.target.position,
        slit_width=slit_width,
        slit_height=slit_height,
        position_angle=general.position_angle,
        wcs=wcs,
    )
    observation_annotation.add_item(longslit_annotation)
    return observation_annotation


def _rss_mos_observation_annotation(
    general: GeneralProperties,
    mos_mask: MosMask,
    reference_star_box_width: Angle,
    wcs: WCS,
) -> GroupAnnotation:
    observation_annotation = _base_annotations(general, wcs)
    mask_annotation = rss.mos_mask_annotation(
        mos_mask=mos_mask,
        wcs=wcs,
        reference_star_box_width=reference_star_box_width,
    )
    observation_annotation.add_item(mask_annotation)
    return observation_annotation


def _rss_fabry_perot_observation_annotation(
    general: GeneralProperties, wcs: WCS
) -> GroupAnnotation:
    return _imaging_annotation(general=general, is_slot_mode=False, wcs=wcs)


def _hrs_observation_annotation(
    general: GeneralProperties, wcs: WCS
) -> GroupAnnotation:
    return _imaging_annotation(general=general, is_slot_mode=False, wcs=wcs)


def _imaging_annotation(
    general: GeneralProperties, is_slot_mode: bool, wcs: WCS
) -> GroupAnnotation:
    observation_annotation = _base_annotations(general, wcs)
    magnitude_range = general.target.magnitude_range
    magnitude_annotation = telescope.magnitude_range_annotation(
        bandpass=magnitude_range.bandpass,
        min_magnitude=magnitude_range.min_magnitude,
        max_magnitude=magnitude_range.max_magnitude,
        fits_center=general.target.position,
        wcs=wcs,
    )
    observation_annotation.add_item(magnitude_annotation)
    if is_slot_mode:
        center = general.target.position
        slot_annotation = telescope.slot_annotation(center, general.position_angle, wcs)
        observation_annotation.add_item(slot_annotation)
    return observation_annotation


def _base_annotations(general: GeneralProperties, wcs: WCS) -> GroupAnnotation:
    return telescope.base_annotations(
        target=general.target.name,
        proposal_code=general.proposal_code,
        pi_family_name=general.pi_family_name,
        position_angle=general.position_angle,
        automated_position_angle=general.automated_position_angle,
        survey=general.survey,
        fits_center=general.target.position,
        wcs=wcs,
    )