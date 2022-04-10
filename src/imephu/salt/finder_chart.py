from __future__ import annotations

import dataclasses
import os
from dataclasses import dataclass
from datetime import datetime
from typing import Any, BinaryIO, Generator, List, Optional, Tuple, Union

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS

from imephu.annotation.general import GroupAnnotation
from imephu.annotation.motion import motion_annotation
from imephu.finder_chart import FinderChart
from imephu.salt.annotation import nir, rss, telescope
from imephu.salt.utils import MosMask
from imephu.utils import (
    Ephemeris,
    MagnitudeRange,
    ephemerides_magnitude_range,
    mid_position,
)

_FINDER_CHART_SIZE = 10 * u.arcmin


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
    magnitude_range: `imephu.utils.MagnitudeRange`, optional
        The magnitude range of the target.
    """

    name: str
    position: SkyCoord
    magnitude_range: Optional[MagnitudeRange]


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
    survey: `str`, optional
        The image survey from which the FITS image is taken.
    """

    target: Target
    position_angle: Angle
    automated_position_angle: bool
    proposal_code: str
    pi_family_name: str
    survey: Optional[str] = None


def salticam_finder_chart(
    fits: Union[str, BinaryIO, os.PathLike[Any]],
    general: GeneralProperties,
    is_slot_mode: bool = False,
) -> FinderChart:
    """Return the finder chart for a Salticam observation.

    Parameters
    ----------
    fits: `str`, `path-like` or `binary file-like`
        FITS file to display.
    general: `GeneralProperties`
        Properties which are not specific to the instrument.
    is_slot_mode: `bool`, default: False
        Whether the observation is a slot mode one.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for a Salticam observation.
    """
    finder_chart = FinderChart(fits)
    annotation = _salticam_observation_annotation(
        general=general, is_slot_mode=is_slot_mode, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def rss_imaging_finder_chart(
    fits: Union[str, BinaryIO, os.PathLike[Any]],
    general: GeneralProperties,
    is_slot_mode: bool = False,
) -> FinderChart:
    """Return the finder chart for an RSS imaging observation.

    Parameters
    ----------
    fits: `str`, `path-like` or `binary file-like`
        FITS file to display.
    general: `GeneralProperties`
        Properties which are not specific to the instrument.
    is_slot_mode: `bool`, default: False
        Whether the observation is a slot mode one.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for an RSS imaging observation.
    """
    finder_chart = FinderChart(fits)
    annotation = _rss_imaging_observation_annotation(
        general=general, is_slot_mode=is_slot_mode, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def rss_longslit_finder_chart(
    fits: Union[str, BinaryIO, os.PathLike[Any]],
    general: GeneralProperties,
    slit_width: Angle,
    slit_height: Angle,
) -> FinderChart:
    """Return the finder chart for an RSS longslit observation.

    Parameters
    ----------
    fits: `str`, `path-like` or `binary file-like`
        FITS file to display.
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
    finder_chart = FinderChart(fits)
    annotation = _rss_longslit_observation_annotation(
        general=general,
        slit_width=slit_width,
        slit_height=slit_height,
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def rss_mos_finder_chart(
    fits: Union[str, BinaryIO, os.PathLike[Any]],
    general: GeneralProperties,
    mos_mask: MosMask,
    reference_star_box_width: Angle = Angle(5 * u.arcsec),
) -> FinderChart:
    """Return the finder chart for an RSS MOS observation.

    Parameters
    ----------
    fits: `str`, `path-like` or `binary file-like`
        FITS file to display.
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
    finder_chart = FinderChart(fits)
    annotation = _rss_mos_observation_annotation(
        general=general,
        mos_mask=mos_mask,
        reference_star_box_width=reference_star_box_width,
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def rss_fabry_perot_finder_chart(
    fits: Union[str, BinaryIO, os.PathLike[Any]],
    general: GeneralProperties,
) -> FinderChart:
    """Return the finder chart for an RSS Fabry-Pérot observation.

    Parameters
    ----------
    fits: `str`, `path-like` or `binary file-like`
        FITS file to display.
    general: `GeneralProperties`
        Properties which are not specific to the instrument.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for an RSS Fabry-Pérot observation.
    """
    finder_chart = FinderChart(fits)
    annotation = _rss_fabry_perot_observation_annotation(
        general=general, wcs=finder_chart.wcs
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def nir_finder_chart(
    fits: Union[str, BinaryIO, os.PathLike[Any]],
    general: GeneralProperties,
    science_bundle_center: SkyCoord,
    bundle_separation: Angle,
) -> FinderChart:
    """Return the finder chart for an NIR observation.

    Parameters
    ----------
    fits: `str`, `path-like` or `binary file-like`
        FITS file to display.
    science_bundle_center: `~astropy.coordinates.SkyCoord`
        The center position of the sky fiber bundle, as a position on the sky, in
        right ascension and declination.
    bundle_separation: ~astropy.coordinates.Angle`
        The separation between the science fiber bundle and the sky fiber bundles, as an
        angle on the sky. The separation is measured between the center of the science
        bundle and the midpoint of the line between the centers of the sky bundles.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for an NIR observation.
    """
    finder_chart = FinderChart(fits)
    annotation = _nir_observation_annotation(
        general=general,
        science_bundle_center=science_bundle_center,
        bundle_separation=bundle_separation,
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(annotation)
    return finder_chart


def hrs_finder_chart(
    fits: Union[str, BinaryIO, os.PathLike[Any]], general: GeneralProperties
) -> FinderChart:
    """Return the finder chart for an HRS observation.

    Parameters
    ----------
    fits: `str`, `path-like` or `binary file-like`
        FITS file to display.
    general: `GeneralProperties`
        Properties which are not specific to the instrument.

    Returns
    -------
    `~imephu.finder_chart.FinderChart`
        The finder chart for an HRS observation.
    """
    finder_chart = FinderChart(fits)
    annotation = _hrs_observation_annotation(general=general, wcs=finder_chart.wcs)
    finder_chart.add_annotation(annotation)
    return finder_chart


def moving_target_finder_charts(
    general: GeneralProperties,
    start: datetime,
    end: datetime,
    ephemerides: List[Ephemeris],
    survey: str,
) -> Generator[Tuple[FinderChart, Tuple[datetime, datetime]], None, None]:
    """Return a generator for the finder charts for a moving target in a time interval.

    Parameters
    ----------
    general: `GeneralProperties`
        Properties which are not specific to the instrument.
    start: `~datetime.datetime`
        Start time of the time interval for which finder charts are generated. This must
        be timezone-aware.
    end: `~datetime.datetime`
        End time of the time interval for which finder charts are generated. This must
        be timezone-aware.
    ephemerides: list of `~imephu.utils.Ephemeris`
        The ephemerides for the moving target.
    survey: `str`
        The name of the survey from which the finder chart image should be taken.

    Yields
    ------
    `~imephu.finder_chart.FinderChart`
        The finder chart.
    """

    def _create_finder_chart(ephemerides_: List[Ephemeris]) -> FinderChart:
        midpoint = mid_position(ephemerides_[0].position, ephemerides_[-1].position)
        magnitude_range = ephemerides_magnitude_range(ephemerides_)
        general_ = dataclasses.replace(
            general,
            target=Target(
                name=general.target.name,
                position=midpoint,
                magnitude_range=magnitude_range,
            ),
        )
        finder_chart = FinderChart.from_survey(survey, midpoint, _FINDER_CHART_SIZE)
        annotation = _non_sidereal_annotation(
            general=general_, ephemerides=ephemerides_, wcs=finder_chart.wcs
        )
        finder_chart.add_annotation(annotation)
        return finder_chart

    return FinderChart.for_time_interval(
        start=start,
        end=end,
        ephemerides=ephemerides,
        max_track_length=0.8 * _FINDER_CHART_SIZE,
        create_finder_chart=_create_finder_chart,
    )


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


def _nir_observation_annotation(
    general: GeneralProperties,
    science_bundle_center: SkyCoord,
    bundle_separation: Angle,
    wcs: WCS,
) -> GroupAnnotation:
    observation_annotation = _base_annotations(general, wcs)
    bundles_annotation = nir.bundles_annotation(
        science_bundle_center=science_bundle_center,
        bundle_separation=bundle_separation,
        position_angle=general.position_angle,
        wcs=wcs,
    )
    observation_annotation.add_item(bundles_annotation)
    return observation_annotation


def _hrs_observation_annotation(
    general: GeneralProperties, wcs: WCS
) -> GroupAnnotation:
    return _imaging_annotation(general=general, is_slot_mode=False, wcs=wcs)


def _non_sidereal_annotation(
    general: GeneralProperties, ephemerides: List[Ephemeris], wcs: WCS
) -> GroupAnnotation:
    imaging_annotation = _imaging_annotation(
        general=general, is_slot_mode=False, wcs=wcs
    )
    track_annotation = motion_annotation(ephemerides=ephemerides, wcs=wcs)
    return GroupAnnotation(items=[imaging_annotation, track_annotation])


def _imaging_annotation(
    general: GeneralProperties, is_slot_mode: bool, wcs: WCS
) -> GroupAnnotation:
    observation_annotation = _base_annotations(general, wcs)
    magnitude_range = general.target.magnitude_range
    if magnitude_range:
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
