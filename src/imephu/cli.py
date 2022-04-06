import io
import json
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, BinaryIO, Dict, Generator, List, Optional, Tuple, Union, cast

import typer
import yaml
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from dateutil.parser import isoparse
from jsonschema import validate

import imephu
import imephu.salt.finder_chart as sfc
from imephu.annotation.motion import motion_annotation
from imephu.finder_chart import FinderChart
from imephu.salt.finder_chart import GeneralProperties, Target
from imephu.salt.utils import MosMask
from imephu.service.horizons import HorizonsService
from imephu.service.survey import load_fits
from imephu.utils import (
    Ephemeris,
    MagnitudeRange,
    ephemerides_magnitude_range,
    mid_position,
)

app = typer.Typer()


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"imephu {imephu.__version__}")
        raise typer.Exit()


@app.command()
def main(
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        file_okay=True,
        exists=True,
        resolve_path=True,
        help="Configuration details for the finder chart(s).",
    ),
    output: Optional[Path] = typer.Option(
        None, "--out", "-o", file_okay=True, resolve_path=True, help="Output file."
    ),
    format: str = typer.Option(
        "png",
        "--format",
        "-f",
        help="Image format (such as png or pdf) for the finder chart(s).",
    ),
    version: Optional[bool] = typer.Option(
        None, "--version", callback=_version_callback, help="Show the version and exit."
    ),
) -> None:
    """A tool for creating finder charts."""  # noqa: D401
    try:
        configuration = _read_configuration(config)
        _validate_configuration(configuration)
        configuration["__config-path"] = config
        _create_finder_charts(configuration, output, format)
    except Exception as e:
        typer.echo(str(e), err=True)
        raise typer.Exit(code=1)


def _read_configuration(config: Optional[Path]) -> Any:
    if config:
        with open(config, "r") as f:
            configuration_yaml = f.read()
    else:
        configuration_yaml = typer.get_text_stream("stdin").read()

    return yaml.safe_load(configuration_yaml.encode("utf-8"))


def _validate_configuration(configuration: Any) -> None:
    data_dir = Path(__file__).parent
    schema_file = data_dir / "schema.json"
    with open(schema_file, "r") as f:
        schema = json.load(f)

    validate(instance=configuration, schema=schema)


def _create_finder_charts(
    configuration: Dict[str, Any], output: Optional[Path], format: Optional[str]
) -> None:
    if configuration["telescope"] == "SALT":
        is_sidereal = "ra" in configuration["target"]
        if is_sidereal:
            finder_chart = _create_sidereal_salt_finder_chart(configuration)
            _save_sidereal_finder_chart(finder_chart, output, format)
        else:
            g = _create_non_sidereal_salt_finder_charts(configuration)
            _save_non_sidereal_finder_charts(g, output, format)
    else:
        raise ValueError(f"Unsupported telescope: {configuration['telescope']}")


def _create_sidereal_salt_finder_chart(configuration: Dict[str, Any]) -> FinderChart:
    # Extract the target details
    target_ = configuration["target"]
    position = SkyCoord(ra=str(target_["ra"]), dec=str(target_["dec"]))
    magnitude_range_ = target_["magnitude-range"]
    magnitude_range = MagnitudeRange(
        bandpass=magnitude_range_["bandpass"],
        min_magnitude=magnitude_range_["minimum"],
        max_magnitude=magnitude_range_["maximum"],
    )
    target = Target(
        name=target_["name"], position=position, magnitude_range=magnitude_range
    )

    # Collect the instrument-independent properties
    general = _general_properties(configuration, target)

    # Get the FITS file
    fits = _fits(configuration, target)

    # Create the finder chart
    instrument = configuration["instrument"]
    return _obtain_sidereal_salt_finder_chart(fits, general, instrument)


def _save_sidereal_finder_chart(
    finder_chart: FinderChart, output: Optional[Path], format: Optional[str]
) -> None:
    if output:
        finder_chart.save(output, format=format)
    else:
        out = io.BytesIO()
        finder_chart.save(out, format=format)
        sys.stdout.buffer.write(out.getvalue())


def _create_non_sidereal_salt_finder_charts(
    configuration: Dict[str, Any]
) -> Generator[Tuple[FinderChart, Tuple[datetime, datetime]], None, None]:
    # Extract the target details
    target_ = configuration["target"]
    name = target_["name"]
    horizons_id = target_["horizons-id"]
    ephemeris_stepsize = (
        u.Quantity(target_["ephemeris-stepsize"])
        if "ephemeris-stepsize" in target_
        else None
    )
    start_time = isoparse(target_["start-time"])
    end_time = isoparse(target_["end-time"])
    if start_time.tzinfo is None or start_time.tzinfo.utcoffset(None) is None:
        raise ValueError(
            "The start time must have a timezone offset. Examples of "
            "valid time strings are 2022-02-25T09:14:58Z and "
            "2022-11-15T16:25:34-08:00."
        )
    if end_time.tzinfo is None or end_time.tzinfo.utcoffset(None) is None:
        raise ValueError(
            "The end time must have a timezone offset. Examples of valid "
            "time strings are 2022-02-25T09:14:58Z and "
            "2022-11-15T16:25:34-08:00."
        )

    horizons_service = HorizonsService(
        object_id=horizons_id,
        location="B31",
        start=start_time,
        end=end_time,
        stepsize=ephemeris_stepsize,
    )
    ephemerides = horizons_service.ephemerides(start=start_time, end=end_time)

    def _create_finder_chart(ephemerides_: List[Ephemeris]) -> FinderChart:
        center = mid_position(ephemerides_[0].position, ephemerides_[-1].position)
        magnitude_range = ephemerides_magnitude_range(ephemerides_)
        target = Target(name=name, position=center, magnitude_range=magnitude_range)
        general = _general_properties(configuration, target)
        fits = _fits(configuration, target)
        finder_chart = _obtain_sidereal_salt_finder_chart(
            fits=fits, general=general, instrument=configuration["instrument"]
        )
        track_annotation = motion_annotation(
            ephemerides=ephemerides_, wcs=finder_chart.wcs
        )
        finder_chart.add_annotation(track_annotation)
        return finder_chart

    return FinderChart.for_time_interval(
        start=start_time,
        end=end_time,
        ephemerides=ephemerides,
        max_track_length=8 * u.arcmin,
        create_finder_chart=_create_finder_chart,
    )


def _save_non_sidereal_finder_charts(
    finder_chart_generator: Generator[
        Tuple[FinderChart, Tuple[datetime, datetime]], None, None
    ],
    output: Optional[Path],
    format: Optional[str],
) -> None:
    if output:
        out: Union[Path, BinaryIO] = output
    else:
        out = io.BytesIO()
    with zipfile.ZipFile(out, "w") as archive:
        for finder_chart, (start, end) in finder_chart_generator:
            start_time = start.astimezone(timezone.utc).strftime("%Y-%m-%d_%Hh%Mm%Ss")
            end_time = end.astimezone(timezone.utc).strftime("%Y-%m-%d_%Hh%Mm%Ss")
            fc = io.BytesIO()
            finder_chart.save(fc, format=format)
            archive.writestr(
                f"finder-chart_{start_time}_{end_time}.{format}", fc.getvalue()
            )
    if not output:
        sys.stdout.buffer.write(cast(io.BytesIO, out).getvalue())


def _general_properties(
    configuration: Dict[str, Any], target: Target
) -> GeneralProperties:
    fits_source = configuration["fits-source"]
    return GeneralProperties(
        target=target,
        position_angle=Angle(configuration["position-angle"]),
        automated_position_angle=False,
        proposal_code=configuration["proposal-code"],
        pi_family_name=configuration["pi-family-name"],
        survey=fits_source["image-survey"] if "image-survey" in fits_source else None,
    )


def _fits(configuration: Dict[str, Any], target: Target) -> Union[BinaryIO, Path]:
    fits_source = configuration["fits-source"]
    if "image-survey" in fits_source:
        survey = fits_source["image-survey"]
        fits_center = target.position
        size = 10 * u.arcmin
        return load_fits(survey, fits_center, size)
    else:
        path = Path(fits_source["file"])
        if path.is_absolute():
            return path
        elif configuration["__config-path"] is not None:
            return cast(Path, configuration["__config-path"]).parent / path
        else:
            raise ValueError(
                "The value of the fits_source.file property in the configuration must "
                "be an absolute path if the configuration is read from stdin."
            )


def _obtain_sidereal_salt_finder_chart(
    fits: Union[BinaryIO, Path], general: GeneralProperties, instrument: Dict[str, Any]
) -> FinderChart:
    if "salticam" in instrument:
        return _create_salticam_finder_chart(fits, general, instrument["salticam"])
    elif "rss" in instrument:
        return _create_rss_finder_chart(fits, general, instrument["rss"])
    elif "hrs" in instrument:
        return _create_hrs_finder_chart(fits, general)
    elif "nir" in instrument:
        return _create_nir_finder_chart(fits, general, instrument["nir"])
    else:
        raise ValueError("No supported instrument found")


def _create_salticam_finder_chart(
    fits: Union[BinaryIO, Path], general: GeneralProperties, salticam: Dict[str, Any]
) -> FinderChart:
    is_slot_mode = salticam.get("slot-mode", False)
    return sfc.salticam_finder_chart(
        fits=fits, general=general, is_slot_mode=is_slot_mode
    )


def _create_rss_finder_chart(
    fits: Union[BinaryIO, Path], general: GeneralProperties, rss: Dict[str, Any]
) -> FinderChart:
    instrument_mode = rss["mode"].lower()
    if instrument_mode == "imaging":
        is_slot_mode = rss.get("slot-mode", False)
        return sfc.rss_imaging_finder_chart(
            fits=fits, general=general, is_slot_mode=is_slot_mode
        )
    elif instrument_mode == "spectroscopy":
        slit_width = Angle(rss["slit-width"])
        slit_height = Angle(rss["slit-height"])
        return sfc.rss_longslit_finder_chart(
            fits=fits, general=general, slit_width=slit_width, slit_height=slit_height
        )
    elif instrument_mode == "mos":
        mos_mask = MosMask.from_file(rss["file"])
        reference_star_box_width = (
            Angle(rss["reference-star-box-width"])
            if "reference-star-box-width" in rss
            else Angle(5 * u.arcsec)
        )
        return sfc.rss_mos_finder_chart(
            fits=fits,
            general=general,
            mos_mask=mos_mask,
            reference_star_box_width=reference_star_box_width,
        )
    elif instrument_mode == "fabry-perot":
        return sfc.rss_fabry_perot_finder_chart(fits=fits, general=general)
    else:
        raise ValueError(f"Unsupported RSS mode: {instrument_mode}")


def _create_hrs_finder_chart(
    fits: Union[BinaryIO, Path], general: GeneralProperties
) -> FinderChart:
    return sfc.hrs_finder_chart(fits=fits, general=general)


def _create_nir_finder_chart(
    fits: Union[BinaryIO, Path], general: GeneralProperties, nir: Dict[str, Any]
) -> FinderChart:
    science_bundle_center = SkyCoord(
        ra=nir["science-bundle"]["ra"], dec=nir["science-bundle"]["dec"]
    )
    bundle_separation = Angle(nir["bundle-separation"])
    return sfc.nir_finder_chart(
        fits=fits,
        general=general,
        science_bundle_center=science_bundle_center,
        bundle_separation=bundle_separation,
    )


if __name__ == "__main__":
    app()
