import io
import json
import sys
from pathlib import Path
from typing import Any, BinaryIO, Dict, Optional, Union

import typer
import yaml
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from jsonschema import validate

import imephu
import imephu.salt.finder_chart as sfc
from imephu.finder_chart import FinderChart
from imephu.salt.finder_chart import GeneralProperties, Target
from imephu.service.survey import load_fits
from imephu.utils import MagnitudeRange

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
    version: Optional[bool] = typer.Option(
        None, "--version", callback=_version_callback, help="Show the version and exit."
    ),
) -> None:
    """A tool for creating finder charts."""  # noqa: D401
    configuration = _read_configuration(config)
    try:
        _validate_configuration(configuration)
        _create_finder_charts(configuration, output)
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
    configuration: Dict[str, Any], output: Optional[Path]
) -> None:
    if configuration["telescope"] == "SALT":
        is_sidereal = "ra" in configuration["target"]
        if is_sidereal:
            finder_chart = _create_sidereal_salt_finder_chart(configuration)
            if output:
                finder_chart.save(output)
            else:
                out = io.BytesIO()
                finder_chart.save(out, "pdf")
                sys.stdout.buffer.write(out.getvalue())
        else:
            raise NotImplementedError("Non-sidereal targets are not supported yet")
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
    fits_source = configuration["fits-source"]
    general = GeneralProperties(
        target=target,
        position_angle=Angle(configuration["position-angle"]),
        automated_position_angle=False,
        proposal_code=configuration["proposal-code"],
        pi_family_name=configuration["pi-family-name"],
        survey=fits_source["image-survey"] if "image-survey" in fits_source else None,
    )

    # Get the FITS file
    if "image-survey" in fits_source:
        survey = fits_source["image-survey"]
        fits_center = position
        size = 10 * u.arcmin
        fits: Union[BinaryIO, Path] = load_fits(survey, fits_center, size)
    else:
        fits = Path(fits_source["file"])

    # Create the finder chart
    instrument = configuration["instrument"]
    instrument_name = instrument["name"].lower()
    if instrument_name == "salticam":
        return _create_salticam_finder_chart(fits, general, instrument)
    elif instrument_name == "rss":
        return _create_rss_finder_chart(fits, general, instrument)
    else:
        raise ValueError(f"Unsupported instrument: {instrument_name}")


def _create_salticam_finder_chart(fits: Union[BinaryIO, Path], general: GeneralProperties, instrument: Dict[str, Any]) -> FinderChart:
    is_slot_mode = instrument.get("slot-mode", False)
    return sfc.salticam_finder_chart(
        fits=fits, general=general, is_slot_mode=is_slot_mode
    )


def _create_rss_finder_chart(fits: Union[BinaryIO, Path], general: GeneralProperties, instrument: Dict[str, Any]) -> FinderChart:
    instrument_mode = instrument["mode"].lower()
    if instrument_mode == "imaging":
        is_slot_mode = instrument.get("slot-mode", False)
        return sfc.rss_imaging_finder_chart(
            fits=fits, general=general, is_slot_mode=is_slot_mode
        )
    elif instrument_mode == "spectroscopy":
        slit_width = Angle(instrument["slit-width"])
        slit_height = Angle(instrument["slit-height"])
        return sfc.rss_longslit_finder_chart(fits=fits, general=general, slit_width=slit_width, slit_height=slit_height)
    else:
        raise ValueError(f"Unsupported RSS mode: {instrument_mode}")


if __name__ == "__main__":
    app()
