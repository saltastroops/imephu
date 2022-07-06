import io
import time
import zipfile
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest import mock

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from typer.testing import CliRunner

import imephu
from imephu.cli import app
from imephu.utils import Ephemeris

runner = CliRunner()


def test_version():
    """Test that the version is output for the --version option."""
    result = runner.invoke(app, ["--version"])
    assert "imephu" in result.stdout
    assert imephu.__version__ in result.stdout


def test_read_configuration_from_file(tmp_path):
    """Test that a finder chart configuration can be read from file."""
    config_file = tmp_path / "config.yaml"
    config_file.write_text("xyz123: configuration")

    result = runner.invoke(app, ["--config", str(config_file)])
    assert "Failed validating" in result.stdout and "xyz123" in result.stdout


def test_read_configuration_from_stdin():
    """Test that a finder chart configuration can be read from stdin."""
    result = runner.invoke(app, [], input="xyz123: configuration")
    assert "Failed validating" in result.stdout and "xyz123" in result.stdout


def test_read_configuration_with_absolute_fits_file_path(tmp_path, check_cli):
    """Test that the configuration may reference an absolute FITS file path."""
    fits_file = Path(__file__).parent / "data" / "ra10_dec-60.fits"
    fits_source_yaml = f"""\
fits-source:
  file: {fits_file}
    """
    instrument_yaml = """\
instrument:
  salticam:
    slot-mode: false
"""
    check_cli(instrument_yaml, fits_source_yaml)


def test_read_configuration_with_relative_fits_file_path(tmp_path, check_image):
    """Test that the configuration may reference a relative FITS file path."""
    fits_file = Path(__file__).parent / "data" / "ra10_dec-60.fits"
    fits_data = fits_file.read_bytes()
    referenced_fits_file = tmp_path / "magrathea.fits"
    referenced_fits_file.write_bytes(fits_data)

    configuration = """
telescope: SALT
pi-family-name: Doe
proposal-code: 2022-1-SCI-042
position-angle: 30d
fits-source:
  file: magrathea.fits
target:
  name: Magrathea
  ra: 0h 40m 00s
  dec: -60d
  magnitude-range:
    bandpass: V
    minimum: 17
    maximum: 17.3
instrument:
  salticam: {}
    """
    config = tmp_path / "config.yaml"
    config.write_text(configuration)
    output = tmp_path / "finder_chart.png"

    np.random.seed(0)
    try:
        runner.invoke(app, ["--config", config, "--out", output])
        finder_chart = output.read_bytes()
        check_image(io.BytesIO(finder_chart))
    finally:
        np.random.seed()


def test_read_configuration_from_stdin_requires_absolute_fits_path(tmp_path):
    """Test that the FITS file path must be absolute when reading config from stdin."""
    configuration = """
telescope: SALT
pi-family-name: Doe
proposal-code: 2022-1-SCI-042
position-angle: 30d
fits-source:
  file: relative-path.fits
target:
  name: Magrathea
  ra: 0h 40m 00s
  dec: -60d
  magnitude-range:
    bandpass: V
    minimum: 17
    maximum: 17.3
instrument:
  salticam: {}
    """
    result = runner.invoke(app, input=configuration)
    assert result.exit_code != 0
    assert "read from stdin" in result.stdout


def test_configuration_must_be_valid():
    """Test that invalid finder chart configurations are rejected."""
    result = runner.invoke(app, [], input="invalid: true")
    assert result.exit_code != 0
    assert "Failed validating" in result.stdout


def test_start_time_must_have_a_timezone_offset(tmp_path):
    """Test that a start time string must include a timezone offset."""
    configuration = """
telescope: SALT
pi-family-name: Doe
proposal-code: 2022-1-SCI-042
position-angle: 30d
fits-source:
  image-survey: POSS2/UKSTU Red
target:
  name: Magrathea
  horizons-id: margathea
  start-time: "2022-02-05T17:00:00"
  end-time: "2022-02-04T17:00:00+02:00"
  ephemeris-stepsize: 1h
instrument:
  salticam: {}
    """
    config = tmp_path / "config.yaml"
    config.write_text(configuration)
    result = runner.invoke(app, ["--config", config])
    assert "start time" in result.stdout and "timezone" in result.stdout


def test_end_time_must_have_a_timezone_offset(tmp_path):
    """Test that an end time string must include a timezone offset."""
    configuration = """
telescope: SALT
pi-family-name: Doe
proposal-code: 2022-1-SCI-042
position-angle: 30d
fits-source:
  image-survey: POSS2/UKSTU Red
target:
  name: Magrathea
  horizons-id: margathea
  start-time: "2022-02-05T17:00:00Z"
  end-time: "2022-02-04T17:00:00"
  ephemeris-stepsize: 1h
instrument:
  salticam: {}
    """
    config = tmp_path / "config.yaml"
    config.write_text(configuration)
    result = runner.invoke(app, ["--config", config])
    assert "end time" in result.stdout and "timezone" in result.stdout


@pytest.mark.parametrize("slot_mode", [None, True, False])
def test_create_salticam_finder_chart(slot_mode, check_cli, mock_from_survey):
    """Test creating a Salticam finder chart with the CLI."""
    if slot_mode is None:
        instrument_yaml = """\
instrument:
  salticam: {}
"""
    elif slot_mode is True:
        instrument_yaml = """\
instrument:
  salticam:
    slot-mode: true
"""
    else:
        instrument_yaml = """\
instrument:
  salticam:
    slot-mode: false
"""
    check_cli(instrument_yaml)


@pytest.mark.parametrize("slot_mode", [None, True, False])
def test_create_rss_imaging_finder_chart(slot_mode, check_cli, mock_from_survey):
    """Test creating a Salticam finder chart with the CLI."""
    if slot_mode is None:
        slot_mode_yaml = ""
    elif slot_mode is True:
        slot_mode_yaml = "slot-mode: true"
    else:
        slot_mode_yaml = "slot-mode: false"
    instrument_yaml = f"""\
instrument:
  rss:
    mode: imaging
    {slot_mode_yaml}
"""
    check_cli(instrument_yaml)


def test_create_rss_spectroscopy_finder_chart(check_cli, mock_from_survey):
    """Test creating an RSS spectroscopy finder chart with the CLI."""
    instrument_yaml = """\
instrument:
  rss:
    mode: spectroscopy
    slit-width: 4 arcsec
    slit-height: 8 arcmin
"""
    check_cli(instrument_yaml)


@pytest.mark.parametrize(
    "file,reference_star_box_width",
    [
        (Path(__file__).parent / "data" / "mos_mask2.xml", "10 arcsec"),
        (Path(__file__).parent / "data" / "mos_mask2.rsmt", None),
    ],
)
def test_create_rss_mos_finder_chart(
    file, reference_star_box_width, check_cli, mock_from_survey
):
    """Test creating an RSS MOS finder chart with the CLI."""
    if reference_star_box_width is not None:
        box_width_yaml = f"reference-star-box-width: {reference_star_box_width}"
    else:
        box_width_yaml = ""
    instrument_yaml = f"""\
instrument:
  rss:
    mode: MOS
    file: {str(file)}
    {box_width_yaml}
"""
    check_cli(instrument_yaml)


def test_create_rss_fabry_perot_finder_chart(check_cli, mock_from_survey):
    """Test creating an RSS Fabry-Perot finder chart with the CLI."""
    instrument_yaml = """\
instrument:
  rss:
    mode: Fabry-Perot
"""
    check_cli(instrument_yaml)


def test_create_hrs_finder_chart(check_cli, mock_from_survey):
    """Test creating an HRS finder chart with the CLI."""
    instrument_yaml = """\
instrument:
  hrs: {}
"""
    check_cli(instrument_yaml)


def test_create_nir_finder_chart(check_cli, mock_from_survey):
    """Test creating an NIR finder chart with the CLI."""
    instrument_yaml = """
instrument:
  nir:
    bundle-separation: 50 arcsec
    science-bundle:
      ra: 0h 40m 0s
      dec: -60d
"""
    check_cli(instrument_yaml)


@pytest.mark.parametrize(
    "finder_chart_file",
    [
        "finder-chart_2022-02-17_00h00m00s_2022-02-17_02h00m00s.png",
        "finder-chart_2022-02-17_02h00m00s_2022-02-17_04h00m00s.png",
    ],
)
def test_create_non_sidereal_salt_finder_charts(
    finder_chart_file, fits_file, fits_file2, tmp_path_factory, check_image
):
    """Test creating non-sidereal SALT finder charts."""
    t = datetime(2022, 2, 17, 0, 0, 0, 0, tzinfo=timezone.utc)
    hour = timedelta(hours=1)
    start = t + 0.5 * hour
    end = t + 3.5 * hour
    start_time = start.astimezone(timezone.utc).isoformat()
    end_time = end.astimezone(timezone.utc).isoformat()
    configuration = f"""\
telescope: SALT
pi-family-name: Doe
proposal-code: 2022-1-SCI-042
position-angle: 30d
fits-source:
  image-survey: POSS2/UKSTU Red
target:
  name: Magrathea
  horizons-id: magrathea
  start-time: "{start_time}"
  end-time: "{end_time}"
  ephemeris-stepsize: 1h
instrument:
  salticam:
    slot-mode: false
"""

    ephemerides = [
        Ephemeris(
            epoch=t,
            position=SkyCoord(ra="0h40m30s", dec=-60 * u.deg),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + hour,
            position=SkyCoord(ra="0h40m00s", dec=-60 * u.deg),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + 2 * hour,
            position=SkyCoord(ra="0h39m30s", dec=-60 * u.deg),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + 3 * hour,
            position=SkyCoord(ra="0h39m00s", dec=-60 * u.deg),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + 4 * hour,
            position=SkyCoord(ra="0h38m30s", dec=-60 * u.deg),
            magnitude_range=None,
        ),
    ]
    with mock.patch.object(
        imephu.cli, "HorizonsService", autospec=True
    ) as mock_horizons:

        mock_horizons.return_value.ephemerides.return_value = ephemerides
        np.random.seed(0)

        with mock.patch.object(
            imephu.cli, "load_fits", autospec=True
        ) as mock_load_fits:
            mock_load_fits.side_effect = [
                io.BytesIO(fits_file.read_bytes()),
                io.BytesIO(fits_file2.read_bytes()),
            ]
            try:
                tmp = tmp_path_factory.mktemp(f"finder-chart-{time.time_ns()}")
                config = tmp / "config.yaml"
                config.write_text(configuration)
                result = runner.invoke(app, ["--config", config])
                zip_content = io.BytesIO(result.stdout_bytes)
                with zipfile.ZipFile(zip_content) as archive:
                    assert len(archive.filelist) == 2
                    finder_chart = archive.read(finder_chart_file)
                    check_image(io.BytesIO(finder_chart))
            finally:
                np.random.seed()


def test_use_format_option_for_sidereal_finder_chart(fits_file, tmp_path_factory):
    """Test using the --format option when creating a non-sidereal finder chart."""
    configuration = """\
telescope: SALT
proposal-code: 2022-1-SCI-042
pi-family-name: Doe
position-angle: 0d
target:
  name: LMC
  ra:  05h 23m 34.5s
  dec: -69d 45m 22s
  magnitude-range:
    bandpass: V
    minimum: 0.9
    maximum: 0.9
fits-source:
  image-survey: POSS2/UKSTU Red
instrument:
  salticam:
    slot-mode: false
"""
    np.random.seed(0)
    try:
        tmp = tmp_path_factory.mktemp(f"finder-chart-{time.time_ns()}")
        config = tmp / "config.yaml"
        config.write_text(configuration)
        output = tmp / "finder_chart.png"
        with mock.patch.object(
            imephu.cli, "load_fits", autospec=True
        ) as mock_load_fits:
            fits = fits_file.read_bytes()
            mock_load_fits.return_value = io.BytesIO(fits)
            runner.invoke(app, ["--config", config, "--out", output, "--format", "pdf"])
            finder_chart = output.read_bytes()
            assert finder_chart.startswith(b"%PDF")
    finally:
        np.random.seed()


def test_use_format_option_for_non_sidereal_finder_chart(fits_file, tmp_path_factory):
    """Test using the --format option when creating a non-sidereal finder chart."""
    t = datetime(2022, 2, 17, 0, 0, 0, 0, tzinfo=timezone.utc)
    hour = timedelta(hours=1)
    start = t + 0.5 * hour
    end = t + 0.6 * hour
    start_time = start.astimezone(timezone.utc).isoformat()
    end_time = end.astimezone(timezone.utc).isoformat()
    configuration = f"""\
    telescope: SALT
    pi-family-name: Doe
    proposal-code: 2022-1-SCI-042
    position-angle: 30d
    fits-source:
      image-survey: POSS2/UKSTU Red
    target:
      name: Magrathea
      horizons-id: magrathea
      start-time: "{start_time}"
      end-time: "{end_time}"
      ephemeris-stepsize: 1h
    instrument:
      salticam:
        slot-mode: false
    """

    ephemerides = [
        Ephemeris(
            epoch=t,
            position=SkyCoord(ra="0h40m30s", dec=-60 * u.deg),
            magnitude_range=None,
        ),
        Ephemeris(
            epoch=t + hour,
            position=SkyCoord(ra="0h40m00s", dec=-60 * u.deg),
            magnitude_range=None,
        ),
    ]
    with mock.patch.object(
        imephu.cli, "HorizonsService", autospec=True
    ) as mock_horizons:

        mock_horizons.return_value.ephemerides.return_value = ephemerides
        np.random.seed(0)

        with mock.patch.object(
            imephu.cli, "load_fits", autospec=True
        ) as mock_load_fits:

            mock_load_fits.return_value = fits_file
            tmp = tmp_path_factory.mktemp(f"finder-chart-{time.time_ns()}")
            config = tmp / "config.yaml"
            config.write_text(configuration)
            result = runner.invoke(app, ["--config", config, "--format", "pdf"])
            zip_content = io.BytesIO(result.stdout_bytes)
            with zipfile.ZipFile(zip_content) as archive:
                finder_chart = archive.read(
                    "finder-chart_2022-02-17_00h00m00s_2022-02-17_01h00m00s.pdf"
                )
                assert finder_chart[:10].startswith(b"%PDF")
