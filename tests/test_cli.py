from pathlib import Path

import pytest
from typer.testing import CliRunner

import imephu
from imephu.cli import app

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


def test_configuration_must_be_valid():
    """Test that invalid finder chart configurations are rejected."""
    result = runner.invoke(app, [], input="invalid: true")
    assert result.exit_code != 0
    assert "Failed validating" in result.stdout


@pytest.mark.parametrize("slot_mode", [None, True, False])
def test_create_salticam_finder_chart(slot_mode, check_cli, mock_salt_load_fits):
    """Test creating a Salticam finder chart with the CLI."""
    if slot_mode is None:
        slot_mode_yaml = ""
    elif slot_mode is True:
        slot_mode_yaml = "slot-mode: true"
    else:
        slot_mode_yaml = "slot-mode: false"
    instrument_yaml = f"""\
instrument:
  name: Salticam
  {slot_mode_yaml}
"""
    check_cli(instrument_yaml)


@pytest.mark.parametrize("slot_mode", [None, True, False])
def test_create_rss_imaging_finder_chart(slot_mode, check_cli, mock_salt_load_fits):
    """Test creating a Salticam finder chart with the CLI."""
    if slot_mode is None:
        slot_mode_yaml = ""
    elif slot_mode is True:
        slot_mode_yaml = "slot-mode: true"
    else:
        slot_mode_yaml = "slot-mode: false"
    instrument_yaml = f"""\
instrument:
  name: RSS
  mode: imaging
  {slot_mode_yaml}
"""
    check_cli(instrument_yaml)


def test_create_rss_spectroscopy_finder_chart(check_cli, mock_salt_load_fits):
    """Test creating an RSS spectroscopy finder chart with the CLI."""
    instrument_yaml = """\
instrument:
  name: RSS
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
    file, reference_star_box_width, check_cli, mock_salt_load_fits
):
    """Test creating an RSS MOS finder chart with the CLI."""
    if reference_star_box_width is not None:
        box_width_yaml = f"reference-star-box-width: {reference_star_box_width}"
    else:
        box_width_yaml = ""
    instrument_yaml = f"""\
instrument:
  name: RSS
  mode: MOS
  file: {str(file)}
  {box_width_yaml}
"""
    check_cli(instrument_yaml)


def test_create_rss_fabry_perot_finder_chart(check_cli, mock_salt_load_fits):
    """Test creating an RSS Fabry-Perot finder chart with the CLI."""
    instrument_yaml = """\
instrument:
  name: RSS
  mode: Fabry-Perot
"""
    check_cli(instrument_yaml)


def test_create_hrs_finder_chart(check_cli, mock_salt_load_fits):
    """Test creating an HRS finder chart with the CLI."""
    instrument_yaml = """\
instrument:
  name: HRS
"""
    check_cli(instrument_yaml)


def test_create_nir_finder_chart(check_cli, mock_salt_load_fits):
    """Test creating an NIR finder chart with the CLI."""
    instrument_yaml = """
instrument:
  name: NIR
  bundle-separation: 50 arcsec
  science-bundle:
    ra: 0h 40m 0s
    dec: -60d
"""
    check_cli(instrument_yaml)
