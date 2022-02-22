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
