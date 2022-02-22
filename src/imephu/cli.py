import json
from pathlib import Path
from typing import Any, Dict, Optional

import typer
import yaml
from jsonschema import validate

import imephu

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
    version: Optional[bool] = typer.Option(
        None, "--version", callback=_version_callback, help="Show the version and exit."
    ),
) -> None:
    """A tool for creating finder charts."""  # noqa: D401
    configuration = _read_configuration(config)
    try:
        _validate_configuration(configuration)
    except Exception as e:
        typer.echo(str(e), err=True)
        raise typer.Exit(code=1)
    _create_finder_charts(configuration)


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


def _create_finder_charts(configuration: Dict[str, Any]) -> None:
    pass


if __name__ == "__main__":
    app()
