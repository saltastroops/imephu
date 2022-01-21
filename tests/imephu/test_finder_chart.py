import io
import pathlib

import numpy as np
import pytest

from imephu.finder_chart import FinderChart


def test_finder_chart_is_generated_from_path(check_finder):
    path = pathlib.Path(__file__).parent.parent / "data" / "ra10_dec-60.fits"
    finder_chart = FinderChart(path)
    check_finder(finder_chart)


def test_finder_chart_is_generated_from_string(check_finder):
    name = str(pathlib.Path(__file__).parent.parent / "data" / "ra10_dec-60.fits")
    finder_chart = FinderChart(name)
    check_finder(finder_chart)


def test_finder_chart_is_generated_from_stream(check_finder):
    path = pathlib.Path(__file__).parent.parent / "data" / "ra10_dec-60.fits"
    with open(path, "rb") as f:
        finder_chart = FinderChart(f)
    check_finder(finder_chart)


# Formats other than jpg or png may produce different files for different runs,
# so that they cannot be tested with pytest-regressions.
@pytest.mark.parametrize("format", ["jpg", "png"])
def test_finder_chart_export_formats(format, file_regression, fits_file):
    np.random.seed(0)
    finder_chart = FinderChart(fits_file)
    contents = io.BytesIO()
    finder_chart.save(contents, format=format)
    file_regression.check(contents.getvalue(), binary=True, extension=f".{format}")
