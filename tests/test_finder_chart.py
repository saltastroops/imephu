import io
import pathlib
from unittest import mock

import numpy as np
import pikepdf
import pytest
from astropy import units as u

import imephu.service.survey
from imephu.finder_chart import FinderChart


def test_finder_chart_is_generated_from_path(check_finder):
    """Test that finder charts can be generated with a `~pathlib.Path object."""
    path = pathlib.Path(__file__).parent / "data" / "ra10_dec-60.fits"
    finder_chart = FinderChart(path)
    check_finder(finder_chart)


def test_finder_chart_is_generated_from_string(check_finder):
    """Test that finder charts can be generated with a string containing the path."""
    name = str(pathlib.Path(__file__).parent / "data" / "ra10_dec-60.fits")
    finder_chart = FinderChart(name)
    check_finder(finder_chart)


def test_finder_chart_is_generated_from_stream(check_finder):
    """Test that finder charts can be generated with a binary stream."""
    path = pathlib.Path(__file__).parent / "data" / "ra10_dec-60.fits"
    with open(path, "rb") as f:
        finder_chart = FinderChart(f)
    check_finder(finder_chart)


# Formats other than jpg or png may produce different files for different runs,
# so that they cannot be tested with pytest-regressions.
@pytest.mark.parametrize("format", ["jpg", "png"])
def test_finder_chart_export_formats(format, file_regression, fits_file):
    """Test that finder charts can be exported to different file formats."""
    np.random.seed(0)
    finder_chart = FinderChart(fits_file)
    contents = io.BytesIO()
    finder_chart.save(contents, format=format)
    file_regression.check(contents.getvalue(), binary=True, extension=f".{format}")


def test_finder_chart_from_survey_returns_finder_chart(
    fits_file, fits_center, check_finder
):
    """Test that the from_survey method returns a finder chart."""
    with mock.patch.object(
        imephu.finder_chart, "load_fits", autospec=True
    ) as mock_load_fits:
        mock_load_fits.return_value = open(fits_file, "rb")
        finder_chart = FinderChart.from_survey(
            "POSS2/UKSTU Red", fits_center, 10 * u.arcmin
        )
        check_finder(finder_chart)


def test_metadata_is_added_to_finder_chart_pdf(fits_file):
    """Test that metadata is added to finder chart pdf files."""
    # save the pdf...
    finder_chart = FinderChart(fits_file)
    pdf = io.BytesIO()
    finder_chart.save(name=pdf, format="pdf")

    # ... and check that the metadata has been added
    document = pikepdf.open(io.BytesIO(pdf.getvalue()))
    meta = document.open_metadata()
    assert meta["dc:title"] == "Finder Chart"
    assert meta["xmp:CreatorTool"] == f"imephu {imephu.__version__}"
