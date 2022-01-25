from imephu.annotation import salt
from imephu.finder_chart import FinderChart


def test_title(fits_file, check_finder):
    """Test the title annotation."""
    finder_chart = FinderChart(fits_file)
    title_annotation = salt.title(
        target="Magrathea",
        proposal_code="2022-1-SCI-042",
        pi_family_name="Adams",
        wcs=finder_chart.wcs,
    )
    finder_chart.add_annotation(title_annotation)
    check_finder(finder_chart)
