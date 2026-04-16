from astropy.coordinates import Angle, SkyCoord

from imephu.annotation.general.arrow import ArrowAnnotation
from imephu.finder_chart import FinderChart


def test_arrow_annotation(fits_file, check_finder):
    """Test arrow annotations."""
    finder_chart = FinderChart(fits_file)
    arrow1_annotation = ArrowAnnotation(
        start=SkyCoord(ra="00h40m30s", dec="-60d04m00"),
        end=SkyCoord(ra="00h40m30s", dec="-59d58m00s"),
        wcs=finder_chart.wcs,
        head_width=Angle("2.5arcmin"),  # noqa
        head_height=Angle("2arcmin"),  # noqa
        head_filled=True,
        color="green",
        linewidth=10,
    )
    arrow2_annotation = ArrowAnnotation(
        start=SkyCoord(ra="00h40m20s", dec="-60d02m00s"),
        end=SkyCoord(ra="00h39m40s", dec="-59d58m00s"),
        wcs=finder_chart.wcs,
        head_width=Angle("60arcsec"),  # noqa
        head_height=Angle("30arcsec"),  # noqa
        color="orange",
        linewidth=4,
    )
    arrow3_annotation = ArrowAnnotation(
        start=SkyCoord(ra="00h40m20s", dec="-60d04m00s"),
        end=SkyCoord(ra="00h39m40s", dec="-60d00m00s"),
        wcs=finder_chart.wcs,
        head_width=Angle("60arcsec"),  # noqa
        head_height=Angle("30arcsec"),  # noqa
        tail_shown=False,
        color="blue",
        linewidth=4,
    )
    arrow4_annotation = ArrowAnnotation(
        start=SkyCoord(ra="00h40m20s", dec="-60d06m00s"),
        end=SkyCoord(ra="00h39m40s", dec="-60d02m00s"),
        wcs=finder_chart.wcs,
        head_width=Angle("20arcsec"),  # noqa
        head_height=Angle("30arcsec"),  # noqa
        head_filled=True,
        tail_shown=False,
        linewidth=4,
    )

    finder_chart.add_annotation(arrow1_annotation)
    finder_chart.add_annotation(arrow2_annotation)
    finder_chart.add_annotation(arrow3_annotation)
    finder_chart.add_annotation(arrow4_annotation)

    check_finder(finder_chart)
