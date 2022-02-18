# Quickstart

*imephu* is the word for *map* in isiXhosa, one of the 11 official languages of South Africa. The package allows you to easily create finder charts, primarily for the [Southern African Large Telescope (SALT)](https://www.salt.ac.za).

## Installation

You can install imephu with pip:

```shell
pip install imephu
```

## Creating ready-made finder charts

We start by creating a barebones finder chart.

```python
from astropy import units as u
from astropy.coordinates import SkyCoord
from imephu.finder_chart import FinderChart

survey = "POSS2/UKSTU Red"
fits_center = SkyCoord(ra="20h17m03s", dec="-35d08m14s")
size = 10 * u.arcmin
finder_chart = FinderChart.from_survey(survey, fits_center, size)
finder_chart.show()
```

As you can see, you need to supply the name of the image survey, the center position and the size. A FITS image is requested from the survey for the give position and size, and a coordinate grid is overlaid. The result is a completely telescope agnostic finder chart.

Of course this finder chart is of limited use. It would be much more practical if the target position (i.e. the center of the image) was highlighted along with whatever instrument fiducials would be of value. We could add al of these manually (and we'll explain how to do that later), but for SALT observations there is an easier option: Use the utility function for the instrument (and mode) under consideration.

For example, you can create a finder chart for an RSS longslit observation as follows.

```python
from astropy import units as u
from astropy.coordinates import SkyCoord
from imephu.salt.finder_chart import (
    rss_longslit_finder_chart,
    GeneralProperties,
    Target,
)
from imephu.utils import MagnitudeRange

survey = "POSS2/UKSTU Red"
fits_center = SkyCoord(ra="20h17m03s", dec="-35d08m14s")
magnitude_range = MagnitudeRange(bandpass="V", min_magnitude=17, max_magnitude=17.3)
target = Target(
    name="My Random Target", position=fits_center, magnitude_range=magnitude_range
)
general = GeneralProperties(
    target=target,
    position_angle=30 * u.deg,
    automated_position_angle=False,
    proposal_code="2022-1-SCI-123",
    pi_family_name="Doe",
    survey=survey,
)
finder_chart = rss_longslit_finder_chart(
    general=general, slit_width=2 * u.arcsec, slit_height=8 * u.arcmin
)
finder_chart.show()
```

