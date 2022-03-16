# Quickstart

*imephu* is the word for *map* in isiXhosa, one of the 11 official languages of South Africa. The package allows you to easily create finder charts, primarily for the [Southern African Large Telescope (SALT)](https://www.salt.ac.za).

## Installation

You can install imephu with pip:

```shell
pip install imephu
```

## Creating off-the-shelf finder charts

imephu ships with a command line interface with command ``imephu``. Let us run it to find out about its options:

```shell
imephu --help
```

We'll get some helpful information:

```text
Usage: imephu [OPTIONS]

  A tool for creating finder charts.

Options:
  -c, --config PATH     Configuration details for the finder chart(s).
                        [required]
  -o, --out PATH        Output file.
  --version             Show the version and exit.
  --install-completion  Install completion for the current shell.
  --show-completion     Show completion for the current shell, to copy it or
                        customize the installation.
  --help                Show this message and exit.
```

We can see that we need at least a configuration file. This must be a YAML file containing all the details necessary for creating the finder chart, such as the telescope and the instrument configuration. Let's start with a SALT setup using the Salticam imaging camera.

```{literalinclude} configuration-examples/salticam.yaml
---
language: yaml
---
```

Save the file as `salticam.yaml` and then use it to create a finder chart:

```shell
imephu --config salticam.yaml --out salticam.pdf
```

We get the following finder chart.

```{image} img/finder-charts/lmc-salticam.png
---
alt: Finder chart for a Salticam setuo
align: center
---
```

Let's also try to create an RSS finder chart for the following configuration.

```{literalinclude} configuration-examples/rss-longslit.yaml
---
language: yaml
---
```

Save the configuration as ``rss-longslit.yaml`` and run the command for creating a finder chart from it:

```shell
imephu --config rss-longslit.yaml --out rss-longslit.png
```

While previously we created a pdf, we now get a png, which looks like this:

```{image} img/finder-charts/ngc6000-rss-longslit.png
---
alt: Finder chart for an RSS longslit setup
align: center
---
```

There may be cases where the standard FITS image obtained from one of the available image servers is not the best one to use. For example, you might want to use a zoomed in version. Let's try this with our RSS finder chart! Here is a setup using a FITS image from file.

```{literalinclude} configuration-examples/rss-longslit-own-fits.yaml
---
language: yaml
---
```

The only difference to the previous configuration is that the `fits-source` now specifies a file rather than an image survey. The file path must be absolute or relative to the configuration file (not the present working directory).

Save the configuration file as ```rss-longslit-own-fits.yaml``` and make sure that you have a FITS file ``ngc6000-zoomed-in.fits`` in the same directory.

```{note}
You can download the FITS file {download}`here <_static/ngc6000-zoomed-in.fits>`.
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

