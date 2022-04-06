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
alt: Finder chart for a Salticam setup
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

### Finder charts for non-sidereal targets

When creating finder charts for non-sidereal targets, you have to specify the asteroid's identifier for the [Horizons database](https://ssd.jpl.nasa.gov/horizons/app.html#/), a start time, an end time and a stepsize between ephemerides instead of the right ascension, declination and magnitude range. The start and end time must be given as an ISO-8601 string with timezone offset.

Here is an example configuration for observing the asteroid Ubuntu with Salticam.

```{literalinclude} configuration-examples/salticam-non-sidereal.yaml
---
language: yaml
---
```

```{warning}
As shown in the example, the values for the start and end date *must* be enclosed in double quotes. Otherwise you will get an error that the value is not of type 'string'.
```

Save the configuration as `salticam-non-sidereal.yaml` and run imephu with this file.

```shell
imephu --config salticam-non-sidereal.yaml --out ubuntu.zip
```

If the path of the target does not fit onto a single finder chart, multiple finder charts are generated, and a zip file of all the charts is created. The name of the finder chart file within this zip file include the start and end of the time interval for which the chart is valid.

Whereas for sidereal finder charts the image type can be inferred from the extension of the chosen output file, this is not the case for non-sidereal targets. However, you can request a specific type (such as ``pdf`` or ``png``) by using the ``--format`` option. For example:

```shell
imephu --config salticam-non-sidereal.yaml --out ubuntu.zip --format png
```

Matplotlib's default format is used if no format is specified. This is usually PNG but can be changed by [setting rcParams("savefig.format")](https://matplotlib.org/3.5.0/tutorials/introductory/customizing.html?highlight=savefig.format#a-sample-matplotlibrc-file).

## Writing the finder chart to stdout

So far we have always used the `--out` option to save the generated file in a file. However, if you want to use imephu as part of a Unix pipeline, you can omit the option so that the file content is written to stdout instead. For example:

```shell
imephu --config salticam.yaml | less
```

Matplotlib's default format is used. This is usually PNG but can be changed by [setting rcParams("savefig.format")](https://matplotlib.org/3.5.0/tutorials/introductory/customizing.html?highlight=savefig.format#a-sample-matplotlibrc-file). You can change this with the `--format` option.

```shell
imephu --config salticam.yaml --format pdf | less
```

You might wonder what happens if you use both the `--out` and the `--format` option. In this case, `--format` sets the image format, and the created file is named as defined by `--out` is used, irrespective of whether the file extension is consistent with the chosen image format. For example, you hence might end up with a pdf file called `finder-chart.png`.

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

