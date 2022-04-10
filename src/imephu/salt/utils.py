import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import List

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from defusedxml.minidom import parseString


@dataclass()
class MosMaskSlit:
    """A slit in a MOS mask.

    The slit is characterised by its center position, width, height and tilt. The tilt
    should not be confused with the position angle of the MOS mask. The slit will be
    rotated by the sum of the position angle and its tilt.

    Attributes
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The position of the slit center, as right ascension and declination.
    width: `~astropy.coordinates.Angle`
        The slit width, as an angle on the sky.
    height: `~astropy.coordinates.Angle`
        The slit height, as an angle on the sky.
    tilt: `~astropy.coordinates.Angle`
        The slit tilt, as an angle on the sky measured from north to east.
    """

    center: SkyCoord
    width: Angle
    height: Angle
    tilt: Angle


class MosMask:
    """A SALT MOS mask.

    Parameters
    ----------
    mask_xml: str
        XML defining the MOS mask content.

    Attributes
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The position of the mask center, as right ascension and declination.
    position_angle: `~astropy.coordinates.Angle`
        The position angle of the mask.
    reference_stars: list of `astropy.coordinates.SkyCoord`
        The positions of the reference stars, as right ascensions and declinations.
    slits: list of `MosMaskSlit`
        The slits on the mask.
    """

    _center: SkyCoord

    _position_angle: Angle

    _reference_stars: List[SkyCoord]

    _slits: List[MosMaskSlit]

    def __init__(self, mask_xml: str) -> None:
        self._init_from_xml(mask_xml)

    @staticmethod
    def from_file(name: Path) -> "MosMask":
        """Return the `~imephu.salt.utils.MosMask defined in a file.

        The file can be an RSMT file, as produced by SALT's RSS Slit Mask Tool, or an
        XML file.

        Parameters
        ----------
        path: `~pathlib.Path`
            Path of the file defining the MOS mask.
        """
        if zipfile.is_zipfile(name):
            with zipfile.ZipFile(name, "r") as archive:
                mask_xml = archive.read("Slitmask.xml").decode("utf-8")
        else:
            with open(name, "r") as f:
                mask_xml = f.read()

        return MosMask(mask_xml)

    @property
    def center(self) -> SkyCoord:
        """Return the mask center."""
        return self._center

    @property
    def position_angle(self) -> Angle:
        """Return the position angle."""
        return self._position_angle

    @property
    def reference_stars(self) -> List[SkyCoord]:
        """Return the list of reference star positions."""
        return self._reference_stars

    @property
    def slits(self) -> List[MosMaskSlit]:
        """Return the list of slits."""
        return self._slits

    def _init_from_xml(self, mask_xml: str) -> None:
        # Get the DOM
        doc = parseString(mask_xml)

        # Extract the mask position and rotation angle
        parameter_elements = doc.getElementsByTagName("parameter")
        parameters = {}
        for par in parameter_elements:
            name = par.getAttribute("name")
            val = par.getAttribute("value")
            parameters[name] = val
        ra = float(parameters["CENTERRA"]) * u.deg
        dec = float(parameters["CENTERDEC"]) * u.deg
        self._center = SkyCoord(ra=ra, dec=dec)
        self._position_angle = float(parameters["ROTANGLE"]) * u.deg

        # Extract the slit details
        slit_elements = doc.getElementsByTagName("slit")
        self._slits = []
        for slit in slit_elements:
            ra = float(slit.attributes["xce"].value) * u.deg
            dec = float(slit.attributes["yce"].value) * u.deg
            width = float(slit.attributes["width"].value) * u.arcsec
            height = float(slit.attributes["length"].value) * u.arcsec
            tilt = 0.0 * u.deg
            if "tilt" in slit.attributes.keys():
                tilt = float(slit.attributes["tilt"].value) * u.deg
            self._slits.append(
                MosMaskSlit(
                    center=SkyCoord(ra=ra, dec=dec),
                    width=width,
                    height=height,
                    tilt=tilt,
                )
            )

        # Extract the reference star details
        reference_star_elements = doc.getElementsByTagName("refstar")
        self._reference_stars = []
        for ref in reference_star_elements:
            ra = float(ref.attributes["xce"].value) * u.deg
            dec = float(ref.attributes["yce"].value) * u.deg
            self._reference_stars.append(SkyCoord(ra=ra, dec=dec))
