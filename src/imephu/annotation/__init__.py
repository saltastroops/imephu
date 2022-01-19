from typing import Protocol

from astropy.coordinates import Angle
from astropy.visualization.wcsaxes import WCSAxes

from imephu.geometry import SkyPosition, SkyVector


"""A finder chart annotation."""


class Annotation(Protocol):
    """An annotation to add to a finder chart."""

    def add_to(self, ax: WCSAxes) -> None:
        """Add this annotation to a plot

        Parameters
        ----------
        ax: `~astropy.visualization.wcsaxes.WCSAxes`
            Plot axes.
        """
        ...

    def rotate(self, center: SkyPosition, angle: Angle) -> "Annotation":
        """Rotate this annotation around a center and return the result.

        Parameters
        ----------
        center: tuple of (`~astropy.units.Quantity` ["angle"], `~astropy.units.Quantity` ["angle"])
            Point around which to rotate the annotation.
        angle: `~astropy.units.Quantity` ["angle"]
            Angle of rotation. A positive value implies an anti-clockwise rotation.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the rotation.
        """
        ...

    def translate(self, displacement: SkyVector) -> "Annotation":
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: `~imephu.geometry.SkyVector`
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the translation.
        """
        ...
