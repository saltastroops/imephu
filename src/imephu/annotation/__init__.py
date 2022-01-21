from typing import Protocol

from astropy.coordinates import Angle, SkyCoord
from astropy.units import Quantity
from astropy.visualization.wcsaxes import WCSAxes


class Annotation(Protocol):
    """An annotation to add to a finder chart."""

    def add_to(self, ax: WCSAxes) -> None:
        """Add this annotation to a finder chart.

        Parameters
        ----------
        ax: `~astropy.visualization.wcsaxes.WCSAxes`
            Plot axes.
        """
        ...

    def rotate(self, center: SkyCoord, angle: Angle) -> "Annotation":
        """Rotate this annotation around a center and return the result.

        The rotation angle is an angle on the sky, measured from north to east.

        Parameters
        ----------
        center: `~astropy.coordinates.SkyCoord`
            Point around which to rotate the annotation.
        angle: `~astropy.coordinates.Angle`
            Angle of rotation, measured from north towards the east.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the rotation.
        """
        ...

    def translate(self, displacement: Quantity) -> "Annotation":
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the translation.
        """
        ...
