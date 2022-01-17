from typing import Protocol

from astropy.coordinates import Angle
from astropy.visualization.wcsaxes import WCSAxes

from imephu.geometry import Point, Vector


"""A finder chart annotation."""


class Annotation(Protocol):
    """An annotation to add to a finder chart."""

    def add_to(self, ax: WCSAxes):
        """Add this annotation to plot axes.

        Parameters
        ----------
        ax: `~astropy.visualization.wcsaxes.WCSAxes`
            Plot axes.
        """
        ...

    def rotate(self, center: Point, angle: Angle) -> "Annotation":
        """Rotate this annotation around a center and return the result.

        Parameters
        ----------
        center: `~imephu.geometry.Point`
            Point around which to rotate the annotation.
        angle: `~astropy.coordinates.Angle`
            Angle of rotation. A positive value implies an anti-clockwise rotation.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the rotation.
        """
        ...

    def translate(self, displacement: Vector) -> "Annotation":
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: `~imephu.geometry.Vector`
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.Annotation`
            The annotation resulting from the translation.
        """
        ...
