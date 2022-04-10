from copy import deepcopy

from astropy.coordinates import Angle, SkyCoord
from astropy.visualization.wcsaxes import WCSAxes

from imephu.annotation import Annotation


class EmptyAnnotation(Annotation):
    """An annotation for adding nothing to a finder chart.

    In other words, this annotation just leaves the finder chart as is.

    This is a convenience class for cases where you formally need to pass an annotation,
    but there is nothing to annotate with.
    """

    def add_to(self, ax: WCSAxes) -> None:
        """Add nothing to a finder chart, leaving it as is.

        ax: `~astropy.visualization.wcsaxes.WCSAxes`
            WCS axes object.
        """
        pass

    def rotate(self, pivot: SkyCoord, angle: Angle) -> "Annotation":
        """Return a copy of this annotation.

        Parameters
        ----------
        pivot: `~astropy.coordinates.SkyCoord`
            Point around which to rotate the annotation.
        angle: `~astropy.coordinates.Angle`
            Angle of rotation, measured from north to east.

        Returns
        -------
        `~imephu.annotation.general.CircleAnnotation`
            A copy of this annotation.
        """
        return deepcopy(self)

    def translate(self, displacement: Angle) -> "Annotation":
        """
        Return a copy of this annotation.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.general.CircleAnnotation`
            A copy of this annotation.
        """
        return deepcopy(self)
