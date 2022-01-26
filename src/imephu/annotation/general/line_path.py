from copy import deepcopy
from typing import Any, Sequence, Tuple, Union

from astropy.coordinates import Angle, SkyCoord
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import WCS
from matplotlib.patches import PathPatch
from matplotlib.path import Path

from imephu.annotation import Annotation
from imephu.geometry import rotate, sky_position_to_pixel, translate


class LinePathAnnotation(Annotation):
    """An annotation consisting of straight lines between vertices.

    The path displayed by this annotation is defined by a set of vertices, which are
    connected by straight lines on the finder chart. (The connecting lines in general
    will refer to curved lines on the sky.)

    By default the path is assumed to be closed, i.e. the last vertex is connected to
    the first by a line. This behavior may be changed with the ``closed`` parameter.

    Parameters
    ----------
    vertices: sequence of `~astropy.coordinates.SkyCoord`
        The path vertices, as an array or sequence of sky coordinates (with right
        ascension abd declination).
    wcs: `~astropy.wcs.WCS`
        WCS object.
    closed: bool, default: True
        Whether the the path is closed, i.e. whether the last vertex is connected to the
        first by a line.
    edgecolor: color, default: "black"
        The edge color.
    facecolor: color, default: "none"
        The filling color.
    **kwargs: dict, optional
        Additional keyword arguments, which will be passed to Matplotlib's
        `~matplotlib.patches.PathPatch` patch constructor.
    """

    def __init__(
        self,
        vertices: Sequence[SkyCoord],
        wcs: WCS,
        closed: bool = True,
        edgecolor: Union[str, Tuple[float, float, float]] = "black",
        facecolor: Union[str, Tuple[float, float, float]] = "none",
        **kwargs: Any,
    ):
        if closed:
            # Matplotlib will ignore the last vertex for a closed path; we avoid this by
            # having that vertex twice
            self._vertices = list(vertices)[:] + [vertices[-1]]
        else:
            self._vertices = list(vertices)
        self._wcs = wcs
        self._closed = closed
        self._kwargs = deepcopy(kwargs)
        self._kwargs["edgecolor"] = edgecolor
        self._kwargs["facecolor"] = facecolor

    def add_to(self, ax: WCSAxes) -> None:
        """Add the line path to a finder chart.

        Parameters
        ----------
        ax: `~astropy.visualization.wcsaxes.WCSAxes`
            Plot axes.
        """
        vertices_px = [
            sky_position_to_pixel(vertex, self._wcs) for vertex in self._vertices
        ]
        path = Path(vertices_px, closed=self._closed)  # noqa
        path_patch = PathPatch(path, **self._kwargs)
        ax.add_patch(path_patch)

    def rotate(self, pivot: SkyCoord, angle: Angle) -> "LinePathAnnotation":
        """Rotate this annotation around a pivot and return the result.

        The rotation angle is an angle on the sky, measured from north to east.

        Parameters
        ----------
        pivot: `~astropy.coordinates.SkyCoord`
            Point around which to rotate the annotation.
        angle: `~astropy.coordinates.Angle`
            Angle of rotation, measured from north towards the east.

        Returns
        -------
        `~imephu.annotation.general.LinePathAnnotation`
            The annotation resulting from the rotation.
        """
        rotated_annotation = deepcopy(self)
        rotated_vertices = [
            rotate(vertex, pivot, angle, self._wcs) for vertex in self._vertices
        ]
        rotated_annotation._vertices = rotated_vertices
        return rotated_annotation

    def translate(self, displacement: Angle) -> "LinePathAnnotation":
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.general.LinePathAnnotation`
            The annotation resulting from the translation.
        """
        translated_annotation = deepcopy(self)
        translated_vertices = [
            translate(vertex, displacement) for vertex in self._vertices
        ]
        translated_annotation._vertices = translated_vertices
        return translated_annotation
