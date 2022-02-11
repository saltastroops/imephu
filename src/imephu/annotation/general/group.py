from typing import List, Sequence

from astropy.coordinates import Angle, SkyCoord
from astropy.visualization.wcsaxes import WCSAxes

from imephu.annotation import Annotation


class GroupAnnotation(Annotation):
    """
    An annotation made of a group of annotations.

    The annotation comprises a list of annotation items. Any of the
    `~impehu.annotation.Annotation` methods are applied to all the items in the list,
    in the order the items appear in the list. So, for example, the ``apply_to`` method
    calls the ``apply_to`` method for all list items, starting from the first item in
    the list and finishing with the last item.

    Annotations can be added to the group when calling the constructor or by using the
    ``add_item`` method afterwards.

    Parameters
    ----------
    items: sequence of `~imephu.annotation.Annotation`
        The initial list of group items.
    """

    def __init__(self, items: Sequence[Annotation]):
        self._items: List[Annotation] = list(items)

    def add_item(self, item: Annotation) -> None:
        """Add an annotation to this group.

        Parameters
        ----------
        item: `~imephu.annotation.Annotation`
            Annotation to add to this group.
        """
        self._items.append(item)

    def add_to(self, ax: WCSAxes) -> None:
        """Add this annotation to a finder chart.

        Parameters
        ----------
        ax: `~astropy.visualization.wcsaxes.WCSAxes`
            Plot axes.
        """
        for item in self._items:
            item.add_to(ax)

    def rotate(self, pivot: SkyCoord, angle: Angle) -> "GroupAnnotation":
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
        `~imephu.annotation.general.GroupAnnotation`
            The annotation resulting from the rotation.
        """
        rotated_items = [item.rotate(pivot, angle) for item in self._items]
        rotated_annotation = GroupAnnotation(items=rotated_items)
        return rotated_annotation

    def translate(self, displacement: Angle) -> "GroupAnnotation":
        """
        Move this annotation along a displacement vector and return the result.

        Parameters
        ----------
        displacement: 2D array of angles
            The displacement by which to move the annotation.

        Returns
        -------
        `~imephu.annotation.general.GroupAnnotation`
            The annotation resulting from the translation.
        """
        translated_items = [item.translate(displacement) for item in self._items]
        translated_annotation = GroupAnnotation(translated_items)
        return translated_annotation
