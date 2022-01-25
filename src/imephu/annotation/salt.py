from astropy.wcs import WCS

from imephu.annotation.general import TextAnnotation


def title(
    target: str, proposal_code: str, pi_family_name: str, wcs: WCS
) -> TextAnnotation:
    """Return an annotation for adding a title to the finder chart.

    The title is of the form "target name (proposal code; surname of the Principal
    Investigator)", such as "Magrathea (2022-1-SCI-042; Adams)".

    Parameters
    ----------
    target: str
        The name of the target for which the finder chart is intended.
    proposal_code: str
        The proposal code for which the finder chart is intended.
    pi_family_name: str
        The family name of the Principal Investigator.
    wcs: `~astropy.wcs.WCS`
        WCS object.

    Returns
    -------
    `~imephu.annotation.general.TextAnnotation`
    """
    title_text = f"{target} ({proposal_code}; {pi_family_name})"
    return TextAnnotation(
        (0.5, 1.03),
        title_text,
        wcs=wcs,
        color="black",
        horizontalalignment="center",
        style="italic",
        weight="bold",
        size="large",
    )
