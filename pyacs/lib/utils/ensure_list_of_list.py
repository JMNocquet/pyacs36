"""Utilities for list normalization."""


def __ensure_list_of_list(ll):
    """Ensure argument is a list of lists.

    [a, b] becomes [[a, b]]; [[a, b]] is returned unchanged.

    Parameters
    ----------
    ll : list
        List or list of lists.

    Returns
    -------
    list
        List of lists.

    Raises
    ------
    TypeError
        If ll is not a list.
    """
    if not isinstance(ll, list):
        raise TypeError('!!! __ensure_list_of_list requires a list or a list of list as argument: ', ll)

    if not isinstance(ll[0], list):
        return [ll]
    else:
        return ll

