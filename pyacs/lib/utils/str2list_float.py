"""String parsing helpers."""


def str2list_float(my_str):
    """Convert a string representation of a list to a list of floats.

    Parameters
    ----------
    my_str : str
        String like '[0, 2.5, 1E9]'.

    Returns
    -------
    list of float
        Parsed values, e.g. [0, 2.5, 1e9].
    """
    return list(map(float, my_str.replace('[', '').replace(']', '').split(',')))

