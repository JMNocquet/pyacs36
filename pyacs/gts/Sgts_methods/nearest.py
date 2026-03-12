def nearest(self, site, n=1):
    """Return the n nearest sites to the given site.

    Parameters
    ----------
    site : str
        Site code to search from.
    n : int, optional
        Number of nearest sites to return. Default is 1.

    Returns
    -------
    list
        List of site codes (nearest first), or None if site not in Sgts.
    """

    # import
    import numpy as np
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug

    # check
    if site not in self.lcode():
        ERROR("site %s not in the current Sgts object" % site)
        return None

    # get the distance matrix
    Dm = self.make_distance_matrix_from_sgts()
    np_sorted_lcode = np.array(sorted(self.lcode()), dtype=str)
    idx = np.argwhere(np_sorted_lcode == site)[0]
    sorted_arg = np.argsort(Dm[idx, :])[0]

    return np_sorted_lcode[sorted_arg][:n + 1]
