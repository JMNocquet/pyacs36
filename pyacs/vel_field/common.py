"""Common-site and matching methods for Velocity_Field."""


def common(self, linGpoint, prefit=10.0, lexclude=[]):
    """Return sites common to this velocity field and a list of GMT_Points.

    Parameters
    ----------
    linGpoint : list
        List of GMT_Point instances (e.g. from SINEX).
    prefit : float, optional
        Maximum coordinate difference (m) to consider a match. Default is 10.0.
    lexclude : list, optional
        Site codes to exclude.

    Returns
    -------
    list
        GMT_Point instances from this field that match by code (coordinates from this field).
    """
    loutGpoint = []
    lcode = self.lcode()
    for M in linGpoint:
        if M.code not in lcode:
            continue
        lGpoint = self.subset(lonly=[M.code])
        for N in lGpoint.sites:
            if (M.spherical_distance(N) > prefit):
                print("! Bad prefit (threshold value %8.3lf) for %s: %10.1lf" % (prefit, M.code, M.xyz_distance(N)))
                break
            else:
                loutGpoint.append(N)
    return loutGpoint
