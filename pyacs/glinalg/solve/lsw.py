###############################################################################
def lsw(G,d,std):
###############################################################################
    """
    Solve least-squares with data uncertainties given as a vector.

    Parameters
    ----------
    G : ndarray
        m x n model matrix (2D).
    d : ndarray
        m observation vector (1D).
    std : ndarray
        Standard deviation vector for d (length m).

    Returns
    -------
    ndarray
        Solution from ordinary LS on G<- (G.T/std).T, d<- d/std.

    Notes
    -----
    System is weighted by 1/std and solved via ls.
    """

    # full solution
    
    GG=(G.T/std).T
    dd=d/std

    return ls(GG,dd)
