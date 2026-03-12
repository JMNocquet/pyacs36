###############################################################################
def ls(G,d, verbose=False):
###############################################################################
    """
    Solve the least-squares problem min (Gm-d).T (Gm-d).

    Parameters
    ----------
    G : ndarray
        m x n model matrix (2D).
    d : ndarray
        m observation vector (1D).
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    ndarray
        Solution m (1D, length n). chi2 available from lstsq internally.

    Notes
    -----
    Solved via numpy.linalg.lstsq.
    """

    # numpy linalg lstsq
    
    import numpy.linalg

    (m,chi2,rank,s)=numpy.linalg.lstsq(G,d,rcond=-1)

    if verbose:
        print('-- ls info:')
        print('-- chi2: ', chi2 )
        print('-- rank: ', rank )
        print('-- s   : ', s )

    return(m)
