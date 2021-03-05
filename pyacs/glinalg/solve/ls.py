###############################################################################
def ls(G,d, verbose=False):
###############################################################################
    """
    Solve the least-squares (LS) problem m so that  (Gm-d).T (Gm-d) is minimum.
    
    :param G: m x n model matrix as 2D numpy array
    :param d: m 1D numpy observation vector
    :param verbose: verbose mode
    :return: x,chi2: m (1D numpy array of dim m), chi2 (chi-square)
    :note: solved through numpy.linalg.lstsq
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
