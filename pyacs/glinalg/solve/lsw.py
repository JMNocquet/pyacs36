###############################################################################
def lsw(G,d,std):
###############################################################################
    """
    Solve the least-squares (LS) with data uncertainties provided as a vector
    
    :param G: m x n model matrix as 2D numpy array
    :param d: m 1D numpy observation vector
    :param std: standard deviation vector for d

    :note: the system is modified to be solved by ordinary LS by the change G<- (G.T/std).T and d<- d/std
    """

    # full solution
    
    GG=(G.T/std).T
    dd=d/std

    return ls(GG,dd)
