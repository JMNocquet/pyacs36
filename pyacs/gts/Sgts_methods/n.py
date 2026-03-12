###################################################################
def n(self,lexclude=[],linclude=[]):
###################################################################
    """Return the number of Gts in this Sgts (after optional filters).

    Parameters
    ----------
    lexclude : list, optional
        Site codes to exclude. Default is [].
    linclude : list, optional
        If non-empty, only these codes are counted. Default is [].

    Returns
    -------
    int
        Number of time series.
    """
     
    
    return( len(self.lcode( linclude=linclude , lexclude=lexclude ) ))
