###################################################################
def lGts(self,lexclude=[],linclude=[]):
###################################################################
    """Return list of Gts instances in this Sgts (after optional filters).

    Parameters
    ----------
    lexclude : list, optional
        Site codes to exclude. Default is [].
    linclude : list, optional
        If non-empty, only these codes are included. Default is [].

    Returns
    -------
    list
        List of Gts instances (sorted by code).
    """
         
    lsite=self.lcode(lexclude=lexclude,linclude=linclude)
    Lgts=[]
    for site in sorted(lsite):Lgts.append(self.__dict__[site])
    return( Lgts )
