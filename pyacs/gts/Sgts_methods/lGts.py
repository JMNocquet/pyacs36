###################################################################
def lGts(self,lexclude=[],linclude=[]):
###################################################################
    """
    Returns the list of Gts in the current Sgts
    :param lexclude: list of sites to be excluded
    :param linclude: list of sites to be included, excluding all other.
     
    """
         
    lsite=self.lcode(lexclude=lexclude,linclude=linclude)
    Lgts=[]
    for site in sorted(lsite):Lgts.append(self.__dict__[site])
    return( Lgts )
