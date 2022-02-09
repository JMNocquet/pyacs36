###################################################################
def n(self,lexclude=[],linclude=[]):
###################################################################
    """
    Returns the number of Gts codes in the current Sgts
    exclude is a list of code to be excluded
     
    :param lexclude: list of sites to be excluded
    :param linclude: list of sites to be included, excluding all other.
    """
     
    
    return( len(self.lcode( linclude=linclude , lexclude=lexclude ) ))
