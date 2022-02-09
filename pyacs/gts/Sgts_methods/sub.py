###################################################################
def sub(self,lexclude=[],linclude=[]):
###################################################################
    """
    Returns a new Sgts instance excluding Gts with code in lexclude and keeping Gts with code in include
    
    :param lexclude: list of sites to be excluded
    :param linclude: list of sites to be included, excluding all other.
     
   """

    from pyacs.gts.Sgts import Sgts
             
    if linclude != []:
        sub_Sgts = Sgts(read=False)
        for code in linclude:
            if ( code not in lexclude) and ( self.has_ts( code ) ):
                sub_Sgts.append(self.__dict__[code].copy())
    else:
        sub_Sgts = self.copy()
        for code in lexclude:
            if self.has_ts( code ):
                sub_Sgts.delts(code)
     
    return( sub_Sgts )
