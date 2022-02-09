###################################################################
def lcode(self,lexclude=[],linclude=[]):
###################################################################
    """
    Returns the list of Gts codes in the current Sgts
    exclude is a list of code to be excluded
     
    :param lexclude: list of sites to be excluded
    :param linclude: list of sites to be included, excluding all other.
     
    """
    # import
    import re 
     
     
    lcode=[]
    if linclude != []:
        for k in list(self.__dict__.keys()):
            if ((re.match('[A-Z0-9_]{4}',k))) and (k not in lexclude) and (k in linclude):lcode.append(k)
    else:
        for k in list(self.__dict__.keys()):
            if ((re.match('[A-Z0-9_]{4}',k))) and (k not in lexclude):lcode.append(k)
    
    return(lcode)
