###################################################################
def append(self, gts):
###################################################################
    """
    Append a Gts to the current Sgts
     
    :param gts: Gts instance to be appended to the current Sgts instance
    """
     
    self.__dict__[gts.code]= gts
