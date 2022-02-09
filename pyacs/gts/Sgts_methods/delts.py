###################################################################
def delts(self,code):
###################################################################
    """
    Delete a time series from an Sgts instance
    
    :param code: code to be excluded
    """
    
    delattr(self, code)
