###################################################################
def copy(self):
###################################################################
    """
    makes a (deep) copy of the Sgts object
     
    :return : a new Sgts instance
    """
    
    import copy
     
    new_Sgts=copy.deepcopy( self )
    return( new_Sgts )
