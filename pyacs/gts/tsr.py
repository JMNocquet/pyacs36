"""
Time series express as a 4D-numpy array D and a separate observation time vector T

    D(i,j,k) would be the displacement
    observation time index i
    for site j
    for component k [0,1,2,3,4,5] [de,dn,du,sde,sdn,sdu]

    no data are NaN

"""
class tsr:

###################################################################
    def __init__ (self):
###################################################################
        pass

    @classmethod        
###################################################################
    def convert(cls, sgts, tol=0.01, verbose=False):
###################################################################
        
        """
        Convert an Sgts object into a tsr
        
        :param sgts: Sgts object
        :param tol : tolerance in decimal day to assign the same date 
        """

        # import
        import numpy as np
        

