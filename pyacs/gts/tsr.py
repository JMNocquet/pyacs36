"""
Time series as 4D array D and observation time vector T.

D(i,j,k): displacement at time i, site j, component k [de,dn,du,sde,sdn,sdu].
No-data entries are NaN.
"""

class tsr:
    """Time series representation (4D array + dates)."""

###################################################################
    def __init__ (self):
###################################################################
        pass

    @classmethod        
###################################################################
    def convert(cls, sgts, tol=0.01, verbose=False):
###################################################################
        """
        Convert an Sgts object into a tsr.

        Parameters
        ----------
        sgts : Sgts
            Sgts object.
        tol : float, optional
            Tolerance in decimal day to assign the same date.
        verbose : bool, optional
            Verbose mode.

        Returns
        -------
        tsr
            tsr instance with NAMES, DATES, D.
        """

        # import
        import numpy as np
        

