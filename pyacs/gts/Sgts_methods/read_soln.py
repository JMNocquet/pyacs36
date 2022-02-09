###################################################################
def read_soln(self,soln,verbose=True):
###################################################################
    """
    read a IGS soln file and add an offsets_dates for any P change in soln.
     
    :param soln: soln.snx IGS file
     
    :note: the method is in place
    """
    
    from pyacs.sol.discontinuity import Discontinuities
     
    fsoln=Discontinuities()
    if verbose:print("-- Reading ",soln)
    fsoln.read_igs_discontinuity(soln)
    
    for gts in self.lGts():
        if verbose:print("-- Adding offsets to ",gts.code)
        if gts.code in fsoln.lsite_offsets_dates:
            gts.offsets_dates=fsoln.lsite_offsets_dates[gts.code]
    
    return()
