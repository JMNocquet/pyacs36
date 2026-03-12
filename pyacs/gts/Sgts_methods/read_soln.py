###################################################################
def read_soln(self,soln,verbose=True):
###################################################################
    """Read IGS soln file and set offsets_dates from solution (P) changes.

    Parameters
    ----------
    soln : str
        Path to IGS soln.snx file.
    verbose : bool, optional
        Verbose mode. Default is True.

    Returns
    -------
    None
        Modifies self in place.

    Notes
    -----
    In-place: populates offsets_dates for each Gts from soln discontinuities.
    """
    
    from pyacs.sol.discontinuity import Discontinuities
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG



    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    fsoln=Discontinuities()
    VERBOSE("Reading %s" % soln)
    fsoln.read_igs_discontinuity(soln)
    
    for gts in self.lGts():
        VERBOSE("Adding offsets for %s" % gts.code)
        if gts.code in fsoln.lsite_offsets_dates:
            gts.offsets_dates=fsoln.lsite_offsets_dates[gts.code]

    return()
