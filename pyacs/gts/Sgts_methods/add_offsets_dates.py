###################################################################
def add_offsets_dates(self,dates,verbose=False):
###################################################################
    """
    add_offsets_dates to every Gts in current Sgts
    """

    # import
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    New_Sgts=Sgts(read=False)
    
    for gts in self.lGts():
        if verbose:print("-- Processing ",gts.code)
        try:
            new_gts=gts
            new_gts.offsets_dates=dates
        except (RuntimeError, TypeError, NameError):
            print("!!! Error processing ",gts.code)
            continue
        if isinstance(new_gts,Gts):
            New_Sgts.append(new_gts)
        else:
            print("!!! Error processing ",gts.code, "!!! No time series created.")
    
    return( New_Sgts )
