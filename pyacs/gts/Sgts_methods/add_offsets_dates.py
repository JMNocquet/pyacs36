###################################################################
def add_offsets_dates(self,dates,verbose=False):
###################################################################
    """
    add_offsets_dates to every Gts in current Sgts
    """

    # import
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts

    
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
