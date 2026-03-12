###################################################################
def sel_rectangle(self,bounds,verbose=True):
###################################################################
    """Select time series for sites inside a rectangular bounds.

    Parameters
    ----------
    bounds : list
        [lon_min, lon_max, lat_min, lat_max] in decimal degrees.
    verbose : bool, optional
        Verbose mode. Default is True.

    Returns
    -------
    Sgts
        New Sgts instance.
    """


    # import    
    from pyacs.gts.Sgts import Sgts
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug
    from tqdm import tqdm


    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    [lon_min,lon_max,lat_min,lat_max]=bounds
    
    new_Sgts=Sgts(read=False)
    for gts in tqdm( self.lGts() ):
        current_lon=gts.lon
        current_lat=gts.lat
         
        if current_lon>=lon_min and current_lon<=lon_max<=lon_max and current_lat>= lat_min and current_lat<=lat_max:
            VERBOSE("%s selected" % gts.code )
            new_Sgts.append(gts)

    return(new_Sgts)
