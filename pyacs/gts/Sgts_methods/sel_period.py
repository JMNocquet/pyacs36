###################################################################
def sel_period(self,period,min_data=2,verbose=True):
###################################################################
    """Select time series that have data in the given period.

    Parameters
    ----------
    period : list
        [start, end] in decimal years.
    min_data : int, optional
        Minimum number of points required to keep a series. Default is 2.
    verbose : bool, optional
        Verbose mode. Default is True.

    Returns
    -------
    Sgts
        New Sgts instance.
    """


    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug



    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # loop on gts
    lsite = []
    for code in sorted( self.lcode() ):
        DEBUG("Testing %s" % code)
        wts = self.__dict__[code].extract_periods( period )
        if wts.data is not None:
            if wts.data.shape[0] >= min_data:
                lsite.append(code)
                DEBUG("Keeping GPS sites:%s " % code )
            else:
                VERBOSE("Less than %d data in requested period: %s" % ( min_data , code ) )
        else:
            VERBOSE("No data for requested period site:%s" %  code )

    MESSAGE("%d over %d sites kept" % ( len(lsite) , len(self.lcode()) ))
    return self.sub(linclude = lsite )            
