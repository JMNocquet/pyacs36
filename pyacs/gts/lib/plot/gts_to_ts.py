
########################################################################################################
def gts_to_ts( gts , date=None, unit='mm' , date_unit='decyear' , date_ref=None, set_zero_at_date=None, center= True ):
########################################################################################################
    """Prepare Gts for plotting (dates, y, yerr).

    Parameters
    ----------
    gts : Gts
        Gts instance.
    date : list, optional
        [start_date, end_date] in decimal year; None = all dates.
    unit : str, optional
        Output unit (e.g. 'mm'). Default is 'mm'.
    date_unit : str, optional
        X-axis date unit. Default is 'decyear'.
    date_ref : float, optional
        Reference date in decimal year. Default is None.
    set_zero_at_date : float, optional
        Date to set zero. Default is None.
    center : bool, optional
        Center series. Default is True.

    Returns
    -------
    np_date : array
        Dates for x-axis.
    y : array
        Values.
    yerr : array
        Errors.
    """
    
    ###########################################################################
    # import
    ###########################################################################

    import pyacs.lib.astrotime
    import pyacs.lib.utils
    import numpy as np
    
    ###########################################################################
    # handle dates
    ###########################################################################
    tmp_gts = gts.copy( )
    
    # case date extraction
    #if date != []:
        #tmp_gts = tmp_gts.extract_periods( pyacs.lib.utils.__ensure_list_of_list( date ) )
    #    tmp_gts = tmp_gts.extract_periods( pyacs.lib.utils.__ensure_list_of_list( date ), no_reset=True )

    # case set_zero_at_date
    if set_zero_at_date is not None:
        tmp_gts = tmp_gts.set_zero_at_date( set_zero_at_date )
    
    # case 'days'
    if date_unit == 'days':
        if (date_ref is None) or (date_ref == 0):
            np_date = pyacs.lib.astrotime.decyear2mjd(tmp_gts.data[:,0]) - pyacs.lib.astrotime.decyear2mjd(tmp_gts.data[0,0])
        else:
            np_date = pyacs.lib.astrotime.decyear2mjd(tmp_gts.data[:,0]) - pyacs.lib.astrotime.decyear2mjd( date_ref )
    
    # case 'decyear'
    if date_unit == 'decyear':
        if (date_ref is None) or (date_ref == 0):
            np_date = tmp_gts.data[:,0]
        else:
            np_date = tmp_gts.data[:,0] - date_ref
    
    # case 'cal'
    if date_unit == 'cal':
        np_date = pyacs.lib.astrotime.decyear2datetime( tmp_gts.data[:,0] )
    
    
    ###########################################################################
    # handle unit
    ###########################################################################

    if unit == 'cm':
        tmp_gts.data[:,1:] = tmp_gts.data[:,1:] * 100.

    if unit == 'mm':
        tmp_gts.data[:,1:] = tmp_gts.data[:,1:] * 1000.

    ###########################################################################
    # center
    ###########################################################################
    
    if center:
        tmp_gts.data[:,1:4] = tmp_gts.data[:,1:4] - np.median( tmp_gts.data[:,1:4] , axis=0 )
    
    return np_date , tmp_gts.data[:,1:]


