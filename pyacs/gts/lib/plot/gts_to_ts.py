
########################################################################################################
def gts_to_ts( gts , date=None, unit='mm' , date_unit='decyear' , date_ref=None, set_zero_at_date=None, center= True ):
########################################################################################################
    """
    prepare a gts for plot.
    
    :param gts: Gts instance
    :param date: a list of dates [start_date , end_date] in decimal degrees. Default is None for all dates.
    :param unit: unit for output
    :param date_unit: date unit for the x-axis
    :param date_ref: reference date in decimal year.
    
    :return np_date, y, yerr
    
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
    if date != []:
        #tmp_gts = tmp_gts.extract_periods( pyacs.lib.utils.__ensure_list_of_list( date ) )
        tmp_gts = tmp_gts.extract_periods( pyacs.lib.utils.__ensure_list_of_list( date ), no_reset=True )

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


