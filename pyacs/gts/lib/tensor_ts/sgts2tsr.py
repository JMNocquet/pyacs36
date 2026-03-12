"""
Time series as 4D array D and observation time vector T.

D(i,j,k): displacement at time index i, site j, component k [0,1,2,3,4,5] = [de,dn,du,sde,sdn,sdu].
No-data entries are NaN.
"""
###################################################################
def sgts2tsr( sgts , tol=0.01, verbose=False ):
###################################################################
        
    """
    Convert an Sgts object into a tsr (4D array + dates).

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
    tuple
        (NAMES, DATES, D). NAMES: 1D string array; DATES: 1D integer array (seconds);
        D(i,j,k): displacement at time i, site j, component k [de,dn,du,sde,sdn,sdu].
    """

    # import
    import numpy as np
        
    
    # get the list of all dates
 
    l_dec_year = []
    
    for gts in sgts.lGts():
        l_dec_year = l_dec_year + gts.data[:,0].tolist()

    
    # work on the dates
    
    np_dec_year = np.unique( np.array( l_dec_year ) )
    
    # eliminate dates that are less than tol
    
    tol_in_dec_year = tol / 365.25
    
    np_diff_dec_year = np.sqrt( np.diff( np_dec_year )**2 )
    
    print( np.min( np_diff_dec_year ) )
    print( np.max( np_diff_dec_year ) )
    
    if np.min( np_diff_dec_year ) < tol_in_dec_year:
        print("-- date conflict to handle")
    else:
        print("-- dates OK.")
    