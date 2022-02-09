"""
Time series express as a 4D-numpy array D and a separate observation time vector T

    D(i,j,k) would be the displacement
    observation time index i
    for site j
    for component k [0,1,2,3,4,5] [de,dn,du,sde,sdn,sdu]

    no data are NaN

"""
###################################################################
def sgts2tsr( sgts , tol=0.01, verbose=False ):
###################################################################
        
    """
    Convert an Sgts object into a tsr
    
    :param sgts: Sgts object
    :param tol : tolerance in decimal day to assign the same date 
    
    :return: NAMES (1D-Numpy string array), DATES (1D-Numpy integer array of seconds), a 4D numpy array
    D(i,j,k) would be the displacement
    observation time index i
    for site j
    for component k [0,1,2,3,4,5] [de,dn,du,sde,sdn,sdu]
    
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
    