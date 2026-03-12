###############################################################################
def sgts2tensor( sgts, rounding='day' , verbose=False ):
###############################################################################
    """
    Return a 3D tensor of all Gts (.data ENU) and site codes.

    Parameters
    ----------
    sgts : Sgts
        Sgts instance.
    rounding : str, optional
        'day', 'hour', or 'second'; dates rounded accordingly (default 'day').
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    tuple
        (T_OBS, np_names, np_date_s). T_OBS[k, i, 0] is East at k-th date for site i (mm).

    Notes
    -----
    T_OBS[k, i, 0] returns the East value at the k-th date for site i in mm.
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime as at

    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.debug_message as DEBUG

    from icecream import ic

# np print option for debug
    #np.set_printoptions(precision=2 , suppress=True)
    
    # np_names
    np_names = np.array( sorted(sgts.lcode() ) )
    
    
    # get all dates
    np_seconds_sol = np.array([],dtype=np.int64)
    for i in np.arange( np_names.shape[0] ):
        # code
        code = np_names[i]
        # get the gts for current site
        wts = sgts.__dict__[code]
        try:
            # update np of all dates        
            np_seconds_sol = np.unique( np.sort( np.append( np_seconds_sol, at.decyear2seconds( wts.data[:,0] , rounding=rounding )) ))
        except:
            ERROR("there is a date problem in gts: %s" % code)
            ERROR("Most probably, it has no data", exit=True)

    # initialize T
    T = np.full( (np_seconds_sol.shape[0] , np_names.shape[0], 6  ) , np.nan )
    
    #DEBUG(("shape of T: %s " % (T.shape,)))
    
    # loop on gts in sgts
    
    for i in np.arange( np_names.shape[0] ):

        # code
        code = np_names[i]
        #DEBUG("%04d / %s " % (i,code) )
        
        # get the gts for current site
        wts = sgts.__dict__[code]

        # date of the current wts in seconds
        np_seconds_site = at.decyear2seconds( wts.data[:,0] , rounding=rounding )
        
        # check for non duplicated dates
        if np_seconds_site.shape[0] > 1:
            if np.min( np.diff( np_seconds_site) ) <= 0:
                #ERROR("there is a date problem in gts: %s" % code)
                #ERROR("most probably, your round option (%s) leads to two successive equal dates" % (rounding))
                #ERROR(("%d" % np.min( np.diff( np_seconds_site) ) ))
                print(np.diff( np_seconds_site) )
                print( np_seconds_site )
                #ERROR("",exit=True)

        # current gts date index in T
        
        lindex= np.intersect1d(np_seconds_sol, np_seconds_site, assume_unique=True, return_indices=True)[1]

        #DEBUG("number of observation in ts: %s:%d" % (code, wts.data.shape[0]) )
        #DEBUG("number of observation in T : %d" % (T.shape[0]) )
        #DEBUG("number of observation to be filled in T from lindex : %d" % (lindex.shape[0]) )
        #DEBUG("size of T slice to be filled: %s" % (T[lindex,i,:].shape,) )
            
        # fill T - ENU - SE SN SU
        T[lindex,i,0] = wts.data[:,2]*1000.
        T[lindex,i,1] = wts.data[:,1]*1000.
        T[lindex,i,2] = wts.data[:,3]*1000.

        T[lindex,i,3] = wts.data[:,5]*1000.
        T[lindex,i,4] = wts.data[:,4]*1000.
        T[lindex,i,5] = wts.data[:,6]*1000.
        
    # return
    return T , np_names , np_seconds_sol

