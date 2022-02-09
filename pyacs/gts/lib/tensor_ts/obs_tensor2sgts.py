
def obs_tensor2sgts( T_OBS_RAW , np_names_t_obs, np_obs_date_s , verbose = True ):
    """
    Converts a obs_tensor object into a Sgts
    """
    

    # import

    import numpy as np
    from pyacs.gts.Gts import Gts
    from pyacs.gts.Sgts import Sgts
    import pyacs.lib.astrotime as at

    # import pyeq.message.message as MESSAGE
    # import pyeq.message.verbose_message as VERBOSE
    # import pyeq.message.error as ERROR
    # import pyeq.message.debug_message as DEBUG

    # initialize Sgts
    
    sgts = Sgts( read=False)

    # loop on sites
    
    for i in np.arange( np_names_t_obs.shape[0] ):
        
        code = np_names_t_obs[i]
        
        #DEBUG('converting %s ' % code)
        
        # get the index of valid dates
        
        lindex = np.where( np.isfinite( T_OBS_RAW[:,i,0] )  )[0]
    
        data = np.zeros( ( lindex.shape[0] ,10 ) )
        
        # obs_tensor is ENU and Gts are NEU

        data[ : , 2 ] = T_OBS_RAW[ lindex , i  , 0 ] * 1E-3
        data[ : , 1 ] = T_OBS_RAW[ lindex , i  , 1 ] * 1E-3    
        data[ : , 3 ] = T_OBS_RAW[ lindex , i  , 2 ] * 1E-3

        data[ : , 5 ] = T_OBS_RAW[ lindex , i  , 3 ] * 1E-3
        data[ : , 4 ] = T_OBS_RAW[ lindex , i  , 4 ] * 1E-3   
        data[ : , 6 ] = T_OBS_RAW[ lindex , i  , 5 ] * 1E-3
        
        data[:,0] = at.datetime2decyear(  at.seconds2datetime(np_obs_date_s[lindex] ) )

        sgts.append( Gts( code=code , data = data ) )
            
    
    #VERBOSE("converted %d time series " % ( len( sgts.lcode() )) )

    return sgts