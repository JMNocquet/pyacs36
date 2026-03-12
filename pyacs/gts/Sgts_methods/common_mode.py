###################################################################
def common_mode(self, lref=[] , detrend_method='detrend_median', method='median' , center=True, verbose=True):
###################################################################
    """Compute and remove common mode from time series.

    Parameters
    ----------
    lref : list, optional
        Site codes for reference (common mode). Default is [].
    detrend_method : str, optional
        'detrend_median' or 'detrend' for reference series. Default is 'detrend_median'.
    method : str, optional
        'median' or 'mean' for common mode. Default is 'median'.
    center : bool, optional
        If True, center the stack. Default is True.
    verbose : bool, optional
        Verbose mode. Default is True.

    Returns
    -------
    Sgts
        New Sgts with filtered series and _CMM for the common mode.

    Notes
    -----
    Assumes daily time series.
    """


    # import  
    import pyacs.gts.lib.tensor_ts.sgts2obs_tensor
    import pyacs.gts.lib.tensor_ts.obs_tensor2sgts
    import numpy as np
    from pyacs.gts.Gts import Gts
    import pyacs.lib.astrotime as at
    import pyacs.lib.coordinates as coor

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect

    from icecream import ic

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # check lref
    if lref==[]:
        ERROR("lref argument must include at least one site.", exit=True)

    # detrend for the reference sites

    VERBOSE("detrending using method %s" % detrend_method)
    if detrend_method is not None:
        dts = self.sub(linclude=lref).gts( detrend_method )
    else:
        dts = self.sub(linclude=lref).copy()

    # filter out None
    ndts = dts.delnone()

    # creates an obs_tensor instance from the detrend time series
    VERBOSE("Converting to obs tensor")
    T_OBS_RAW_REF , np_names_t_obs_ref, np_obs_date_s_ref = pyacs.gts.lib.tensor_ts.sgts2obs_tensor.sgts2tensor( ndts, rounding='day' , verbose=verbose )
    
    # center time series
    VERBOSE("Centering time series")
    if center:
        T_OBS_RAW_REF[:,:,0] = T_OBS_RAW_REF[:,:,0] - np.nanmedian( T_OBS_RAW_REF[:,:,0] , axis=0 )
        T_OBS_RAW_REF[:,:,1] = T_OBS_RAW_REF[:,:,1] - np.nanmedian( T_OBS_RAW_REF[:,:,1] , axis=0 )
        T_OBS_RAW_REF[:,:,2] = T_OBS_RAW_REF[:,:,2] - np.nanmedian( T_OBS_RAW_REF[:,:,2] , axis=0)
    
    # get the index of sites used for the common mode
    
    #lidx_code = []
    
    #for i in np.arange(np_names_t_obs.shape[0]):
    #    if np_names_t_obs[i] in lref:
    #         lidx_code.append(i)
    
    #np_idx_code = np.array( sorted( lidx_code ) )
    
    # compute common mode

    VERBOSE("Computing common mode")

    #T_OBS_RAW_REF = T_OBS_RAW[:, np_idx_code, :]
    
    if method == 'median':
        CMM = np.nanmedian(T_OBS_RAW_REF,axis=1)
    if method == 'mean':
        CMM = np.nanmean(T_OBS_RAW_REF,axis=1)

    # creates an obs_tensor instance from the original time series
    
    T_OBS_RAW , np_names_t_obs, np_obs_date_s = pyacs.gts.lib.tensor_ts.sgts2obs_tensor.sgts2tensor( self, rounding='day' , verbose=verbose )
    
    # check that dates are the same
    if not np.all( np_obs_date_s == np_obs_date_s_ref ):
        WARNING("The list of reference sites for common mode does not provide a common mode estimate at every date")
        ERROR("Try interpolation", exit=True)

    
    # remove common mode


    VERBOSE("Removing common mode")

    T_OBS_RAW[:,:,0] = T_OBS_RAW[:,:,0] - CMM[:,0].reshape(-1,1)
    T_OBS_RAW[:,:,1] = T_OBS_RAW[:,:,1] - CMM[:,1].reshape(-1,1)
    T_OBS_RAW[:,:,2] = T_OBS_RAW[:,:,2] - CMM[:,2].reshape(-1,1)
    
    # dispersion of the common mode
    # TBD
    
    
    # converts obs_tensor object to Sgts


    VERBOSE("Creating filtered Sgts")

    filtered_sgts = pyacs.gts.lib.tensor_ts.obs_tensor2sgts.obs_tensor2sgts( T_OBS_RAW , np_names_t_obs, np_obs_date_s , verbose=verbose)
    
    # adds the common mode time series
    
    cmm_ts = Gts( code='_CMM' )
    
    # get index where values are Nan
    
    lindex = np.where( np.isfinite( CMM[:,0] )  )[0]
        
    data = np.zeros( ( lindex.shape[0] , 10 ) )
    
    # obs_tensor is ENU and Gts are NEU
    data[ : , 2 ] = CMM[ lindex , 0 ] * 1E-3
    data[ : , 1 ] = CMM[ lindex , 1 ] * 1E-3    
    data[ : , 3 ] = CMM[ lindex , 2 ] * 1E-3
    # set std as 1 mm to avoid singularity
    data[ : , 4 ] =  1E-3
    data[ : , 5 ] =  1E-3    
    data[ : , 6 ] =  1E-3
    
    
    data[:,0] = at.datetime2decyear( at.seconds2datetime(np_obs_date_s[lindex] ) )

    # populate X0, Y0, Z0, lon, lat, he
    X_cmm = 0
    Y_cmm = 0
    Z_cmm = 0
    for code in filtered_sgts.lcode():
        filtered_sgts.__dict__[code].X0 = self.__dict__[code].X0
        filtered_sgts.__dict__[code].Y0 = self.__dict__[code].Y0
        filtered_sgts.__dict__[code].Z0 = self.__dict__[code].Z0

        filtered_sgts.__dict__[code].lon = self.__dict__[code].lon
        filtered_sgts.__dict__[code].lat = self.__dict__[code].lat
        filtered_sgts.__dict__[code].h   = self.__dict__[code].h

        X_cmm = X_cmm + self.__dict__[code].X0
        Y_cmm = Y_cmm + self.__dict__[code].Y0
        Z_cmm = Z_cmm + self.__dict__[code].Z0


    cmm_ts.data = data
    
    X_cmm = X_cmm / len( filtered_sgts.lcode() )
    Y_cmm = Y_cmm / len( filtered_sgts.lcode() )
    Z_cmm = Z_cmm / len( filtered_sgts.lcode() )
    
    cmm_ts.X0 = X_cmm / len( filtered_sgts.lcode() )
    cmm_ts.Y0 = Y_cmm / len( filtered_sgts.lcode() )
    cmm_ts.Z0 = Z_cmm / len( filtered_sgts.lcode() )
    
    ( cmm_ts.lon , cmm_ts.lat , cmm_ts. h ) = coor.xyz2geo(cmm_ts.X0, cmm_ts.Y0, cmm_ts.Z0, unit='dec_deg')
    
    filtered_sgts.append( cmm_ts )
    
    return filtered_sgts