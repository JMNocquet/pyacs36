def compute_common_mode_l1trend(self, already_ivel=True, max_ivel=40):
    """
    This approach estimates a common mode from instantaneous velocity

    parameter already_ivel: boolean, if True ts_ref are already results from l1trend
    """


    # import
    from icecream import ic

    from pyacs.gts.lib.tensor_ts.sgts2obs_tensor import sgts2tensor
    from pyacs.gts.Gts import Gts
    import numpy as np
    import pyacs.lib.astrotime as at

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # if not ivel, then computes ivel
    if not already_ivel:
        MESSAGE("Computing instantaneous velocity")
        ts_ref = self.gts('disp2vel')
    else:
        ts_ref = self

    # convert to observation tensor
    MESSAGE("Converting to observation tensor")
    T, np_name, np_s = sgts2tensor(ts_ref, rounding='day')
    MESSAGE("Read %d time series" % T.shape[1])

    # set nan unrealistic velocities
    for i in np.arange( T.shape[1] ):
        # detect too high velocities on East & North only
        np_idx = np.where( np.fabs(T[:,i,:1])>max_ivel )[0]
        T[np_idx,i,:] = np.nan



    # center velocity
    T = T - np.nanmedian(T, axis=0)

    # removes rows with only NaN
    lidx = ~np.isnan(T[:,:,0]).all(axis=1)
    T = T[lidx]
    np_s = np_s[lidx]


    # computes the median velocity at every date
    MESSAGE("Computing instantaneous median velocity")
    np_median_velocity = np.nanmedian(T / 365.25, axis=1)

    # interpolate to avoid nan
    from scipy import interpolate

    for i in np.arange(3):
        lidx = np.where(~np.isnan(np_median_velocity[:, i]))[0]

        if lidx.size > 0:
            f = interpolate.interp1d(np_s[lidx], np_median_velocity[lidx, i])
            np_median_velocity[:, i] = f(np_s)

    # do the integration to get a time series in mm
    MESSAGE("Integrating velocity to displacement")
    ts_cmm = np.cumsum(np_median_velocity, axis=0) * 1.E-3

    # creates a time series
    cmm = Gts(code='CMM_')
    cmm.data = np.zeros((ts_cmm.shape[0], 10))
    cmm.data[:, 0] = at.datetime2decyear(at.seconds2datetime(np_s))
    cmm.data[:, 1] = ts_cmm[:, 1]
    cmm.data[:, 2] = ts_cmm[:, 0]
    cmm.data[:, 3] = ts_cmm[:, 2]

    cmm.data[:, 4:] = 1.E-3

    # detrend
    MESSAGE("Removing trend")

    vel = cmm.detrend_seasonal().velocity
    cmm_detrend = cmm.remove_velocity(vel)

    # return

    return cmm_detrend

