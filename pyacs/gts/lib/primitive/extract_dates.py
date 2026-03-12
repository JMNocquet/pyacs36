###################################################################
def extract_dates(self,dates, tol= 0.05 , rounding=None, in_place=False, verbose=True):
###################################################################
    """
    Return a time series extracted for a given list of dates.

    Parameters
    ----------
    dates : list or ndarray
        Dates as a list or 1D numpy array of decimal dates.
    tol : float, optional
        Date tolerance in days to assert that two dates are equal (default 0.05 day).
    rounding : str, optional
        Rounding for date comparison ('second', 'minute', 'hour', 'day'). Inferred from tol if None.
    in_place : bool, optional
        If True, will make change in place; if False, returns a new time series.
    verbose : bool, optional
        Verbose mode.
    """

    # import 
    import inspect
    import numpy as np
    import pyacs.gts
    import pyacs.lib.astrotime as at

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        ERROR( error )
        return( self )

    # check that data_xyz and data are consistent
    if self.data_xyz is not None:
        if self.data_xyz.shape != self.data.shape:
            ERROR(".data and .data_xyz are not consistent for time series %s. Returning time series" % self.code)
            return self.copy()

    # handles rounding
    if rounding is None:
        rounding = 'day'
        if tol < 0.05:
            rounding = 'hour'
        if tol < 0.0007:
            rounding = 'minute'
        if tol < 1.E-5:
            rounding = 'second'

    # converts dates

    np_date_s = at.decyear2seconds( self.data[:,0],rounding=rounding)
    np_date_s_ref = at.decyear2seconds( np.array(dates),rounding=rounding)

    # find date index

    xy, np_idx, y_ind = np.intersect1d(np_date_s, np_date_s_ref, return_indices=True)

    if np_idx.size ==0:
        WARNING("No data available for time series %s at requested dates. Check rounding. Returning input time series" % self.code )
        return self.copy()

    # working gts
    new_gts = self.copy()

    # case .data_xyz
    
    if new_gts.data_xyz is not None:
        new_gts.data_xyz = self.data_xyz[np_idx,:]
        # re-generate NEU time series
        new_gts.xyz2neu(corr=False)

    else:
        new_gts.data = self.data[np_idx, :]

    # handles outliers
    new_gts.outliers = []

    # handles offsets_date
    new_gts.offsets_dates = []

    return(new_gts)
