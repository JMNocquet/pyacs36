###################################################################
## detrend_seasonal
###################################################################

def detrend_seasonal(self, method='L2', periods=None, exclude_periods=None):
    """
    Estimate trend + annual + semi-annual terms and remove them from the time series.

    Velocity, annual and semi_annual are saved in Gts.velocity, Gts.annual, Gts.semi_annual.

    Parameters
    ----------
    method : str, optional
        Estimation method (e.g. 'L2').
    periods : list, optional
        Periods used for estimation.
    exclude_periods : list, optional
        Periods to exclude from estimation.

    Returns
    -------
    Gts
        Detrended time series (new instance).

    Notes
    -----
    Outliers (Gts.outliers) are omitted; offsets (Gts.offsets_dates) are estimated simultaneously.
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # after this method .data  and .data_xyz are not consistent so .data_xyz is set to None
    #self.data_xyz = None

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        ERROR(error)
        return (self)
    ###########################################################################

    import copy
    outliers = copy.deepcopy(self.outliers)
    tmp_ts = self.remove_outliers()

    if periods:
        tmp_ts = tmp_ts.extract_periods(periods)
    if exclude_periods:
        tmp_ts = tmp_ts.exclude_periods(periods)

    detrended = tmp_ts.make_model(option='detrend_seasonal', method=method)

    vel = detrended.velocity
    offsets_values = detrended.offsets_values
    annual = detrended.annual
    semi_annual = detrended.semi_annual

    new_gts = self.copy()
    new_gts.outliers = outliers
    new_gts.offsets_values = offsets_values
    new_gts.velocity = vel
    new_gts.annual = annual
    new_gts.semi_annual = semi_annual

    model = new_gts.mmodel()

    new_gts.data[:, 1:4] = new_gts.data[:, 1:4] - model.data[:, 1:4]

    return (new_gts)
