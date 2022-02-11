###################################################################
def detrend(self, method='L2', in_place=False, periods=[], exclude_periods=[]):
    ###################################################################
    """
    detrends a time series and save velocity estimates in velocity attribute

    :param periods         : periods used to estimate the velocity
    :param exclude_periods : periods to be excluded for the velocity estimate
    :param in_place        : if True then replace the current time series
    :return                : the detrended time series
    :note                  : outliers from Gts.outliers are omitted in the estimation and
    offsets given Gts.offsets_dates are estimated simultaneously

    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect

    # after this method .data  and .data_xyz are not consistent so .data_xyz is set to None
    self.data_xyz = None

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print(error)
        return (self)
    ###########################################################################

    import copy
    outliers = copy.deepcopy(self.outliers)
    tmp_ts = self.remove_outliers()

    if periods != []:
        tmp_ts = tmp_ts.extract_periods(periods)
    if exclude_periods != []:
        tmp_ts = tmp_ts.exclude_periods(periods)

    detrended = tmp_ts.make_model(option='detrend', method=method)

    vel = detrended.velocity
    offsets_values = detrended.offsets_values

    new_gts = self.copy()
    new_gts.outliers = outliers
    new_gts.offsets_values = offsets_values
    new_gts.velocity = vel

    model = new_gts.mmodel()

    new_gts.data[:, 1:4] = new_gts.data[:, 1:4] - model.data[:, 1:4]

    if in_place:
        self.data = new_gts.data
        return (self)
    else:
        return (new_gts)
