###################################################################
##  MODEL A TIME SERIES WITH A SPLINE
###################################################################

def spline(self, smoothing=1, degree=5, date=None):
    """
    :param smoothing: Positive smoothing factor used to choose the number of knots. Number of knots will be increased
    until the smoothing condition is satisfied:
        sum((w[i] * (y[i]-spl(x[i])))**2, axis=0) <= s
    :param degree: Degree of the smoothing spline. Must be <= 5. Default is k=3, a cubic spline.
    :param date: 1D array of interpolation dates in decimal year, or 'day' for every day. defualt None will interpolate
    at data date only.
    :return: new gts instance
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect

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

    # import
    from scipy.interpolate import UnivariateSpline
    import numpy as np

    Spline_ts = self.copy()
    s_e = UnivariateSpline(self.data[:, 0], self.data[:, 1], s=smoothing * 1E-3, k=degree)
    s_n = UnivariateSpline(self.data[:, 0], self.data[:, 2], s=smoothing * 1E-3, k=degree)
    s_u = UnivariateSpline(self.data[:, 0], self.data[:, 3], s=smoothing * 1E-3, k=degree)
    if date is None:
        Spline_ts.data[:, 1] = s_e(self.data[:, 0])
        Spline_ts.data[:, 2] = s_n(self.data[:, 0])
        Spline_ts.data[:, 3] = s_u(self.data[:, 0])
    if isinstance(date, np.ndarray):
        Spline_ts.data[:, 1] = s_e(date)
        Spline_ts.data[:, 2] = s_n(date)
        Spline_ts.data[:, 3] = s_u(date)
    if date == 'day':
        import pyacs.lib.astrotime as at
        np_date = at.mjd2decyear(np.arange(at.decyear2mjd(self.data[0, 0]), at.decyear2mjd(self.data[-1, 0])))
        Spline_ts.data = np.zeros((np_date.shape[0], 10))
        Spline_ts.data[:, 0] = np_date
        Spline_ts.data[:, 1] = s_e(np_date)
        Spline_ts.data[:, 2] = s_n(np_date)
        Spline_ts.data[:, 3] = s_u(np_date)

    return (Spline_ts)
