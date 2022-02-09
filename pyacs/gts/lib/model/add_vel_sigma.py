###################################################################
def add_vel_sigma(self, in_place=False, b_fn=4, verbose=True):
###################################################################
    """
    calculates realistic sigma on velocity components assuming white &
    flicker using eq (19) & (23) from Williams (J. of Geodesy, 2003)
    b_fn is the value for flicker noise, taken as 4 mm/yr^1/4
    model can be detrend, detrend_annual, detrend_seasonal
    if in_place = True then replace the current time series
    """

    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect

    # white noise
    # white noise is estimated from the time series

    if not isinstance(self.velocity, np.ndarray):
        print(
            "!!!ERROR: Can't estimate velocity sigma before velocity because a residual time series is required to estimate noise components")
        print("!!!ERROR: use detrend, detrend_annual or detrend_seasonal first.")
        return (None)

    n = self.data.shape[0]
    if n > 5:
        (a_n, a_e, a_u) = np.std(np.diff(self.data, n=1, axis=0), axis=0)[1:4]
    else:
        (a_n, a_e, a_u) = (2.0, 2.0, 5.0)

    def sigma2_wn(ti, a):
        ti = ti - ti[0]
        n = ti.shape[0]
        sigma2_wn = n * a ** 2 / (n * np.sum(ti ** 2) - (np.sum(ti)) ** 2)
        return (sigma2_wn)

    # flicker noise
    def sigma2_fn(ti, b_fn):

        n = ti.shape[0]
        delta_t = np.mean(np.diff(ti, n=1))

        sigma2_fn = 9 * b_fn ** 2 / (16 * delta_t ** 2 * (n ** 2 - 1))
        #        print 'duration ',delta_t**2 * n**2
        return (sigma2_fn)

    ti = self.data[:, 0]

    if verbose:
        print('fn: ', sigma2_fn(ti, b_fn))

    new_ts = self.copy()
    min_rms = np.min(np.array([a_n, a_e, a_e]))

    sigma_vn = np.sqrt(sigma2_wn(ti, a_n) + sigma2_fn(ti, b_fn)) * (a_n / min_rms)
    sigma_ve = np.sqrt(sigma2_wn(ti, a_e) + sigma2_fn(ti, b_fn)) * (a_e / min_rms)
    sigma_vu = np.sqrt(sigma2_wn(ti, a_u) + sigma2_fn(ti, b_fn)) * (a_u / min_rms)
    # to meters
    new_ts.velocity[3] = sigma_vn / 1000.0
    new_ts.velocity[4] = sigma_ve / 1000.0
    new_ts.velocity[5] = sigma_vu / 1000.0

    if in_place:
        self.data = new_ts.data.copy()
    return (new_ts)
