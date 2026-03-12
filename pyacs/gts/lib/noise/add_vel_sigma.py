###################################################################
def add_vel_sigma(self, b_fn=4, verbose=True):
###################################################################
    """
    Estimate velocity uncertainties from residuals (white + flicker noise).

    Velocity and residual time series must exist: run detrend(), detrend_annual(),
    or detrend_seasonal() first. Uncertainty on velocity components (N, E, U) is
    computed using white noise from the residual scatter and flicker noise
    (Williams, 2003, eq. 19 and 23). Returns a new Gts with velocity sigmas set.

    Parameters
    ----------
    b_fn : float, optional
        Flicker noise parameter in mm/yr^0.25. Default is 4.
    verbose : bool, optional
        If True, print flicker noise variance. Default is True.

    Returns
    -------
    Gts or None
        New Gts with velocity uncertainties set in ``velocity[3:6]`` (m/yr),
        or None if ``self.velocity`` is not set.

    Notes
    -----
    White noise is estimated from the first differences of the residual time
    series. The combined white + flicker sigma is scaled by component RMS so
    that the minimum component (N, E, or U) keeps the nominal flicker level.

    References
    ----------
    Williams, S. D. P. (2003). The effect of coloured noise on the uncertainties
    of rates from geodetic time series. Journal of Geodesy, 76(9-10), 483-494.
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
    import pyacs.debug

    

    # white noise
    # white noise is estimated from the time series

    if not isinstance(self.velocity, np.ndarray):
        ERROR("Cannot estimate velocity uncertainties: Gts.velocity is not set. "
              "Run detrend(), detrend_annual(), or detrend_seasonal() first to obtain velocity and residuals for noise estimation.")
        return None

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

    return new_ts
