"""
Functions for checking and analyzing L1-trend model quality.
"""

import numpy as np
import pyacs.lib.astrotime as at

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def check_l1_trend(ts, l1ts, component='ENU', min_samples_per_segment=4, threshold_bias_res_detect=40, threshold_vel=8, plot=False):
    """
    Inspect the result from l1trend model of a time series.
    Returns a list periods ([sdate,edate]) where bad modelling is suspected.

    Parameters
    ----------
    ts : pyacs.gts.Gts.Gts
        Raw time series
    l1ts : pyacs.gts.Gts.Gts
        L1-trend time series
    component : str
        Components to be analyzed. Default: 'EN'
    min_samples_per_segment : int
        Minimum number of samples for a segment to be inspected (default: 6)
    threshold_bias_res_detect : float
        Threshold to detect bias residuals (default: 40)
    threshold_vel : float
        Instantaneous velocity threshold for a segment to be considered in detection (default: 8 mm/yr)
    plot : bool
        Plot the results (default: False)

    Returns
    -------
    tuple
        (H_period, H_cp, H_cp_pb)
        - H_period: Dictionary of suspicious periods for each component
        - H_cp: Dictionary of breakpoint indices for each component
        - H_cp_pb: Dictionary of problematic breakpoint indices for each component
    """

    # initialize
    H_period = {}
    H_cp = {}
    H_cp_pb = {}


    VERBOSE("Searching suspicious periods for %s" % (l1ts.code))

    # define n_median_filter
    n_median_filter = np.min([11, ts.data.shape[0]-1])
    if n_median_filter % 2 == 0:
        n_median_filter = n_median_filter - 1

    # East component
    if 'E' in component:
        VERBOSE('Component E')
        t = ts.data[:, 0]
        y = ts.data[:, 2] * 1.E3
        l1y = l1ts.data[:, 2] * 1.E3
        fy = ts.median_filter(n=n_median_filter).data[:, 2] * 1.E3
        lidx = make_node_ts(t, l1y)
        H_cp['E'] = lidx
        VERBOSE("Searching suspicious periods for component E ")
        H_period['E'], H_cp_pb['E'] = compute_measure_bias(t, l1y, fy - l1y, lidx,
                                                           min_samples_per_segment=min_samples_per_segment,
                                                           threshold_bias_res_detect=threshold_bias_res_detect,
                                                           threshold_vel=threshold_vel)
        VERBOSE("#%d periods suspected for component E" % len(H_period['E']))
    
    # North component
    if 'N' in component:
        VERBOSE('Component N')
        t = ts.data[:, 0]
        y = ts.data[:, 1] * 1.E3
        l1y = l1ts.data[:, 1] * 1.E3
        fy = ts.median_filter(n=n_median_filter).data[:, 1] * 1.E3
        lidx = make_node_ts(t, l1y)
        H_cp['N'] = lidx
        VERBOSE("Searching suspicious periods for component N ")
        H_period['N'], H_cp_pb['N'] = compute_measure_bias(t, l1y, fy - l1y, lidx,
                                                           min_samples_per_segment=min_samples_per_segment,
                                                           threshold_bias_res_detect=threshold_bias_res_detect,
                                                           threshold_vel=threshold_vel)
        VERBOSE("#%d periods suspected for component N" % len(H_period['N']))

    # Up component
    if 'U' in component:
        VERBOSE('Component U')
        t = ts.data[:, 0]
        y = ts.data[:, 3] * 1.E3
        l1y = l1ts.data[:, 3] * 1.E3
        fy = ts.median_filter(n=n_median_filter).data[:, 3] * 1.E3
        lidx = make_node_ts(t, l1y)
        H_cp['U'] = lidx
        VERBOSE("Searching suspicious periods for component U ")
        H_period['U'], H_cp_pb['U'] = compute_measure_bias(t, l1y, fy - l1y, lidx,
                                                           min_samples_per_segment=min_samples_per_segment,
                                                           threshold_bias_res_detect=threshold_bias_res_detect,
                                                           threshold_vel=threshold_vel)
        VERBOSE("#%d periods suspected for component U" % len(H_period['U']))

    if plot:
        ts.plot(superimposed=l1ts, lperiod=H_period)

    return H_period, H_cp, H_cp_pb


def measure_bias(t, y):
    """
    Provide a measure of potential bias in a residual time series (0-100).

    Parameters
    ----------
    t : ndarray
        Time array.
    y : ndarray
        Residual values.

    Returns
    -------
    tuple
        (res, resv): bias score and median absolute residual.
    """
    nsample = y.shape[0]
    mid_date = (t[-1] + t[0]) / 2

    lidx_1 = np.where(t < mid_date)[0]
    lidx_2 = np.where(t > mid_date)[0]

    n_positif1 = np.where(y[lidx_1] >= 0)[0].shape[0]
    n_negatif2 = np.where(y[lidx_2] <= 0)[0].shape[0]

    n_negatif1 = lidx_1.shape[0] - n_positif1
    n_positif2 = lidx_2.shape[0] - n_negatif2

    res = np.max([np.fabs(((n_positif1 + n_negatif2) / nsample * 100.) - 50.) * 2,
                  np.fabs(((n_negatif1 + n_positif2) / nsample * 100.) - 50.) * 2])

    resv = np.median(np.fabs(y))

    return res, resv


def make_node_ts(x, y, threshold=.5):
    """
    Find breakpoints assuming y is piecewise linear.

    Parameters
    ----------
    x : ndarray
        Time/abscissa.
    y : ndarray
        Values (piecewise linear).
    threshold : float, optional
        Threshold on change in slope to declare breakpoint.

    Returns
    -------
    ndarray
        Indices of breakpoints in the time series.
    """
    dy = np.diff(y) / (np.diff(x))
    cp = np.where(np.fabs(np.diff(dy)) > threshold)[0]
    return cp + 1


def compute_measure_bias(t, l1y, y, cp, min_samples_per_segment=6, threshold_bias_res_detect=40, threshold_vel=8,
                         verbose=True):
    """
    Compute fit indicators from residual (obs minus piecewise linear).

    Parameters
    ----------
    t : ndarray
        Time array.
    l1y : ndarray
        L1-trend values.
    y : ndarray
        Residual (obs - piecewise_linear).
    cp : ndarray
        Breakpoint indices.
    min_samples_per_segment : int, optional
        Minimum samples per segment (default 6).
    threshold_bias_res_detect : float, optional
        Threshold for bias detection (default 40).
    threshold_vel : float, optional
        Velocity threshold in mm/yr (default 8).
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    tuple
        (list of suspicious periods [sdate, edate], list of problematic breakpoint indices).
    """
    datetime_t = at.decyear2datetime(t)
    rms_l1 = np.median(np.fabs(y))

    # computes average velocity
    ivel = np.diff(l1y) / np.diff(t)
    median_vel = np.median(ivel)
    median_vel_change = np.median(np.fabs(ivel - median_vel))
    lperiod = []
    lcp = []
    for i in np.arange(cp.shape[0] - 1):
        wt = t[cp[i]:cp[i + 1] + 1]
        wy = y[cp[i]:cp[i + 1] + 1]
        if wt.shape[0] < min_samples_per_segment:
            continue

        mb, resv = measure_bias(wt, wy)

        # threshold for velocity 4 * median scatter
        rvel = np.fabs(np.mean(np.diff(l1y[cp[i]:cp[i + 1] + 1]) / np.diff(wt)) - median_vel) / median_vel_change

        if mb > threshold_bias_res_detect and rvel > threshold_vel:
            VERBOSE((" %s-%s nsample: %3d velocity_bias_measure: %.1lf rms_l1: %.1lf vs %.1lf rvel %.1lf " % \
                      (datetime_t[cp[i]].isoformat()[:10], datetime_t[cp[i + 1]].isoformat()[:10], wt.shape[0], mb,
                       resv, rms_l1, rvel)))

            lperiod.append([t[cp[i]], t[cp[i + 1]]])
            lcp.append(cp[i])
    return lperiod, lcp
