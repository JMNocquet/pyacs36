"""
Statistics and model evaluation functions for L1-trend analysis.
"""

import numpy as np
from scipy.stats import chi2
import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def get_stats_l1_model(x, y, fy, alpha, component_mask=None):
    """
    Calculate statistics for L1-trend model evaluation.

    Parameters
    ----------
    x : numpy.ndarray
        Input time as 1D array
    y : numpy.ndarray
        Input raw data as 1D array
    fy : numpy.ndarray
        Input l1trend filtered data as 1D array
    alpha : float
        Hyperparameter used for l1 filtering
    component_mask : numpy.ndarray, optional
        Boolean mask indicating which data points should be considered for statistics.
        If None, all points are considered.

    Returns
    -------
    tuple
        (alpha, cchi2, sigma, cp, cchi2_probability, Cp, AICc, BIC)
    """
    # tolerance to detect velocity change
    #tol = .1 # mm/yr note JMN 10/08/2024: this seems to be too tight
    tol = 2. # mm/yr note JMN 10/08/2024: this seems to be too tight

    # Apply component mask if provided
    if component_mask is not None:
        # Use only the data points where the component is active (non-zero)
        active_indices = component_mask
        if np.sum(active_indices) == 0:
            # If no active points, return default values
            return alpha, np.inf, np.inf, 2, 0.0, np.inf, np.inf, np.inf
        
        x_active = x[active_indices]
        y_active = y[active_indices]
        fy_active = fy[active_indices]
        n = y_active.shape[0]
    else:
        x_active = x
        y_active = y
        fy_active = fy
        n = y.shape[0]

    sigma = 1.5 * np.median(np.fabs(np.diff(y_active)))
    cchi2 = np.sum( (fy_active - y_active)**2 / sigma**2 )

    # avoid division by zero
    #if np.sqrt(cchi2/n) < .5:
    #    c_cchi2 = n / 4
    #    WARNING("median of residuals gets too small. Set manually to .5. %.1lf %.1lf" % (cchi2, c_cchi2))
    #    cchi2 = c_cchi2

    dy = np.diff(fy_active) / (np.diff(x_active))
    # number of parameters = number of change points + 2 for start & end
    cp = np.where(np.fabs(np.diff(dy)) > tol )[0].shape[0] + 2
    AIC = n * np.log(cchi2 / n) + 2 * cp
    AICc = n * np.log(cchi2 / n) + 2 * cp + 2 * cp * (cp + 1) / (n - cp - 1)
    BIC = n * np.log(cchi2 / n) + np.log(n) * cp
    cdof = x_active.shape[0] - cp
    cchi2_probability = chi2.cdf(cchi2, cdof) * 100.
    Cp = cchi2 + cp
    DEBUG("alpha: %10f log10_alpha: %6.2lf chi2 : %10.2f fit %8.2lf n_param %05d chi2_proba: %5.2lf Cp: %10.2lf AIC: %10.2lf AICc: %10.2lf BIC: %8.2lf " % (
     alpha, np.log10(alpha), cchi2, np.sqrt(cchi2 / x_active.shape[0]), cp, cchi2_probability, Cp, AIC, AICc, BIC))
    return alpha, cchi2, np.sqrt(cchi2 / x_active.shape[0]), cp, cchi2_probability, Cp, AICc, BIC
