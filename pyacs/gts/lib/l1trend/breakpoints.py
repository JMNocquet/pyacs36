"""
Breakpoint detection and extraction from L1-trend filtered time series.
"""

import numpy as np
#from pyacs.gts.Gts import Gts

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def l1trend_to_breakpoints(self, tol='auto', threshold=[1.,1., 5.]):
    """
    Convert a Gts resulting from a L1-trend-filtering to a dictionary of breakpoints.
    The breakpoints are computed by looking for significant changes in the slope of the time series.

    Parameters
    ----------
    tol : float or 'auto'
        Tolerance for detecting breakpoints in mm/yr. If 'auto', finds the greatest tol
        such that the max difference between interpolated breakpoints and the time series is below a threshold for each component.
    threshold : list of floats
        List of thresholds in mm for each component ('N', 'E', 'U') to determine the best tolerance if tol is 'auto'.
        Default is [1., 1., 5.] for 'N', 'E', and 'U' respectively.
        If tol is a float, this parameter is ignored.
    Returns
    -------
    bp : record
        Dictionary with keys as component names ('E', 'N', 'U') and dates/values as
        bp[component][0], H_bp[component][1]
    """

    lcomponent = ['N', 'E', 'U']
    bp  = {comp: [[],[]] for comp in lcomponent}   
    time_decyear = self.data[:,0]

    def compute_bp(tol_vals, lcomponent=['N', 'E', 'U']):
        bp_local = {comp: [[], []] for comp in lcomponent}
        licomponent = np.array([['N', 'E', 'U'].index(comp) for comp in lcomponent])
        for i, icomponent in enumerate(licomponent):
            y = self.data[:, icomponent + 1]
            dy = np.diff(y) / (np.diff(time_decyear))
            cp = np.where(np.fabs(np.diff(dy)) > tol_vals[i])[0]
            cp = cp + 1
            bp_local[lcomponent[i]][0] = time_decyear[cp]
            bp_local[lcomponent[i]][1] = y[cp]
            bp_local[lcomponent[i]][0] = np.insert(bp_local[lcomponent[i]][0], 0, time_decyear[0])
            bp_local[lcomponent[i]][1] = np.insert(bp_local[lcomponent[i]][1], 0, y[0])
            bp_local[lcomponent[i]][0] = np.append(bp_local[lcomponent[i]][0], time_decyear[-1])
            bp_local[lcomponent[i]][1] = np.append(bp_local[lcomponent[i]][1], y[-1])
        return bp_local

    if tol == 'auto':
        tol_candidates = np.array([0.5, 1, 2, 3, 4, 5, 10, 20, 50, 100])
        best_tol = [tol_candidates[0]] * 3
        for i, component in enumerate(lcomponent):
            for t in tol_candidates:
                tol_vals = [0.0005, 0.0005, 0.0005]
                tol_vals[i] = t / 1000.0
                bp_test = compute_bp(tol_vals, lcomponent=[lcomponent[i]])
                bp_times = bp_test[component][0]
                bp_values = bp_test[component][1]
                icomponent = np.array(['N', 'E', 'U'].index(component)) + 1
                values = self.data[:,icomponent]
                interp_values = np.interp(time_decyear, bp_times, bp_values)
                max_mm = np.max(np.abs(interp_values - values)) * 1E3
                if max_mm < threshold[i]:
                    best_tol[i] = t
                else:
                    break
        tol_vals = [t / 1000.0 for t in best_tol]
        bp = compute_bp(tol_vals)
    else:
        tol_val = float(tol) / 1000.0
        bp = compute_bp([tol_val, tol_val, tol_val])
    return bp
