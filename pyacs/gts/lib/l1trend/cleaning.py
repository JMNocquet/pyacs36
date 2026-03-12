"""
Cleaning functions for an already L1-trend filtered time series.
"""

import numpy as np
#from pyacs.gts.Gts import Gts

import logging
import pyacs.message.message as MESSAGE
import pyacs.message.verbose_message as VERBOSE
import pyacs.message.error as ERROR
import pyacs.message.warning as WARNING
import pyacs.message.debug_message as DEBUG


def clean_l1trend(self, raw_gts, threshold='auto'):
    """
    Cleans breakpoints from an already L1-trend filtered time series.
    Removes breakpoints that are too close in value (<0.5 mm/yr).
    Returns a new Gts object interpolated on the dates from l1trend_gts.

    Parameters
    ----------
    raw_gts : pyacs.gts.Gts.Gts
        Original (position) time series.
    threshold : float or 'auto'
        Threshold for cleaning. If 'auto', computed from median difference.

    Returns
    -------
    pyacs.gts.Gts.Gts
        Cleaned/interpolated Gts object.
    """

    if isinstance(threshold, str) and threshold == 'auto':
        # Compute threshold from median difference (in mm)
        diff = raw_gts.data[:,1:4] - self.data[:,1:4]
        threshold = np.median(np.abs(diff), axis=0) * 1E3  # mm

    # Get breakpoints using tol='auto' and threshold
    bp = self.l1trend_to_breakpoints(tol='auto', threshold=threshold)

    # Interpolate cleaned breakpoints on l1trend_gts dates
    new_data = np.zeros_like(self.data)
    new_data[:,0] = self.data[:,0]
    for i, comp in enumerate(['N', 'E', 'U']):
        # Interpolate breakpoints on l1trend_gts dates (no cleaning)
        bp_dates = bp[comp][0]
        bp_values = bp[comp][1]
        new_data[:,i+1] = np.interp(new_data[:,0], bp_dates, bp_values)
    # Copy uncertainties from l1trend_gts
    new_data[:,4:8] = self.data[:,4:8]

    # Create new Gts object
    cleaned_gts = self.copy()
    cleaned_gts.data = new_data

    return cleaned_gts
