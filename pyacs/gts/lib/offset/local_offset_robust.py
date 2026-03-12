"""
Estimate a local offset with a robust method (median of before/after differences).
"""

import numpy as np


def local_offset_robust(self, date, n, verbose=False, debug=False):
    """Estimate a local offset (no velocity) with a robust method.

    Parameters
    ----------
    date : float
        Offset date in decimal year.
    n : int
        Number of dates before/after used in estimation.
    verbose : bool, optional
        Verbose mode. Default is False.
    debug : bool, optional
        Debug output. Default is False.

    Returns
    -------
    numpy.ndarray
        1D array [date, north, east, up, s_north, s_east, s_up].
    """
    tmp_gts = self.extract_ndates_around_date(date, n)
    array_of_offsets_estimates = None
    for i in np.arange(1, n):
        before = np.median(tmp_gts.extract_ndates_before_date(date, i).data[:, 1:4].reshape(-1, 3), axis=0)
        after = np.median(tmp_gts.extract_ndates_after_date(date, i).data[:, 1:4].reshape(-1, 3), axis=0)
        offsets_value = after - before
        if array_of_offsets_estimates is None:
            array_of_offsets_estimates = offsets_value.reshape(-1, 3)
        else:
            array_of_offsets_estimates = np.vstack((array_of_offsets_estimates, offsets_value))
    if debug:
        print('all offsets value')
        print(array_of_offsets_estimates * 1.E3)
    final_values = np.median(array_of_offsets_estimates, axis=0)
    if verbose:
        print("-- combined robust offset estimates (mm): %.1lf %.1lf %.1lf" % (final_values[0]*1.E3, final_values[1]*1.E3, final_values[2]*1.E3))
    return np.array([date, final_values[0], final_values[1], final_values[2], 0.0, 0.0, 0.0])
