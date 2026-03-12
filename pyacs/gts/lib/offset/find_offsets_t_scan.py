"""
Find suspected offsets using t-scan (step detection) statistics.
"""

import numpy as np

from ._utils import __fmt_date


def find_offsets_t_scan(self, threshold=0.8, window=250, in_place=False, lcomponent='NE', verbose=True, debug=True):
    """Find suspected offsets using t-scan step detection."""
    from pyacs.gts.lib import step_detect
    tmp_ts = self.copy()
    lindex_step_north = []
    lindex_step_east = []
    lindex_step_up = []

    if 'N' in lcomponent:
        t_stat = step_detect.t_scan(tmp_ts.data[:, 1], window=window)
        t_stat_max = t_stat.max()
        t_stat /= np.abs(t_stat).max()
        lindex_step_north = step_detect.find_steps(np.abs(t_stat), threshold)
        if debug:
            print("-- North suspected: %d significance: %10.1lf" % (len(lindex_step_north), t_stat_max))

    if 'E' in lcomponent:
        t_stat = step_detect.t_scan(tmp_ts.data[:, 2], window=window)
        t_stat_max = t_stat.max()
        t_stat /= np.abs(t_stat).max()
        lindex_step_east = step_detect.find_steps(np.abs(t_stat), threshold)
        if debug:
            print("-- East  suspected: %d significance: %10.1lf" % (len(lindex_step_east), t_stat_max))

    if 'U' in lcomponent:
        t_stat = step_detect.t_scan(tmp_ts.data[:, 3], window=window)
        t_stat_max = t_stat.max()
        t_stat /= np.abs(t_stat).max()
        lindex_step_up = step_detect.find_steps(np.abs(t_stat), threshold)
        if debug:
            print("-- Up    suspected: %d significance: %10.1lf" % (len(lindex_step_up), t_stat_max))

    if verbose:
        print('-- merging the information and selecting most obvious suspect offsets')

    suspects_all_date = lindex_step_north + lindex_step_east + lindex_step_up
    lindex = list(set(sorted(suspects_all_date)))

    if lindex != []:
        lindex_1 = (np.array(lindex) - 1).tolist()
        loffsets_dates = (tmp_ts.data[lindex, 0] + tmp_ts.data[lindex_1, 0]) / 2.
    else:
        loffsets_dates = []

    if verbose:
        print("-- Suspected offsets (t-statistics) at dates")
        for date in loffsets_dates:
            print("%10.6lf %s" % (date, __fmt_date(date)))

    if in_place:
        self.offsets_dates = loffsets_dates
        return self
    new_Gts = self.copy()
    new_Gts.offsets_dates = loffsets_dates
    return new_Gts
