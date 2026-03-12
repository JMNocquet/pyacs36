"""
Find suspected offsets from day-to-day differences (differentiate time series).
"""

from ._utils import get_suspected_dates


def suspect_offsets(self, threshold=3, verbose=True, lcomponent='NE', n_max_offsets=10, in_place=False):
    """Try to find offsets in a time series from day-to-day differences."""
    import numpy as np
    dates = self.data[:, 0]
    diff_data = np.diff(self.data, n=1, axis=0)
    new_dates = dates[:-1] + diff_data[:, 0] / 2.
    diff_data[:, 0] = new_dates

    lldates = get_suspected_dates(diff_data, threshold, verbose=verbose)

    if len(lldates) > n_max_offsets:
        if verbose:
            print('-- Too many offsets detected. Keeping only the largest for analysis and will see later the others...')
        while len(lldates) > n_max_offsets:
            threshold = threshold * 1.5
            lldates = get_suspected_dates(diff_data, threshold)

    if in_place:
        self.offsets_dates = lldates
        return self
    new_Gts = self.copy()
    new_Gts.offsets_dates = lldates
    return new_Gts
