"""
Estimate local offset amplitudes using window around each offset date.
"""

import numpy as np

import pyacs.gts.lib.gts_estimators


def estimate_local_offset(self, window_length=4, in_place=False):
    """Estimate local offset amplitudes using window_length positions before and after each offset.

    Returns
    -------
    Gts
        New Gts (or self if in_place) with offsets_values set.
    """
    clean_data = self.remove_outliers().data
    local_offsets_values = np.zeros((len(self.offsets_dates), 7))
    local_offsets_values[:, 0] = self.offsets_dates
    noffset = len(self.offsets_dates)
    list_index_t0_i = np.searchsorted(clean_data[:, 0], local_offsets_values[:, 0])

    index_local = []
    for i in list_index_t0_i:
        index_local.append(list(range(i - window_length, i + window_length)))
    if len(index_local) != 0:
        index_local = np.hstack(index_local)
        index_local = np.delete(index_local, np.argwhere(index_local < 0))
        index_local = np.delete(index_local, np.argwhere(index_local >= len(clean_data[:, 0])))
    else:
        local_offsets_values = None
        return self

    data_local = np.take(clean_data, index_local, axis=0)
    ndate = len(data_local[:, 0])
    for k in range(1, 4):
        A = np.ones((ndate, (1 + noffset)), float)
        for i in range(ndate):
            for j in range(noffset):
                if data_local[i, 0] <= local_offsets_values[j, 0]:
                    A[i, (1 + j)] = 0.
        X, s_X = pyacs.gts.lib.gts_estimators.least_square(A, data_local[:, k])[0:2]
        if X is not None:
            local_offsets_values[:, k] = X[1:]
            local_offsets_values[:, 3 + k] = s_X[1:]
        else:
            local_offsets_values = None

    if in_place:
        self.offsets_values = local_offsets_values
        return self
    else:
        new_Gts = self.copy()
        new_Gts.offsets_values = local_offsets_values
        return new_Gts
