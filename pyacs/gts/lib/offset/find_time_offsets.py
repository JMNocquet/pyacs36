"""
Find time of suspected offsets by RMS time series over a sliding window.
"""

import numpy as np


def find_time_offsets(self, option=None, ndays=7, th_detection_rms=3, th_detection_offset=3):
    """Find the time of suspected offsets by RMS time series calculated over ndays."""
    if not option:
        option = 'detrend'

    data_saved = self.copy().data
    outliers_saved = self.outliers
    if self.outliers:
        self.remove_outliers(in_place=True)
    data_rms = self.rms(ndays)

    lindex_offset = []
    for i in (1, 2, 3):
        threshold_rms = th_detection_rms * np.median(data_rms[:, i])
        lindex_offset_i = []
        if max(data_rms[:, i]) > threshold_rms:
            index_out_rms = np.argwhere(data_rms[:, i] > threshold_rms)
            for j in index_out_rms:
                differentiate_ndays = np.hstack(np.fabs(np.diff(self.data[j:(j+ndays), i])))
                index = np.argwhere(differentiate_ndays > th_detection_offset * np.median(differentiate_ndays))
                if len(index) == 1:
                    lindex_offset_i.append(index[0] + j)
            if len(lindex_offset_i) > 0:
                lindex_offset_i = np.hstack(lindex_offset_i)
                for idx in np.unique(lindex_offset_i):
                    if len(np.argwhere(lindex_offset_i == idx)) >= 2:
                        lindex_offset.append(idx)

    if len(lindex_offset) > 0:
        self.offsets = list(self.data[np.unique(lindex_offset), 0])
    else:
        self.offsets = []

    self.data = data_saved
    self.test_offsets()
    self.outliers = outliers_saved
