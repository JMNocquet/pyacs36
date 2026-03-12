"""
Find suspected offsets using median filter and differentiation.
"""

import numpy as np

from ._utils import __fmt_date, __nargmax


def __print_suspects(D):
    for i in np.arange(D.shape[0]):
        print("%10.6lf %04d %8.1lf " % (D[i, 0], D[i, 1], D[i, 2]))


def suspect_offsets_mf(self, threshold=3, verbose=True, lcomponent='NE', n_max_offsets=5, in_place=False, debug=False):
    """Try to find offsets in a time series using a median filter."""
    n = n_max_offsets
    lwindow = [3, 5, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33]
    suspects_all_date = None

    for wl in lwindow:
        tmp_ts = self.median_filter(wl).differentiate()
        abs_data = np.abs(tmp_ts.data[:, 1:4])
        norm_abs_data = abs_data / np.median(abs_data, axis=0)

        if 'N' in lcomponent:
            suspects_N = __nargmax(norm_abs_data[:, 0], n, threshold)
            if debug:
                print("-- North wl=%d " % wl)
            if suspects_N.size > 0:
                ldate = tmp_ts.data[suspects_N[:, 0].astype(int), 0]
                suspects_N_date = np.insert(suspects_N, 0, ldate, axis=1)
                if suspects_all_date is None:
                    suspects_all_date = suspects_N_date
                else:
                    suspects_all_date = np.vstack((suspects_all_date, suspects_N_date))
                if debug:
                    __print_suspects(suspects_N_date)
            else:
                if debug:
                    print('  -- no suspect found for component North')

        if 'E' in lcomponent:
            suspects_E = __nargmax(norm_abs_data[:, 1], n, threshold)
            if debug:
                print("-- East wl=%d " % wl)
            if suspects_E.size > 0:
                ldate = tmp_ts.data[suspects_E[:, 0].astype(int), 0]
                suspects_E_date = np.insert(suspects_E, 0, ldate, axis=1)
                if suspects_all_date is None:
                    suspects_all_date = suspects_E_date
                else:
                    suspects_all_date = np.vstack((suspects_all_date, suspects_E_date))
                if debug:
                    __print_suspects(suspects_E_date)
            else:
                if debug:
                    print('  -- no suspect found for component East')

        if 'U' in lcomponent:
            suspects_U = __nargmax(norm_abs_data[:, 2], n, threshold)
            np.save('test.npy', norm_abs_data[:, 2])
            if debug:
                print("-- Up wl=%d " % wl)
            if suspects_U.size > 0:
                ldate = tmp_ts.data[suspects_U[:, 0].astype(int), 0]
                suspects_U_date = np.insert(suspects_U, 0, ldate, axis=1)
                if suspects_all_date is None:
                    suspects_all_date = suspects_U_date
                else:
                    suspects_all_date = np.vstack((suspects_all_date, suspects_U_date))
                if debug:
                    __print_suspects(suspects_U_date)
            else:
                if debug:
                    print('  -- no suspect found for component Up')

    if verbose:
        print('-- merging the information and selecting most obvious suspect offsets')

    if suspects_all_date is not None:
        suspects_all_date = suspects_all_date[np.argsort(-suspects_all_date[:, 2])]
        lindex = []
        for i in np.arange(suspects_all_date.shape[0]):
            if int(suspects_all_date[i, 1]) not in lindex:
                lindex.append(int(suspects_all_date[i, 1]))
                if len(lindex) > n:
                    break
        loffsets_dates = tmp_ts.data[lindex, 0]
    else:
        loffsets_dates = []

    if verbose:
        print('-- Suspected offsets (median filter - differentiation) at dates ')
        for date in loffsets_dates:
            print("%10.6lf %s" % (date, __fmt_date(date)))

    if in_place:
        self.offsets_dates = loffsets_dates
        return self
    new_Gts = self.copy()
    new_Gts.offsets_dates = loffsets_dates
    return new_Gts
