"""
Apply given offsets to a Gts time series.
"""

import copy
import numpy as np

from ._utils import __ensure_list_of_list


def apply_offsets(self, np_offset, opposite=False, in_place=False, verbose=False):
    """Apply given offsets to the time series.

    Parameters
    ----------
    np_offset : numpy.ndarray or list
        1D or 2D array (or list of lists): offset_dates, north, east, up, s_north, s_east, s_up.
    opposite : bool, optional
        If True, apply opposite of offsets. Default is False.
    in_place : bool, optional
        If True, modify self; else return new Gts. Default is False.
    verbose : bool, optional
        Verbose mode. Default is False.

    Returns
    -------
    Gts
        self (if in_place) or new Gts.
    """
    new_Gts = copy.deepcopy(self)

    if not isinstance(np_offset, np.ndarray):
        np_offset = np.array(__ensure_list_of_list(np_offset))
        if np_offset.ndim != 2:
            raise TypeError("!!! ERROR: Could not understand offset parameter type:" % np_offset)
            return self
        else:
            if np_offset.shape[1] < 4:
                raise TypeError("!!! ERROR: Could not understand offset parameter type:" % np_offset)
                return self

    if opposite:
        np_offset[:, 1:4] = -np_offset[:, 1:4]

    new_data = self.data.copy()
    cumulated_offset = np.zeros(new_data.shape)

    for i in np.arange(np_offset.shape[0]):
        date = np_offset[i, 0]
        if verbose:
            print("-- adding %10.5lf %10.2lf %10.2lf %10.2lf (mm)" % (date, np_offset[i, 1]*1000.0, np_offset[i, 2]*1000.0, np_offset[i, 3]*1000.0))
        lindex = np.where(new_data[:, 0] > date)
        cumulated_offset[:, 1][lindex] = cumulated_offset[:, 1][lindex] + np_offset[i, 1]
        cumulated_offset[:, 2][lindex] = cumulated_offset[:, 2][lindex] + np_offset[i, 2]
        cumulated_offset[:, 3][lindex] = cumulated_offset[:, 3][lindex] + np_offset[i, 3]
        cnew_data = new_data + cumulated_offset
        new_Gts.data = cnew_data.copy()
        new_Gts.data_xyz = None

    if in_place:
        self.data = new_Gts.data.copy()
        new_Gts.data_xyz = None
        return self
    else:
        return new_Gts
