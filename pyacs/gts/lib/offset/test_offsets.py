"""
Test offsets: delete small ones, F-ratio test, re-check small offsets.
"""

import numpy as np


def test_offsets(self, verbose=False, debug=True, window_length=None):
    """Test offsets: delete small, F-ratio test, re-check small offsets."""
    ouliers_saved = self.outliers
    self.find_outliers_by_residuals()
    self.offsets = self.delete_small_offsets(self.offsets)
    if self.offsets is not None:
        if not window_length:
            window_length = 180
        data_saved = self.copy().data
        data = self.remove_outliers().data
        t0 = self.offsets
        list_index = np.hstack([0, np.searchsorted(data[:, 0], t0), len(data[:, 0])])
        t0_tested = []
        for i in range(len(t0)):
            if verbose:
                print("  -> Testing offset time %10.3lf" % t0[i])
            if window_length:
                if (list_index[i+1] - list_index[i] < window_length/2):
                    index_start = list_index[i]
                else:
                    index_start = list_index[i+1] - window_length/2
                if (list_index[i+2] - list_index[i+1] < window_length/2):
                    index_end = list_index[i+2]
                else:
                    index_end = list_index[i+1] + window_length/2
            else:
                index_start = list_index[i]
                index_end = list_index[i+2]
            data_i = data[index_start:index_end]
            t0_i = self.Ftest_4_offsets(data_i, t0[i])
            if t0_i is not None:
                if verbose:
                    print("    -> Significant. Goes to next test.")
                t0_tested.append(t0[i])
            else:
                if verbose:
                    print("    -> No nignificant. Goes to next test.")
        if len(t0_tested) > 0:
            self.offsets = list(np.unique(t0_tested))
        else:
            self.offsets = None
        self.offsets_values = None
        self.data = data_saved

    if self.offsets is not None:
        self.find_outliers_by_residuals()
        self.offsets = self.delete_small_offsets(self.offsets, del_by_pricise=True)

    self.outliers = ouliers_saved
