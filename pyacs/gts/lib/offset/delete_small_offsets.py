"""
Delete small offsets (below threshold) from a list of offset dates.
"""

import numpy as np


def delete_small_offsets(self, offsets, del_by_pricise=False):
    """Estimate offsets with clean data, then remove offsets that are too small.

    Parameters
    ----------
    offsets : list
        List of offset dates (decimal year).
    del_by_pricise : bool, optional
        If True, also require offset > sigma. Default is False.

    Returns
    -------
    list or None
        Filtered list of offset dates, or None.
    """
    min_offset_NE = 0.002
    #if self.H_conf and 'threshold_smallest_offset' in self.H_conf:
    #    min_offset_NE = self.H_conf['threshold_smallest_offset']
    min_offset_UP = 3. * min_offset_NE
    self.offsets = offsets
    new_Gts = self.make_model(loutlier=self.outliers)
    if (new_Gts is not None) and (new_Gts.offsets_values is not None):
        offsets_values = np.fabs(new_Gts.offsets_values)
        choose_index = []
        for i in range(len(self.offsets)):
            if del_by_pricise:
                if (offsets_values[i, 1] > min_offset_NE) and (offsets_values[i, 1] > offsets_values[i, 4]):
                    choose_index.append(i)
                if (offsets_values[i, 2] > min_offset_NE) and (offsets_values[i, 2] > offsets_values[i, 5]):
                    choose_index.append(i)
                if (offsets_values[i, 3] > min_offset_UP) and (offsets_values[i, 3] > offsets_values[i, 6]):
                    choose_index.append(i)
            else:
                if offsets_values[i, 1] > min_offset_NE:
                    choose_index.append(i)
                if offsets_values[i, 2] > min_offset_NE:
                    choose_index.append(i)
                if offsets_values[i, 3] > min_offset_UP:
                    choose_index.append(i)
        if len(choose_index) > 0:
            choose_index = np.unique(choose_index)
            offsets = list(np.take(self.offsets, choose_index, axis=0))
        else:
            offsets = None
    self.offsets = offsets
    self.offsets_values = None
    return offsets
