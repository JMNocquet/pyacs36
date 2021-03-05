###############################################################################
def remove_outliers(self, periods=None, in_place=False):
###############################################################################
    """
    removes outliers provided in self.outliers
    return a new Gts without the outliers
    if in_place = True then self has the outliers removed as well (in _place)
    """

    # import
    import numpy as np

    if self.outliers:
        if periods == None:
            data = np.delete(self.data, self.outliers, axis=0)
        else:
            lindex = self.lindex_in_periods(periods)
            ldelete = np.intersect1d(lindex, self.outliers)
            data = np.delete(self.data, ldelete, axis=0)


    else:
        data = np.copy(self.data)

    new_Gts = self.copy()
    new_Gts.outliers = []
    new_Gts.data = data

    if in_place:
        self.data = new_Gts.data.copy()
        del new_Gts
        self.outliers = []

        return (self)
    else:
        return (new_Gts)
