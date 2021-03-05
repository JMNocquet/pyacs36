###############################################################################
def find_outliers_percentage(self, percentage=0.03, in_place=False, verbose=False, component='NEU', periods=None,
                             excluded_periods=None):
###############################################################################
    """
    detrend a time series and
    ranks the residuals by increasing absolute value
    populate the outliers with the x % largest ones on each component
    """

    # import
    import numpy as np


    n_to_remove = int(percentage * self.data.shape[0])
    if n_to_remove < 1: n_to_remove = 1

    if verbose: print("-- Removing the ", n_to_remove, " largest outliers (each time series)")

    new_ts = self.detrend()
    lindex_north = list(np.argsort(new_ts.data[:, 1] ** 2))[-n_to_remove:]
    lindex_east = list(np.argsort(new_ts.data[:, 2] ** 2))[-n_to_remove:]
    lindex_up = list(np.argsort(new_ts.data[:, 3] ** 2))[-n_to_remove:]

    lindex = []
    if 'N' in component: lindex = lindex + lindex_north
    if 'E' in component: lindex = lindex + lindex_east
    if 'U' in component: lindex = lindex + lindex_up

    if periods and excluded_periods:
        print("!!!! periods and excluded_periods provided. Possible overlap not checked.")
        print("!!!! The program will run first on periods and then will exclude outliers in excluded_periods.")
    if periods:
        lkept_index = []
        for index in lindex:
            for period in periods:
                start_date_period = period[0]
                end_date_period = period[1]
                if self.data[index, 0] >= start_date_period and self.data[index, 0] <= end_date_period:
                    lkept_index.append(index)
                    break
        lindex = lkept_index

    if excluded_periods:
        lexcluded_index = []
        for index in lindex:
            for period in excluded_periods:
                start_date_period = period[0]
                end_date_period = period[1]
                if self.data[index, 0] >= start_date_period and self.data[index, 0] <= end_date_period:
                    lexcluded_index.append(index)
                    break
        lkept_index = []
        for index in lindex:
            if index not in lexcluded_index: lkept_index.append()
        lindex = lkept_index

    new_Gts = self.copy()

    if in_place:
        self.outliers = list(set(lindex))
        return (self)
        del new_Gts
    else:
        new_Gts = self.copy()
        new_Gts.outliers = list(set(lindex))
        return (new_Gts)
