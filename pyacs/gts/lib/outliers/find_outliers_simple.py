def find_outliers_simple(self, threshold=100, window_length=10, in_place=False, verbose=False, component='NEU',
                         periods=None, excluded_periods=None):

    import numpy as np

    lindex_outlier = []
    for i in range(0, self.data.shape[0] - window_length):
        window = self.data[i:i + window_length, 1:4] * 1000.
        # print 'window',window
        # print 'median ',np.median(window,axis=0)
        residuals = np.abs(window - np.median(window, axis=0))
        # print 'residuals',residuals
        median_of_residuals = np.median(residuals, axis=0)
        # print 'median of res ',median_of_residuals
        for index in range(0, residuals.shape[0]):
            # print residuals[index,0],median_of_residuals[0]

            if residuals[index, 0] > threshold * median_of_residuals[0]:
                if (i + index) not in lindex_outlier: print('outlier at ', self.data[i + index, 0],
                                                            ' N');lindex_outlier.append(i + index)
            if residuals[index, 1] > threshold * median_of_residuals[1]:
                if (i + index) not in lindex_outlier: print('outlier at ', self.data[i + index, 0],
                                                            ' E');lindex_outlier.append(i + index)
            if residuals[index, 2] > threshold * median_of_residuals[2]:
                if (i + index) not in lindex_outlier: print('outlier at ', self.data[i + index, 0],
                                                            ' U');lindex_outlier.append(i + index)

    if periods and excluded_periods:
        print("!!!! periods and excluded_periods provided. Possible overlap not checked.")
        print("!!!! The program will run first on periods and then will exclude outliers in excluded_periods.")
    if periods:
        lkept_index = []
        for index in lindex_outlier:
            for period in periods:
                start_date_period = period[0]
                end_date_period = period[1]
                if self.data[index, 0] >= start_date_period and self.data[index, 0] <= end_date_period:
                    lkept_index.append(index)
                    break
        lindex_outlier = lkept_index

    if excluded_periods:
        lexcluded_index = []
        for index in lindex_outlier:
            for period in excluded_periods:
                start_date_period = period[0]
                end_date_period = period[1]
                if self.data[index, 0] >= start_date_period and self.data[index, 0] <= end_date_period:
                    lexcluded_index.append(index)
                    break
        lkept_index = []
        for index in lindex_outlier:
            if index not in lexcluded_index: lkept_index.append()
        lindex_outlier = lkept_index

    new_Gts = self.copy()

    if in_place:
        self.outliers = list(set(lindex_outlier))
        return (self)
        del new_Gts
    else:
        new_Gts = self.copy()
        new_Gts.outliers = list(set(lindex_outlier))
        return (new_Gts)
