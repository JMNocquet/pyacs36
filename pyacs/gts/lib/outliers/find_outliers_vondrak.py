###############################################################################
def find_outliers_vondrak(self, threshold=10, fc=2., in_place=False, verbose=True, \
                          periods=[[]], excluded_periods=[[]], component='NE'):
###############################################################################
    """
    Find outliers using a Vondrak filter
    """
    # import
    import numpy as np
    from pyacs.gts.Gts import get_index_from_dates

    # init
    loutliers_dates = []
    lindex_north = []
    lindex_east = []
    lindex_up = []

    # keep selected period
    tmp_ts = self.extract_periods(periods).exclude_periods(excluded_periods).detrend()

    # calculates vondrak filter
    vondrak_ts = tmp_ts.vondrak(fc, component=component, verbose=verbose)

    # calculates the residual time series

    residual_ts = tmp_ts.copy()

    residual_ts.data[:, 1] = tmp_ts.data[:, 1] - vondrak_ts.data[:, 1]
    residual_ts.data[:, 2] = tmp_ts.data[:, 2] - vondrak_ts.data[:, 2]
    residual_ts.data[:, 3] = tmp_ts.data[:, 3] - vondrak_ts.data[:, 3]

    # get the median
    [median_north, median_east, median_up] = np.median(np.abs(residual_ts.data[:, 1:4]), axis=0)

    # get the outliers
    if 'N' in component:
        lindex_north = np.where(np.abs(residual_ts.data[:, 1]) > threshold * median_north)[0].tolist()
    if 'E' in component:
        lindex_east = np.where(np.abs(residual_ts.data[:, 2]) > threshold * median_east)[0].tolist()
    if 'U' in component:
        lindex_up = np.where(np.abs(residual_ts.data[:, 3]) > threshold * median_up)[0].tolist()

    loutliers = list(set(lindex_north + lindex_east + lindex_up))

    if verbose: print(("-- Outliers detection using Vondrak filter with fc=%.2lf : %03d new outliers detected" % (
    fc, len(loutliers))))

    # get the outliers dates
    loutliers_dates += tmp_ts.data[loutliers, 0].tolist()

    # get outliers index in original time series
    if loutliers != []:
        loutliers_index = get_index_from_dates(loutliers_dates, self.data, tol=0.25)


    else:
        loutliers_index = self.outliers

    # return
    new_Gts = self.copy()

    if in_place:
        self.outliers = loutliers_index
        return (self)
        del new_Gts
    else:
        new_Gts = self.copy()
        new_Gts.outliers = loutliers_index
        return (new_Gts)
