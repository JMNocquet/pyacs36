###############################################################################
def find_outliers_sliding_window(self, \
                                 threshold=3, in_place=False, verbose=True, \
                                 periods=[[]], excluded_periods=[[]], component='NE', window_len=15, automatic=True):
    ###############################################################################
    """
    Find outliers using sliding windows
    """

    import numpy as np
    from pyacs.gts.Gts import get_index_from_dates


    lindex_north = []
    lindex_east = []
    lindex_up = []

    if self.data.shape[0] > window_len:

        itermax = 5

        lindex_north = []
        lindex_east = []
        lindex_up = []

        OK = True
        loutliers = []
        loutliers_dates = []
        i = 0

        smooth = self.extract_periods(periods).exclude_periods(excluded_periods).smooth(window_len=window_len)
        new_ts = self.extract_periods(periods).exclude_periods(excluded_periods)
        residual_ts = self.extract_periods(periods).exclude_periods(excluded_periods)
        residual_ts.data[:, 1:4] = new_ts.data[:, 1:4] - smooth.data[:, 1:4]

        diff_data = np.diff(self.data[:, 1:4], n=1, axis=0)
        [median_north, median_east, median_up] = np.median(np.abs(diff_data), axis=0)

        while OK:

            if 'N' in component:
                lindex_north = np.where(np.abs(residual_ts.data[:, 1]) > threshold * median_north)[0].tolist()
            if 'E' in component:
                lindex_east = np.where(np.abs(residual_ts.data[:, 2]) > threshold * median_east)[0].tolist()
            if 'U' in component:
                lindex_up = np.where(np.abs(residual_ts.data[:, 3]) > threshold * median_up)[0].tolist()

            loutliers = list(set(lindex_north + lindex_east + lindex_up))

            if verbose: print(("-- Outliers detection pass #%02d : %03d new outliers detected" % (i, len(loutliers))))

            # print loutliers_dates,new_ts.data[loutliers,0].tolist()
            loutliers_dates += new_ts.data[loutliers, 0].tolist()

            if loutliers == []: OK = False

            i += 1
            if i > itermax: OK = False

            smooth = self.extract_periods(periods).exclude_periods([[]]).smooth(window_len=window_len)

            new_ts.outliers = loutliers
            new_ts = new_ts.remove_outliers()

            smooth = new_ts.smooth(window_len=window_len)

            residual_ts = new_ts.copy()
            residual_ts.data[:, 1:4] = new_ts.data[:, 1:4] - smooth.data[:, 1:4]

            diff_data = np.diff(self.data[:, 1:4], n=1, axis=0)
            [median_north, median_east, median_up] = np.median(np.abs(diff_data), axis=0)

        if verbose: print("-- ", len(loutliers_dates), " outliers found")

        loutliers_index = get_index_from_dates(loutliers_dates, self.data, tol=0.25)


    else:
        loutliers_index = self.outliers

    new_Gts = self.copy()

    if in_place:
        self.outliers = loutliers_index
        return (self)
        del new_Gts
    else:
        new_Gts = self.copy()
        new_Gts.outliers = loutliers_index
        return (new_Gts)

