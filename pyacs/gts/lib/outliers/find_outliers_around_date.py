###################################################################
def find_outliers_around_date(self, date, conf_level=95, n=3, lcomponent='NE', verbose=True):
###################################################################
    """
    Find an outlier around a given date using F-ratio test.

    Returns the index of the outlier, or [] if none found.

    Parameters
    ----------
    date : float
        Date in decimal year.
    conf_level : float, optional
        Confidence level for F-ratio test (default 95).
    n : int, optional
        Number of dates on each side of date (default 3).
    lcomponent : str, optional
        Components to test: 'N', 'E', 'U', 'NE', 'NEU' (default 'NE').
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    Gts or list
        self with outlier flagged, or [] if no significant outlier.
    """

    # import
    import numpy as np
    from pyacs.gts.Gts import get_index_from_dates

    if verbose:
        print((
                          "-- Searching outlier around date %10.5lf on components %s with confidence level %6.1lf and %02d samples" % (
                  date, lcomponent, conf_level, 2 * n)))

    # self.info()
    tmp_gts = self.detrend().remove_outliers().extract_ndates_around_date(date, n)
    nn = tmp_gts.data.shape[0]

    score = {}
    llindex = {}

    # F_ratio test
    ###############################################
    def f_ratio(chi_square_1, p1, chi_square_2, p2, n):
        ###############################################
        """
        Return the F-ratio test CDF value (probability).

        Parameters
        ----------
        chi_square_1, chi_square_2 : float
            Chi-square values.
        p1, p2 : int
            Number of parameters.
        n : int
            Sample size.

        Returns
        -------
        float
            F CDF value.
        """
        F = ((chi_square_1 - chi_square_2) / (p2 - p1)) / (chi_square_2 / (n - p2))

        from scipy.stats import f
        return (f.cdf(F, p2 - p1, n - p2))

    #
    H_component = {1: 'North', 2: 'East', 3: 'Up'}
    # find outlier
    li = []
    if 'N' in lcomponent: li.append(1)
    if 'E' in lcomponent: li.append(2)
    if 'U' in lcomponent: li.append(3)

    for i in sorted(li):
        #        if verbose:
        #            if i==1:print "  => Testing component: North"
        #            if i==2:print "  => Testing component: East"
        #            if i==3:print "  => Testing component: Up"

        index = np.where(np.abs(tmp_gts.data[:, i] - np.median(tmp_gts.data[:, i])) == np.max(
            np.abs(tmp_gts.data[:, i] - np.median(tmp_gts.data[:, i]))))

        if verbose:
            print(
                ("-- suspected outlier at date %10.4lf on component %s " % (tmp_gts.data[index, 0][0], H_component[i])))

        tmp_gts_no_outlier = tmp_gts.copy()
        tmp_gts_no_outlier.outliers = [index]
        tmp_gts_no_outlier.remove_outliers(in_place=True)

        chi_square_1 = nn * np.std(tmp_gts.data[:, i]) ** 2
        chi_square_2 = (nn - 1) * np.std(tmp_gts_no_outlier.data[:, i]) ** 2

        score[i] = f_ratio(chi_square_1, 1, chi_square_2, 2, 2 * n) * 100.0
        print(("-- probability of outlier (F_ratio) %5.2lf%% " % (score[i])))
        llindex[i] = index

    # make decision

    if np.max(list(score.values())) < conf_level:
        if verbose: print("-- No significant outlier found")
        return (self)
    else:
        # choose the outlier as the maximum probability
        component = li[0]
        current_score = score[component]
        del li[0]
        for i in li:
            if score[i] > current_score:
                current_score = score[i]
                component = i

        date = tmp_gts.data[llindex[component], 0]
        # return the index in the original time series
        if verbose:
            print("=> Getting index for date ", date)
        returned_index = get_index_from_dates([date], self.data, tol=0.25)
        self.outliers += returned_index

        return (self)
