def find_l1trend(self,
                           lam,
                           threshold,
                           period=None,
                           gap=10,
                           components= 'NE' ,
                           plot=False,
                           verbose=False,
                           in_place=False):
    """

    :param self: Gts instance
    :param lam: lambda parameter for L1 trend filtering
    :param threshold: All residuals with threshold * standard deviation will be flagged as outliers
    :param period: period(s) for searching outliers. Could be a single of a list of periods.
    :param gap: number of days to consider that there is a gap. Default is gap=10.
    :param components: components used for outliers detection
    :param plot: boolean. If True, will plot the filter result and the flagged outliers
    :param verbose: boolean. Verbose mode.
    :param in_place: boolean. if True, apply to the original Gts. Default is False, returning a new Gts
    :return: a new Gts instance if in_place is False or the current Gts
    """

    # import
    import pyacs.lib.astrotime as at
    import numpy as np

    # initialize the working time series

    if period is None:
        ts = self.copy()
    else:
        ts = self.extract_periods(period)
    ts.outliers = []

    dd = at.decyear2seconds(ts.data[:, 0])

    # identify outliers dates
    #if dd.shape[0] > 4:
    l1 = ts.l1_trend(lam, gap=gap,component=components, verbose=verbose)
    #else:
    #    l1 = ts.detrend(method='L1')

    # residuals and mean
    r = (ts.data[:, 1:4] - l1.data[:, 1:4]) * 1.E3
    mr = np.mean(np.fabs(r), axis=0)
    if verbose:
        print("-- median of |l1trend - obs| NEU: %.1lf %.1lf %.1lf threshold %.1lf" % (mr[0], mr[1], mr[2], threshold))

    lindexh = None
    lindexu = None

    # horizontal
    if ('N' in components) and ('E' in components):
        threshold = np.sqrt(mr[0] ** 2 + mr[1] ** 2) * threshold
        lindexh = np.argwhere((np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2) > threshold)).flatten()

    else:
        if ('N' in components):
            threshold = mr[0] * threshold
            lindexh = np.argwhere((np.sqrt(r[:, 0] ** 2) > threshold)).flatten()

        if ('E' in components):
            threshold = mr[1] * threshold
            if verbose:
                print(
                    "-- median to filter for E: %.1lf %.1lf %.1lf threshold %.1lf" % (mr[0], mr[1], mr[2], threshold))
            lindexh = np.argwhere((np.sqrt(r[:, 1] ** 2) > threshold)).flatten()

    if ('U' in components):
        threshold = mr[2] * threshold
        lindexu = np.argwhere((np.sqrt(r[:, 0] ** 2) > threshold)).flatten()

    if lindexh is not None:
        lindex = lindexh.tolist()
    else:
        lindex = []

    if lindexu is not None:
        lindex = lindex + lindexu

    outlier_dates = np.array(dd[lindex], dtype=int)

    if verbose:
        print("-- found %d outliers" % (outlier_dates.shape[0]))

    if plot:
        ts.plot(superimposed=[l1], date_unit='cal', center=False, plot_size=(10, 15), title=("%s data vs l1_trend" % ts.code ))


    ts.outliers = np.searchsorted(at.decyear2seconds(ts.data[:, 0]), outlier_dates).tolist()

    if verbose:
        print("-- Number of outliers: %d" % len(ts.outliers))

    if in_place:
        self.outliers = ts.outliers
        return(self)
        del ts
    else:
        return ts

