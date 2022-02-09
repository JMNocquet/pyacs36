def split_gap(self,gap=10,verbose=False):
    """

    :param gap: gap in number of days to split the time series
    :param verbose: verbose mode
    :return: a list a gts split from the original
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime as at


    # convert time to seconds
    dd = at.decyear2seconds(self.data[:, 0])

    # find whether the time series has data gaps
    dd_diff = np.diff(dd) / (24 * 60 * 60)
    gaps = np.argwhere(dd_diff > gap).flatten()
    gaps = np.sort(np.append(gaps, np.append(gaps + 1, [0, dd.shape[0] - 1])))

    periods = []
    for i in np.arange(int(gaps.shape[0] / 2)):
        periods.append([dd[gaps[2 * i]], dd[gaps[2 * i + 1]]])

    np_outliers = None

    # return lGts
    R_gts = []

    for period in periods:
        datetime1 = at.seconds2datetime(period[0]-1)
        datetime2 = at.seconds2datetime(period[1]+1)
        decyr1 = at.datetime2decyear(datetime1)
        decyr2 = at.datetime2decyear(datetime2)
        if verbose:
            print("---- period %s -- %s " % (datetime1.isoformat(), \
                                             datetime2.isoformat()))
        period_decyear=[at.datetime2decyear(datetime1),at.datetime2decyear(datetime2)]
        R_gts.append( self.extract_periods(period_decyear,no_reset=True) )

    if verbose:
        print("---- returning %d gts time series" % len(R_gts))

    return R_gts