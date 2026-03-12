def disp2vel(self, alpha='auto',ndays_mf=1):
    """Compute time series derivative via L1 trend filter (BIC-optimal parameter).

    Uses https://pypi.org/project/trendfilter/

    Parameters
    ----------
    alpha : str or float, optional
        Filter parameter; 'auto' for BIC-optimal. Default is 'auto'.
    ndays_mf : int, optional
        Median filter window (days) before L1-trend. Default is 1.

    Returns
    -------
    Gts
        New Gts instance; units m/yr.
    """


    fts = self.l1trend(alpha=alpha,ndays_mf=ndays_mf)

    new_gts = fts.ivel()

    return new_gts