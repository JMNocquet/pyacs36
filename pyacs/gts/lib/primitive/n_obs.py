def n_obs(self,date):
    """
    Return the number of available observations for the given date or period.

    Parameters
    ----------
    date : str or list
        Date or period. Formats: decimal year or pandas-style (e.g. '2018-01-01',
        ['2018-01-01','2018-02-01'], or [2019., 2019.5] for decimal year).

    Returns
    -------
    int
        Number of observations in the current time series for the given date/period.

    Examples
    --------
    A specific day: ts.RIOP.n_obs('2018-01-01')
    A period: ts.RIOP.n_obs(['2018-01-01','2018-02-01'])
    A period in decimal year: ts.RIOP.n_obs([2019., 2019.5])
    """
    # import
    import pyacs.lib.astrotime as at
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG


    print(date)
    # init n
    n = None

    # convert to pandas df
    df = self.to_pandas_df()

    # decipher input
    if isinstance(date,str):
        try:
            n = len( df.loc[date] )
        except:
            n = 0
            pass

    if isinstance(date,list):
        # ensures that it is a list of list
        if not isinstance(date[0], list):
            date = [date]

        # loop on periods
        n = 0
        for period in date:
            if isinstance(period[0],str):
                sd = period[0]
            else:
                print(at.seconds2datetime( at.decyear2seconds(period[0]) ))
                sd = at.seconds2datetime( at.decyear2seconds(period[0]) ).isoformat()
            if isinstance(period[1],str):
                ed = period[1]
            else:
                ed = at.seconds2datetime( at.decyear2seconds(period[1]) ).isoformat()
            n = n + len(df.loc[sd:ed])

    if n is None:
        WARNING("Something went wrong during Gts.n_obs. Argument passed was %s" % str(date))

    return n