def interpolate(self, date='day', kind='linear', gap=10, in_place=False, verbose=False):
    """

    :param self: Gts instance
    :param date: 'day' will perform daily interpolation, alternatively date is a 1D numpy array with either datetime or decimal year
    :param method: scipy.interpolate.interp1d kind argument
    :param gap: maximum gap for interpolation
    :param in_place: boolean.
    :param verbose: verbose mode
    :return:
    """

    # import
    import pyacs.lib.astrotime as at
    import numpy as np
    import scipy.interpolate

    # init
    one_day_in_sec = 60 * 60 *24
    one_sec_dec_year = 1./ (365*24*60*60)
    # working gts

    wts = self.copy()
    wts.data_xyz = None

    # split time series according to gap
    lgts = wts.split_gap(gap,verbose=verbose)

    wts.data =  None

    for wwts in lgts:

        # handle dates
        np_seconds_data = at.decyear2seconds(wwts.data[:, 0])

        if date =='day':
            np_seconds = np.arange( np_seconds_data[0], np_seconds_data[-1], one_day_in_sec )

        else:

            lindex = np.argwhere( (date>wwts.date[0,0]- 10*one_sec_dec_year) & (date < wwts.date[0,0] + 10*one_sec_dec_year) )[0]
            wdate = date[lindex]

            if not isinstance(date, np.ndarray):
                print("-- ERROR: date argument must be day or a numpy array")
                import sys
                sys.exit()

            if date.dtype == 'float':
                # decyear case
                np_seconds = at.decyear2seconds( wdate )

            if np.issubdtype(date.dtype, np.datetime64):
                # decyear case
                np_seconds = at.datetime2seconds( wdate )

        # performs interpolation
        if np_seconds_data.shape[0] == 1:
            new_data = np.zeros((np_seconds.shape[0], 10))
            new_data[:, 0] = at.datetime2decyear(at.seconds2datetime(np_seconds))
            new_data[:, 1:4] = np.atleast_2d(wwts.data[:,1:4])
        else:
            f = scipy.interpolate.interp1d(np_seconds_data, wwts.data[:, 1:4],axis=0)
            yi = f(np_seconds)

            new_data = np.zeros((np_seconds.shape[0],10))
            new_data[:,0] = at.datetime2decyear( at.seconds2datetime( np_seconds ) )
            new_data[:,1:4] = yi

        # adds interpolation to the time series
        if wts.data is None:
            wts.data = new_data
        else:
            wts.data = np.vstack((wts.data,new_data))


    return wts