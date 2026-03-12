def get_values_at_date(self, date, data_type='data_xyz'):
    """
    return the values of the time series at a given date.
    date can be provided as a datetime, string, decimal year, modified julian day, or (doy,year)
    data_type can be 'data' or 'data_xyz'

    parameters:
    -----------
    date: datetime, string, float, tuple
        the date to look for
    data_type: string
        'data' or 'data_xyz' depending on the type of data to return
    """

    # import
    from datetime import datetime
    import pyacs.lib.astrotime as at
    import numpy as np

    # converts to seconds
    dates_ts = at.decyear2seconds(self.data[:, 0])

    # handle different date formats
    if type(date) == datetime:
        date = at.datetime2seconds(date)
    if type(date) == str:
        date = at.datetime2seconds(datetime.strptime(date, '%Y-%m-%d'))
    if type(date) == float:
        if date < 10000:
            date = at.decyear2seconds(date)
        else:
            date = at.mjd2seconds(date)
    if type(date) == tuple:
        date = at.datetime2seconds(at.dayno2datetime(date[0], date[1]))

    # find the closest date
    idx = np.argmin(np.abs(dates_ts - date))

    if data_type == 'data':
        return self.data[idx, :]
    else:
        return self.data_xyz[idx, :]




