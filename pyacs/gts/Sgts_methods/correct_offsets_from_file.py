def correct_offsets_from_file(self, offset_file, fill_offsets_dates=False):
    """Read offsets from a text file and apply them to the time series.

    Parameters
    ----------
    offset_file : str
        Path to the offset file.
    fill_offsets_dates : bool, optional
        If True, fill the offsets_dates attribute of each Gts. Default False.

    Returns
    -------
    Sgts
        Copy of the current Sgts with offsets applied (one Gts per site that
        appears in the file and in the collection).

    Notes
    -----
    Offset file is whitespace-separated; lines starting with ``#`` are
    ignored. Expected columns: lon, lat, de, dn, du, sde, sdn, sdu, code,
    offset_date, offset_info. Displacements (de, dn, du, sde, sdn, sdu) are
    in mm in the file. Example::

        #      lon        lat       de       dn       du      sde      sdn      sdu    code       offset_date        offset_info
          -79.1000     0.1000  -103.51     3.67     0.00     2.37     2.60     0.00    ARSH  2016-04-16T23:58  Mw7.8/-80.25/-0.12/22
    """
    # import
    import pandas as pd
    import pyacs.lib.astrotime as at
    import numpy as np
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.debug

    # copy the time series
    ts = self.copy()

    # reads offset file
    # for compatibility with pandas new release - JMN 06/11/2025
    data = pd.read_csv(offset_file, sep='\s+',
                       names=['lon', 'lat', 'de', 'dn', 'du', 'sde', 'sdn', 'sdu', 'code', 'offset_date',
                              'offset_info'], comment='#')
    # old 
    #data = pd.read_csv(offset_file, delim_whitespace=True,
    #                   names=['lon', 'lat', 'de', 'dn', 'du', 'sde', 'sdn', 'sdu', 'code', 'offset_date',
    #                          'offset_info'], comment='#')
    # convert to meters
    data['de'] = data['de'] / 1000.
    data['dn'] = data['dn'] / 1000.
    data['du'] = data['du'] / 1000.
    data['sde'] = data['sde'] / 1000.
    data['sdn'] = data['sdn'] / 1000.
    data['sdu'] = data['sdu'] / 1000.

    # convert offset_info to datetime
    data['date'] = pd.to_datetime(data['offset_date'], format='%Y-%m-%dT%H:%M')

    # loop over the code
    for code in data['code'].unique():
        if not ts.has_ts(code):
            continue
        # get the offset
        offset = data[data['code'] == code]
        # get dates for code
        np_datetime = at.decyear2datetime(ts.__dict__[code].data[:, 0])
        # loop over the time series
        for index, row in offset.iterrows():
            lidx = np.where(np_datetime > row['date'])[0]
            if len(lidx) == 0:
                continue
            else:
                VERBOSE('Correcting offset for code %s at date %s %8.2lf %8.2lf %8.2lf'
                        % (code, row['date'].isoformat()[:16], row['de']*1E3, row['dn']*1E3, row['du']*1E3))
                ts.__dict__[code].data[lidx, 1] = ts.__dict__[code].data[lidx, 1] - row['dn']
                ts.__dict__[code].data[lidx, 2] = ts.__dict__[code].data[lidx, 2] - row['de']
                ts.__dict__[code].data[lidx, 3] = ts.__dict__[code].data[lidx, 3] - row['du']
                if fill_offsets_dates:
                    try:
                        ts.__dict__[code].offsets_dates.append(at.datetime2decyear(row['date']))
                    except:
                        ts.__dict__[code].offsets_dates = [at.datetime2decyear(row['date'])]
    return ts