def insert_ts(self, ts, rounding='day',data='xyz',overlap=True):
    """

    :param ts: Gts to be inserted
    :param rounding: data rounding, used to decide whether an entry should be replaced.
        Choose among ['second','minute','hour','day]
    :param data: Gts attribute to be updated. 'xyz' for .data_xyz or None for .data
    :param overlap: if True, update occurs only on dates. If False, then ts overwrites the current Gts over the ts period
    :return: a new gts
    :note: The returned gts will have .data or .data_xyz will be set to None according to data argument
    """

    ###########################################################################
    # import
    ###########################################################################

    import numpy as np
    from pyacs.lib import astrotime as at
    from datetime import datetime
    import sys
    from pyacs.gts.Gts import Gts

    ###########################################################################
    def round_np_datetime(np_datetime, rounding):
    ###########################################################################

        from datetime import timedelta

        for i in np.arange(np_datetime.shape[0]):
            if rounding == 'day':
                np_datetime[i] = np_datetime[i].replace(hour=12, minute=0, second=0, microsecond=0)
            if rounding == 'hour':
                if np_datetime[i].minute >= 30:
                    np_datetime[i] = np_datetime[i] + timedelta(minutes=30)
                np_datetime[i] = np_datetime[i].replace(minute=0, second=0, microsecond=0)
            if rounding == 'minute':
                if np_datetime[i].second >= 30:
                    np_datetime[i] = np_datetime[i] + timedelta(seconds=30)
                np_datetime[i] = np_datetime[i].replace(second=0, microsecond=0)
            if rounding == 'second':
                np_datetime[i] = np_datetime[i].replace(microsecond=0)

        return np_datetime

    ###########################################################################

    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)

    # converts dates to seconds
    ###########################################################################

    if data is None:
        np_datetime_master_ts = round_np_datetime( at.decyear2datetime( self.data[:,0] ) , rounding )
        np_datetime_slave_ts = round_np_datetime( at.decyear2datetime( ts.data[:,0] ) , rounding )
    if data=='xyz':
        np_datetime_master_ts = round_np_datetime( at.decyear2datetime( self.data_xyz[:,0] ) , rounding )
        np_datetime_slave_ts = round_np_datetime( at.decyear2datetime( ts.data_xyz[:,0] ) , rounding )


    np_sec_master_ts = at.datetime2seconds( np_datetime_master_ts )
    np_sec_slave_ts  = at.datetime2seconds( np_datetime_slave_ts )

    # check that there is no duplicated dates
    ###########################################################################
    if not np.diff( np_sec_master_ts ).all():
        print("!!!ERROR: %s time series has duplicated dates after rounding" % self.code )
        sys.exit()

    if not np.diff( np_sec_slave_ts ).all():
        print("!!!ERROR: %s time series has duplicated dates after rounding" % ts.code )
        sys.exit()

    # check that dates are in increasing order
    ###########################################################################
    if not np.all( np.diff( np_sec_master_ts )>0 ):
        print("!!!ERROR: dates in %s time series are not in ascending order" % self.code )
        sys.exit()

    if not np.all( np.diff( np_sec_slave_ts )>0 ):
        print("!!!ERROR: dates in %s time series are not in ascending order" % ts.code )
        sys.exit()

    # if overlap is False, then removes all data corresponding to the period of ts
    ###########################################################################

    if not overlap:
        lindex = np.where( (np_sec_master_ts > np_sec_slave_ts[0]) &  (np_sec_master_ts < np_sec_slave_ts[-1]) )
        if data is None:
            master_data = np.delete(self.data,lindex,axis=0)
        if data=='xyz':
            master_data = np.delete(self.data_xyz,lindex,axis=0)

        np_sec_master_ts  = np.delete(np_sec_master_ts,lindex)

    else:
        if data is None:
            master_data = np.copy(self.data)
        if data=='xyz':
            master_data = np.copy(self.data_xyz)


    # removes dates in master that are in the inserted ts
    ###########################################################################

    values, xi, yi = np.array(np.intersect1d(np_sec_master_ts, np_sec_slave_ts, assume_unique=True, return_indices=True), dtype=int)
    master_data = np.delete(master_data,xi,axis=0)

    # concatenate arrays
    if data is None:
        slave_data = ts.data
    if data=='xyz':
        slave_data = ts.data_xyz

    new_data = np.vstack((master_data,slave_data))

    # sort array
    new_data = new_data[new_data[:, 0].argsort()]

    # creates the new gts
    ###########################################################################

    new_gts = self.copy()
    new_gts.data = None
    new_gts.data_xyz = None

    if data is None:
        new_gts.data = new_data
    if data=='xyz':
        new_gts.data_xyz = new_data

    return new_gts