
def get_unr( self , site , verbose=False ):

    """Fetch a time series from UNR (IGS20 txyz) and return as Gts.

    Parameters
    ----------
    site : str
        4-letter site code.
    verbose : bool, optional
        If True, print progress. Default is False.

    Returns
    -------
    Gts or None
        Loaded time series or None on failure.
    """

    # import    
    import urllib.request
    from urllib.error import HTTPError, URLError
    import socket
    from pyacs.gts.Gts import Gts
    import numpy as np
    import os
    import pyacs.lib.astrotime as at
    import datetime
    from datetime import timedelta

    delta_12h = timedelta(hours=12)

    # url
    url = ("https://geodesy.unr.edu/gps_timeseries/IGS20/txyz/%s.txyz2" % site.upper())
    # template on 11/09/2024
    #https://geodesy.unr.edu/gps_timeseries/IGS20/txyz/BCES.txyz2
    # get data
    
    try:
        urllib.request.urlretrieve(url = url, filename = "test.dat")
    except HTTPError as error:
        print('Data not retrieved because %s\nURL: %s', error, url)
    except URLError as error:
        if isinstance(error.reason, socket.timeout):
            print('socket timed out - URL %s', url)
        else:
            print('some other error happened')

    # creates Gts
    
    # get code
    code = np.genfromtxt('test.dat' , usecols=0, dtype=str )[0]
    gts = Gts( code = code )
    
    # get data
    gts.data_xyz = np.genfromtxt('test.dat' , usecols=(2,3,4,5,6,7,8,9,10,11) )

    # decimal year dates in UNR files only have 4 digits, making the day time very approximate
    # we round the dates at 12:00 and prefer the date string
    str_date =  np.genfromtxt('test.dat' , usecols=(1) , dtype=str )
    np_datetime = np.array([datetime.datetime.strptime(x, "%y%b%d") for x in str_date]) + delta_12h

    #unr_dates_decyear = gts.data_xyz[:,0]
    #np_year = np.array(unr_dates_decyear, dtype=int)
    #(np_doy,_np_ut) = at.decyear2dayno( unr_dates_decyear )
    #gts.data_xyz[:,0] = at.dayno2decyear( np_doy , np_year )
    gts.data_xyz[:, 0] = at.datetime2decyear( np_datetime )
    # convert data
    gts.xyz2neu(corr=True,ref_xyz=None, verbose=verbose)
    
    # remove 'test.dat'
    
    os.remove('test.dat')
    
    # return
    return gts
