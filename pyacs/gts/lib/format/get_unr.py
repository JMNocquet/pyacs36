
def get_unr( self , site , verbose=False ):

    """
    Get a time series from http://geodesy.unr.edu/gps_timeseries/txyz/IGS14/ in PYACS
    
    :param site: 4-letters code
    :param verbose: verbose mode
 
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
    url = ("http://geodesy.unr.edu/gps_timeseries/txyz/IGS14/%s.txyz2" % site.upper())
    
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
