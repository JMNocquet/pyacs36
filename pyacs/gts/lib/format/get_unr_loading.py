def get_unr_loading(self, site, verbose=False):
    """
    Get NTAL, NTOL, HYDL, MASC, RESE loading predictions from UNR.

    See http://geodesy.unr.edu/gps_timeseries/README_tenv3load.txt for format.

    Parameters
    ----------
    site : str
        Station 4-letter code.
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    tuple of 5 Gts
        (NTAL, NTOL, HYDL, MASC, RESE) loading time series.

    Examples
    --------
    NTAL_CSEC, NTOL_CSEC, HYDL_CSEC, MASC_CSEC, RESE_CSEC = Sgts().get_unr_loading('CSEC')

    Notes
    -----
    UNR .tenv3 files include load prediction triplets (NTAL, NTOL, HYDL, MASC, RESE, WUS, WUSA).
    Columns 21-23: NTAL; 24-26: NTOL; 27-29: HYDL; 30-32: MASC; 33-35: RESE.
    Prior to 2002 GRACE mascon data are unavailable (extrapolation used).
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
    url = ("http://geodesy.unr.edu/gps_timeseries/tenv3_loadpredictions/%s.tenv3" % site.upper())

    # get data

    try:
        urllib.request.urlretrieve(url=url, filename="test.dat")
    except HTTPError as error:
        print('Data not retrieved because %s\nURL: %s', error, url)
    except URLError as error:
        if isinstance(error.reason, socket.timeout):
            print('socket timed out - URL %s', url)
        else:
            print('some other error happened')

    # get code
    code = np.genfromtxt('test.dat', usecols=0, dtype=str)[0]

    # get dates
    np_mjd = np.genfromtxt('test.dat', usecols=(3), dtype=int, skip_header=1)
    np_datetime = at.mjd2datetime(np_mjd) + delta_12h
    # instanciate Gts
    gts_ntal = Gts(code=site+'_ntal')
    gts_ntol = Gts(code=site+'_ntol')
    gts_hydl = Gts(code=site+'_hydl')
    gts_masc = Gts(code=site+'_masc')
    gts_rese = Gts(code=site+'_rese')

    # get data
    # col 21-23: east, north and up component of non-tidal atmospheric loading displacements  (NTAL)
    # col 24-26: east, north and up components of non-tidal ocean loading displacements  (NTOL)
    # col 27-29: east, north and up components of hydrological model load displacements (HYDL)
    # col 30-32: east, north, and up components from GRACE mascon-based loading displacements (MASC)
    # col 33-35: east, north, and up components from of reservoir induced loading displacements (RESE)
    # Get the NTAL, NTOL, HYDL, MASC, RESE loading predictions from UNR

    gts_ntal.data = np.zeros((np_datetime.shape[0],10))
    gts_ntal.data[:,1:4] = np.genfromtxt('test.dat', usecols=(20, 21, 22), skip_header=1)
    gts_ntal.data[:,0] = at.datetime2decyear( np_datetime )

    gts_ntol.data = np.zeros((np_datetime.shape[0],10))
    gts_ntol.data[:,1:4] = np.genfromtxt('test.dat', usecols=(23, 24, 25), skip_header=1)
    gts_ntol.data[:,0] = at.datetime2decyear( np_datetime )

    gts_hydl.data = np.zeros((np_datetime.shape[0],10))
    gts_hydl.data[:,1:4] = np.genfromtxt('test.dat', usecols=(26, 27, 28), skip_header=1)
    gts_hydl.data[:,0] = at.datetime2decyear( np_datetime )

    gts_masc.data = np.zeros((np_datetime.shape[0],10))
    gts_masc.data[:,1:4] = np.genfromtxt('test.dat', usecols=(29, 30, 31), skip_header=1)
    gts_masc.data[:,0] = at.datetime2decyear( np_datetime )

    gts_rese.data = np.zeros((np_datetime.shape[0],10))
    gts_rese.data[:,1:4] = np.genfromtxt('test.dat', usecols=(32, 33, 34), skip_header=1)
    gts_rese.data[:,0] = at.datetime2decyear( np_datetime )

    # remove 'test.dat'
    os.remove('test.dat')

    # return
    return (gts_ntal, gts_ntol, gts_hydl, gts_masc, gts_rese)
