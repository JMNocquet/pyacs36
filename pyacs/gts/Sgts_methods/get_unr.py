def get_unr(self,lcode=[], center=None,radius=None):
    """Download time series from UNR and append to this Sgts.

    Parameters
    ----------
    lcode : list, optional
        Site codes to download. Default is [].
    center : list, optional
        [lon, lat] in decimal degrees for center of search.
    radius : list, optional
        [min_radius, max_radius] in km.

    Returns
    -------
    None
        Time series are appended to self (in-place).
    """

    # import
    from pyacs.gts.Gts import Gts
    from pyacs.gts.Sgts import Sgts

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import urllib.request
    from bs4 import BeautifulSoup
    import numpy as np
    import pandas as pd
    import io
    import pyacs.lib.coordinates as coo
    import pyacs.lib.astrotime as at
    import os

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # initialize Sgts
    ts = Sgts(read=False)

    # case center and radius
    if center is not None and radius is not None:

        url = "http://geodesy.unr.edu/NGLStationPages/gpsnetmap/GPSNetMap.html"

        request = urllib.request.Request(
            url=url,
            # set a real User-Agent header
            # to avoid 403 errors
            headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36'
            }
        )
        # perform the HTTP GET request and automatically freeing
        # the resources allocated
        with urllib.request.urlopen(request) as response:
            # get the HTML content of the target web page
            body = response.read()

        html = body.decode('utf-8')

        # decipher html and read as panda frame
        soup = BeautifulSoup(html,features="lxml")
        LIST_STRING = soup.script.prettify().split(';')[1].split('=')[-1].split("\n")
        LIST_STRING_CLEAN = [line.replace('"', "").replace("[", "").replace("]", "").strip(',').replace(",", " ") for
                             line in LIST_STRING]

        DF = pd.read_csv(io.StringIO('\n'.join(LIST_STRING_CLEAN)), sep='\s+', header=None,
                         names=["code", "latitude", "longitude", "dummy1", "dummy2"],
                         usecols=["code", "latitude", "longitude"])

        # convert to numpy
        np_code = DF['code'].to_numpy(dtype=str)
        np_lon = DF['longitude'].to_numpy()
        np_lat = DF['latitude'].to_numpy()

        lidx = np.where(np_lon < -180)[0]
        np_lon[lidx] = np_lon[lidx] + 360.

        lidx = np.where(np_lon > 180)[0]
        np_lon[lidx] = np_lon[lidx] - 360.

        # save a dictionary of code with coordinates
        H_code = {}
        for i in np.arange(np_code.shape[0]):
            H_code[np_code[i]] = [np_lon[i], np_lat[i]]

        # convert EQ lon, lat into X,Y,Z coordinates
        (X_CENTER, Y_CENTER, Z_CENTER) = coo.geo2xyz(center[0], center[1], 0., unit='dec_deg')

        # convert GPS lon, lat into X,Y,Z coordinates
        (X, Y, Z) = coo.geo2xyz(np_lon, np_lat, 0, unit='dec_deg')

        # get list of
        D = np.sqrt((X - X_CENTER) ** 2 + (Y - Y_CENTER) ** 2 + (Z - Z_CENTER) ** 2) * 1.E-3  # to get results in km
        lidx = np.where( (D <= radius[1]) & (D >= radius[0]) )
        lcode = lcode + np_code[lidx].tolist()

    if lcode == []:
        WARNING("No site to download. Returning empty Sgts")
        return ts

    # loop on sites

    for code in lcode:
        VERBOSE("downloading %s from UNR" % code )
        try:
            ts.append( Gts( ).get_unr( code ) )
        except:
            WARNING("Could not download time series for %s from UNR" % code )

    # return
    return ts
