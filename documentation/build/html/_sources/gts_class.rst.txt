Gts class
=============


``Gts`` class description
*****************************
.. autoclass:: pyacs.gts.Gts.Gts

    The individual Geodetic Time Series Class

    The Gts implemented in PYACS has the following attributes:

    Mandatory attributes:
        - data:        a 2D numpy array with 10 columns: dec.year, N, E, U, S_N, S_E, S_U, S_NE, S_NU, S_EU
        - code:            station 4-letters code
    Coordinates attributes
        - lon,lat,h:       approximate longitude, latitude (geodetic, deg.dec) and ellipsoidal height (m)
        - X0,Y0,Z0         XYZ reference position in the Geocentric Frame. N,E,U are considered with respect to X0,Y0,Z0
        - t0               reference date in decimal year for X0,Y0,Z0
    Not persisting attributes
        - data_xyz:        a 2D numpy array with 10 columns: dec.year, X, Y, Z, SX, SY, SZ, S_XY, S_XZ, S_YZ
            because many Gts methods are applied on NEU components, .data_xyz is often set to None.
            it can however be rebuilt using the neu2xyz method
    Attributes populated after some analysis
        - outliers:        list of index of outliers in a time series (all components)
        - offsets_values:  a 2D numpy array with 7 columns: dec.year N, E, U, S_N, S_E, S_U
        - offsets_dates:   a list of dates for offsets
        - velocity:        a 1D numpy array with 6 columns: vel_N, vel_E, vel_U, S_vel_N, S_vel_E, S_vel_U
        - annual:          a 1D numpy array with 6 columns: Amplitude_N, Phase_N, Amplitude_E, Phase_E, Amplitude_U, Phase_U
        - semi_annual:     a 1D numpy array with 6 columns: Amplitude_N, Phase_N, Amplitude_E, Phase_E, Amplitude_U, Phase_U
    Metadata attributes
        - ifile:          original input file of the time series
        - log:            log of operations
        - metadata:       any information the analyst would like to be recorded

    Units Conventions
        - dates are in decimal year
        - coordinates are in meters
        - phases are in radians
 