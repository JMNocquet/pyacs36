###############################################################################
def detrend_seasonal_median(self, wl=11, in_place=False, verbose=False):
###############################################################################
    """
    Calculates a velocity using the median of pair of displacements exactly separated by one year, inspired from MIDAS and then removes repeating yearly signal
    If the time series has less than three years of data, then the time series is kept untouched.

    """


    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect


    # after this method .data  and .data_xyz are not consistent so .data_xyz is set to None
    self.data_xyz = None

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print(error)
        return (self)
    ###########################################################################

    import pyacs.lib.astrotime

    duration_in_year = self.data[-1, 0] - self.data[0, 0]

    if duration_in_year >= 3.:

        H_year = {}

        # creates H_year arrays

        for year in np.arange(int(self.data[0, 0]), int(self.data[-1, 0]) + 1):
            H_year[year] = np.zeros((365, 3))
            H_year[year][:, :] = np.nan

        # deal with dates

        np_mjd = pyacs.lib.astrotime.decyear2mjd(self.data[:, 0])
        (np_doy, np_year, _np_ut) = pyacs.lib.astrotime.mjd2dayno(np_mjd)

        np_doy = np_doy.astype(int)
        np_year = np_year.astype(int)

        # fill H_year arrays

        for i in np.arange(np_doy.shape[0]):
            if np_doy[i] != 366:
                H_year[np_year[i]][np_doy[i] - 1, :] = self.data[i, 1:4]

        # stack the velocities

        np_vel = np.zeros((365 * (len(H_year) - 1), 3))
        np_vel[:, :] = np.nan

        i = 0
        for year in sorted(H_year.keys())[0:-1]:
            np_vel[i * 365:(i + 1) * 365] = H_year[year + 1] - H_year[year]
            i = i + 1
        # calculates the median velocity

        med_vel = np.nanmedian(np_vel, axis=0)

        # return detrended time series

        detrended = self.remove_velocity(med_vel)

        H_year = {}

        # creates H_year arrays

        for year in np.arange(int(detrended.data[0, 0]), int(detrended.data[-1, 0]) + 1):
            H_year[year] = np.zeros((365, 3))
            H_year[year][:, :] = np.nan

        # deal with dates

        np_mjd = pyacs.lib.astrotime.decyear2mjd(detrended.data[:, 0])
        (np_doy, np_year, _np_ut) = pyacs.lib.astrotime.mjd2dayno(np_mjd)

        np_doy = np_doy.astype(int)
        np_year = np_year.astype(int)

        # fill H_year arrays

        for i in np.arange(np_doy.shape[0]):
            if np_doy[i] != 366:
                H_year[np_year[i]][np_doy[i] - 1, :] = detrended.data[i, 1:4]

        # center all H_year arrays

        for year in sorted(H_year.keys()):
            H_year[year] = H_year[year] - np.nanmedian(H_year[year], axis=0)
        #            plt.plot(H_year[year][:,1])
        # create the median daily signal

        A = np.array(list(H_year.values()))

        np_doy_median_signal = np.nanmedian(A[:, :, :], axis=0)

        #        plt.plot(np_doy_median_signal[:,1],'ro')

        # run a median filter on it

        np_doy_median_signal_3 = np.vstack((np_doy_median_signal, np_doy_median_signal, np_doy_median_signal))

        import scipy.signal
        np_doy_median_signal[:, 0] = scipy.signal.medfilt(np_doy_median_signal_3[:, 0], wl)[365:2 * 365]
        np_doy_median_signal[:, 1] = scipy.signal.medfilt(np_doy_median_signal_3[:, 1], wl)[365:2 * 365]
        np_doy_median_signal[:, 2] = scipy.signal.medfilt(np_doy_median_signal_3[:, 2], wl)[365:2 * 365]

        #        plt.plot(np_doy_median_signal[:,1],'bo')

        # remove it from the detrended time series

        detrended_seasonal = detrended.copy()

        # loop on np_doy

        for i in np.arange(detrended_seasonal.data.shape[0]):
            if np_doy[i] != 366:
                detrended_seasonal.data[i, 1:4] = detrended_seasonal.data[i, 1:4] - np_doy_median_signal[np_doy[i] - 1,
                                                                                    :]



    else:
        # time series is shorter than minimum

        detrended_seasonal = self.copy()

    return (detrended_seasonal)
