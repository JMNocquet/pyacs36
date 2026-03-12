###############################################################################
def detrend_median(self, delta_day=None, periods=[], exclude_periods=[], verbose=False, auto=False, log_dir=None):
###############################################################################
    """Compute velocity from median of 1-year pair displacements (MIDAS-style).

    If the time series is shorter than one year, it is left unchanged.

    Parameters
    ----------
    delta_day : int or None, optional
        Time delta in days for pairs. None = one year, 0 = relax mode for campaign data. Default is None.
    periods : list, optional
        Periods to include for trend. Default is [].
    exclude_periods : list, optional
        Periods to exclude. Default is [].
    verbose : bool, optional
        If True, print progress. Default is False.
    auto : bool, optional
        If True, try delta_day=None then 0, else fallback to regular detrend. Default is False.
    log_dir : str, optional
        Directory for log file {self.code}_detrend_median.log. Default is None.

    Returns
    -------
    Gts or None
        Detrended Gts, or None if series shorter than 1 year.
    """

    # import
    import numpy as np
    from pyacs.gts.Gts import Gts
    import inspect
    import logging
    import os
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    # Setup Gts-specific logging if log_dir is provided
    gts_logger = None
    if log_dir is not None:
        # Create log directory if it doesn't exist
        os.makedirs(log_dir, exist_ok=True)
        
        # Create Gts-specific log file
        log_file = os.path.join(log_dir, f"{self.code}.log")
        
        # Configure Gts-specific logger
        gts_logger = logging.getLogger(f'pyacs.gts.{self.code}')
        if not gts_logger.handlers:
            gts_logger.setLevel(logging.INFO)
            # Prevent propagation to parent loggers to avoid duplicate logging
            gts_logger.propagate = False
            file_handler = logging.FileHandler(log_file, mode='a')  # append mode
            file_handler.setLevel(logging.INFO)
            # Custom formatter with [Gts.code] prefix
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - [%(name)s] %(message)s')
            file_handler.setFormatter(formatter)
            gts_logger.addHandler(file_handler)
        
        # Log method start with [Gts.code] prefix
        gts_logger.info(f"[{self.code}] [detrend_median] Starting method")
        gts_logger.info(f"[{self.code}] [detrend_median] Parameters: delta_day={delta_day}, auto={auto}")
        gts_logger.info(f"[{self.code}] [detrend_median] Time series: {len(self.data)} points from {self.data[0, 0]:.3f} to {self.data[-1, 0]:.3f}")

    # after this method .data  and .data_xyz are not consistent so .data_xyz is set to None
    #self.data_xyz = None

    ###########################################################################
    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone

    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3], __name__, self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        if gts_logger:
            gts_logger.error(f"[{self.code}] [detrend_median] Input data is None")
        print(error)
        return (self)
    ###########################################################################

    ###########################################################################
    # 1-yr threshold
    ###########################################################################

    if (self.data[-1, 0] - self.data[0, 0]) < 1.:
        warning_msg = f"time series shorter than 1 year for side: {self.code}"
        if gts_logger:
            gts_logger.warning(f"[{self.code}] [detrend_median] Time series shorter than 1 year")
        WARNING(warning_msg)
        return (None)

    ###########################################################################
    # run the appropriate estimator using autoreference if auto is True
    ###########################################################################

    if auto:
        if gts_logger:
            gts_logger.info(f"[{self.code}] [detrend_median] Auto mode enabled, trying different delta_day values")

        tts = self.detrend_median(delta_day=None, \
                                  periods=periods, exclude_periods=exclude_periods, \
                                  verbose=verbose, auto=False, log_dir=log_dir)

        if tts is None:
            if gts_logger:
                gts_logger.info(f"[{self.code}] [detrend_median] delta_day=None failed, trying delta_day=0")
            tts = self.detrend_median(delta_day=0, \
                                      periods=periods, exclude_periods=exclude_periods, \
                                      verbose=verbose, auto=False, log_dir=log_dir)

        if tts is None:
            if gts_logger:
                gts_logger.info(f"[{self.code}] [detrend_median] delta_day=0 failed, trying regular detrend")
            tts = self.detrend(periods=periods, exclude_periods=exclude_periods)

        if gts_logger:
            if tts is not None:
                gts_logger.info(f"[{self.code}] [detrend_median] Auto mode completed successfully")
            else:
                gts_logger.warning(f"[{self.code}] [detrend_median] Auto mode failed")

        return (tts)

    ###########################################################################
    # start main detrend_median
    ###########################################################################

    import pyacs.lib.astrotime

    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend_median] Starting main processing")

    tmp_ts = self.remove_outliers()

    if periods != []:
        tmp_ts = tmp_ts.extract_periods(periods)
    if exclude_periods != []:
        tmp_ts = tmp_ts.exclude_periods(periods)

    # minimum number of pair velocity estimates to consider the result

    min_vel_estimate = 100

    ###########################################################################
    # MIDAS-like method, using 1-year pair data for velocity estimate
    ###########################################################################

    if delta_day is None:
        if gts_logger:
            gts_logger.info(f"[{self.code}] [detrend_median] Using 1-year pair method")

        duration_in_year = tmp_ts.data[-1, 0] - tmp_ts.data[0, 0]

        if duration_in_year >= 1.:

            H_year = {}

            # creates H_year arrays

            for year in np.arange(int(tmp_ts.data[0, 0]), int(tmp_ts.data[-1, 0]) + 1):
                H_year[year] = np.zeros((365, 3))
                H_year[year][:, :] = np.nan

            # deal with dates

            np_mjd = pyacs.lib.astrotime.decyear2mjd(tmp_ts.data[:, 0])
            (np_doy, np_year, _np_ut) = pyacs.lib.astrotime.mjd2dayno(np_mjd)

            np_doy = np_doy.astype(int)
            np_year = np_year.astype(int)

            # fill H_year arrays

            for i in np.arange(np_doy.shape[0]):
                if np_doy[i] != 366:
                    H_year[np_year[i]][np_doy[i] - 1, :] = tmp_ts.data[i, 1:4]

            # stack the velocities

            np_vel = np.zeros((365 * (len(H_year) - 1), 3))
            np_vel[:, :] = np.nan

            i = 0
            for year in sorted(H_year.keys())[0:-1]:
                np_vel[i * 365:(i + 1) * 365] = H_year[year + 1] - H_year[year]
                i = i + 1

            # test whether there are at least 3 non-Nan numbers

            if np.count_nonzero(~np.isnan(np_vel)) < min_vel_estimate:
                if gts_logger:
                    gts_logger.warning(f"[{self.code}] [detrend_median] Insufficient velocity estimates: {np.count_nonzero(~np.isnan(np_vel))} < {min_vel_estimate}")
                return (None)

            # calculates the median velocity

            med_vel = np.nanmedian(np_vel, axis=0)

            [med_vel_north, med_vel_east, med_vel_up] = med_vel

            # calculate uncertainty and refined value of vel using step 2 from Blewitt et al. (2016) p. 2057
            # remove Nan

            np_vel_north = np_vel[:, 0][np.logical_not(np.isnan(np_vel[:, 0]))]
            np_vel_east = np_vel[:, 1][np.logical_not(np.isnan(np_vel[:, 1]))]
            np_vel_up = np_vel[:, 2][np.logical_not(np.isnan(np_vel[:, 2]))]

            # calculates sigma

            fabs_vn = np.fabs(np_vel_north - med_vel_north)
            fabs_ve = np.fabs(np_vel_east - med_vel_east)
            fabs_vu = np.fabs(np_vel_up - med_vel_up)

            sigma_vn = 1.4826 * np.median(fabs_vn)
            sigma_ve = 1.4826 * np.median(fabs_ve)
            sigma_vu = 1.4826 * np.median(fabs_vu)

            # removes all vel > 1.4826 * med_vel

            np_vel_north_cln = np_vel_north[np.where(fabs_vn <= 2 * sigma_vn)]
            np_vel_east_cln = np_vel_east[np.where(fabs_ve <= 2 * sigma_ve)]
            np_vel_up_cln = np_vel_up[np.where(fabs_vu <= 2 * sigma_vu)]

            # get the new vel estimates

            med_vn = np.median(np_vel_north_cln)
            med_ve = np.median(np_vel_east_cln)
            med_vu = np.median(np_vel_up_cln)

            # get the new sigma

            med_sigma_vn = 1.4826 * np.median(np.fabs(np_vel_north_cln - med_vn))
            med_sigma_ve = 1.4826 * np.median(np.fabs(np_vel_east_cln - med_ve))
            med_sigma_vu = 1.4826 * np.median(np.fabs(np_vel_up_cln - med_vu))

            # empirical factor for correlated noise
            ef = 3.

            sigma_vn = 1.2533 * med_sigma_vn / np.sqrt(fabs_vn.shape[0]) * ef
            sigma_ve = 1.2533 * med_sigma_ve / np.sqrt(fabs_ve.shape[0]) * ef
            sigma_vu = 1.2533 * med_sigma_vu / np.sqrt(fabs_vu.shape[0]) * ef

            # final med_vel

            med_vel = np.array([med_vn, med_ve, med_vu, sigma_vn, sigma_ve, sigma_vu])

        else:
            # time series has less than a year of data
            # velocity is set to 0.
            if gts_logger:
                gts_logger.warning(f"[{self.code}] [detrend_median] Time series has less than 1 year of data, setting velocity to 0")
            med_vel = np.array([0., 0., 0.])

    # end case 1-yr

    elif (isinstance(delta_day, int) and delta_day == 0):

        ###########################################################################
        # case campaign, close to relaxed method in MIDAS
        ###########################################################################
        if gts_logger:
            gts_logger.info(f"[{self.code}] [detrend_median] Using campaign mode (delta_day=0)")

        H_year = {}
        lvel = []

        # creates H_year arrays

        for year in np.arange(int(tmp_ts.data[0, 0]), int(tmp_ts.data[-1, 0]) + 1):
            lindex = np.where(tmp_ts.data.astype(int)[:, 0] == year)

            #            tt = tmp_ts.extract_periods([float(year),float(year)+0.999999])
            #            if tt.data is not None:
            #            H_year[year] = tt.data[:,:4]
            H_year[year] = tmp_ts.data[lindex[0], :4]

        # removes years with no available data
        # H_year = {key:val for key, val in H_year.items() if val is not None}
        # calculates velocity values for all data pairs separated by more than one year
        #        while H_year != {}:

        for year_start in sorted(list(H_year.keys())):
            #            np_start = H_year[year_start]

            for year_end in sorted(list(H_year.keys())):

                ok = False

                for i in np.arange(H_year[year_end].shape[0]):
                    for j in np.arange(H_year[year_start].shape[0]):

                        # check wether the pair date includes a given offset date

                        include_offset_date = False

                        for odate in self.offsets_dates:

                            if (odate >= H_year[year_start][j, 0]) and (odate <= H_year[year_end][i, 0]):
                                include_offset_date = True

                        if not include_offset_date:

                            # displacement
                            v = (H_year[year_end][i] - H_year[year_start][j])[1:]
                            # delta_year
                            delta_year = H_year[year_end][i, 0] - H_year[year_start][j, 0]
                            # append to lvel
                            if delta_year > 0.90:
                                lvel.append(v / delta_year)
                                ok = True
                                if verbose:
                                    print(
                                        "-- adding %.5lf - %.5lf " % (H_year[year_start][j, 0], H_year[year_end][i, 0]))
                                    print(i, j, v / delta_year)
                                    print(H_year[year_end][i], H_year[year_start][j])
                # if some pairs have already been found, stop here to mitigate the effect of offsets
                # if ok:
                #    break

            del H_year[year_start]

        # calculates the median velocity

        np_vel = np.array(lvel)

        np_vel_north = np_vel[:, 0]
        np_vel_east = np_vel[:, 1]
        np_vel_up = np_vel[:, 2]

        med_vel = np.median(np_vel, axis=0)

        [med_vel_north, med_vel_east, med_vel_up] = med_vel

        # calculates sigma

        fabs_vn = np.fabs(np_vel_north - med_vel_north)
        fabs_ve = np.fabs(np_vel_east - med_vel_east)
        fabs_vu = np.fabs(np_vel_up - med_vel_up)

        sigma_vn = 1.4826 * np.median(fabs_vn)
        sigma_ve = 1.4826 * np.median(fabs_ve)
        sigma_vu = 1.4826 * np.median(fabs_vu)

        # removes all vel > 1.4826 * med_vel

        np_vel_north_cln = np_vel_north[np.where(fabs_vn <= 2 * sigma_vn)]
        np_vel_east_cln = np_vel_east[np.where(fabs_ve <= 2 * sigma_ve)]
        np_vel_up_cln = np_vel_up[np.where(fabs_vu <= 2 * sigma_vu)]

        # get the new vel estimates

        med_vn = np.median(np_vel_north_cln)
        med_ve = np.median(np_vel_east_cln)
        med_vu = np.median(np_vel_up_cln)

        # get the new sigma

        med_sigma_vn = 1.4826 * np.median(np.fabs(np_vel_north_cln - med_vn))
        med_sigma_ve = 1.4826 * np.median(np.fabs(np_vel_east_cln - med_ve))
        med_sigma_vu = 1.4826 * np.median(np.fabs(np_vel_up_cln - med_vu))

        # empirical factor for correlated noise
        ef = 3.

        sigma_vn = 1.2533 * med_sigma_vn / np.sqrt(fabs_vn.shape[0]) * ef
        sigma_ve = 1.2533 * med_sigma_ve / np.sqrt(fabs_ve.shape[0]) * ef
        sigma_vu = 1.2533 * med_sigma_vu / np.sqrt(fabs_vu.shape[0]) * ef

        # final med_vel

        med_vel = np.array([med_vn, med_ve, med_vu, sigma_vn, sigma_ve, sigma_vu])

    else:

        ###########################################################################
        # case delta_day with integer value
        ###########################################################################
        if gts_logger:
            gts_logger.info(f"[{self.code}] [detrend_median] Using custom delta_day={delta_day}")

        def delta(n):
            """
            Simple delta function for integer
            """
            if not n == 0:
                return (1)
            else:
                return (0)

        # first form the day vector

        np_mjd_data = list(map(int, pyacs.lib.astrotime.decyear2mjd(tmp_ts.data[:, 0])))
        # so the index are
        i_np_mjd_data = np_mjd_data - np.min(np_mjd_data)

        # the void array filled with np.nan

        void = np.zeros((np.max(np_mjd_data) - np.min(np_mjd_data) + 1, 3)) * np.nan

        # I fill void with the time series at existing dates

        void[i_np_mjd_data] = tmp_ts.data[:, 1:4]

        # now reshape the vector for easy pair differentiation

        # if TS = void[ : , (0,1,2) ] easy diff would be achieved by working on a array W = TS.reshape(-1,delta_day)
        # however, this requires to complete the number of lines of TS with np.nan to be able to reshape it

        CTS = np.zeros(
            (delta_day - np.mod(void.shape[0], delta_day)) * delta(np.mod(void.shape[0], delta_day))) * np.nan

        # init med_vel

        med_vel = np.array([0., 0., 0.])

        to_year = 365.25 / float(delta_day)

        for i in [0, 1, 2]:
            med_vel[i] = np.nanmedian(np.diff(np.append(void[:, i], CTS).reshape(-1, delta_day).T, axis=1)) * to_year

    ###########################################################################
    # return detrended time series
    ###########################################################################


    new_gts = self.remove_velocity(med_vel)
    new_gts.outliers = self.outliers
    new_gts.offsets_values = self.offsets_values
    new_gts.velocity = med_vel

    VERBOSE("Velocity estimates (N,E,U, mm/yr) from detrend_median for %s: [%.3f, %.3f, %.3f]" % (self.code, med_vel[0]*1.E3, med_vel[1]*1.E3, med_vel[2]*1.E3))

    # Log final results
    if gts_logger:
        gts_logger.info(f"[{self.code}] [detrend_median] Velocity estimates (N,E,U, mm/yr): [{med_vel[0]*1.E3:.3f}, {med_vel[1]*1.E3:.3f}, {med_vel[2]*1.E3:.3f}]")
        if len(med_vel) > 3:
            gts_logger.info(f"[{self.code}] [detrend_median] Velocity uncertainties (N,E,U, mm/yr): [{med_vel[3]*1.E3:.3f}, {med_vel[4]*1.E3:.3f}, {med_vel[5]*1.E3:.3f}]")
        gts_logger.info(f"[{self.code}] [detrend_median] Method completed successfully")

    return (self.remove_velocity(med_vel))
