###############################################################################
def to_obs_tensor(self, rounding='day', verbose=False):
###############################################################################
    """Return ENU data as a 3D tensor and site code array.

    Parameters
    ----------
    rounding : str, optional
        'day', 'hour', or 'second' for date rounding. Default is 'day'.
    verbose : bool, optional
        Verbose mode. Default is False.

    Returns
    -------
    tuple
        (T_OBS, np_names, np_date_s); T_OBS[k,i,0] = East at date k, site i (mm).

    Notes
    -----
    T_OBS is 3D (dates x sites x components), values in mm.
    """
# TODO : it looks like it is not always working. But from pyeq.obs_tensor.sgts2obs_tensor import sgts2tensor works
    # import
    import numpy as np
    import pyacs.lib.astrotime as at

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.debug_message as DEBUG


    from tqdm import tqdm


    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)

    # np print option for debug
    # np.set_printoptions(precision=2 , suppress=True)

    if verbose:
        # DEBUG MODE
        logging.getLogger("my_logger").setLevel(logging.DEBUG)

    # remove time series with less than 3 observation
    for code in self.lcode():
        if (self.__dict__[code].data) is None or (self.__dict__[code].data.shape[0]<3):
            MESSAGE("removing time series for %s because <3 epochs available" % code)
            self.delts( code )


    # np_names
    np_names = np.array(sorted(self.lcode()))

    # get all dates
    np_seconds_sol = np.array([], dtype=np.int64)
    for i in np.arange(np_names.shape[0]):
        # code
        code = np_names[i]
        # get the gts for current site
        wts = self.__dict__[code]
        # update np of all dates
        np_seconds_sol = np.unique(
            np.sort(np.append(np_seconds_sol, at.decyear2seconds(wts.data[:, 0], rounding=rounding))))

    # initialize T
    T = np.full((np_seconds_sol.shape[0], np_names.shape[0], 6), np.nan)

    DEBUG(("shape of T: %s " % (T.shape,)))

    # loop on gts in sgts


    #for i in np.arange(np_names.shape[0]):
    for site in tqdm(np.arange(np_names.shape[0]), desc='to_obs_tensor'):

        # code
        code = np_names[i]

        DEBUG("%04d / %s " % (i, code))

        # get the gts for current site
        wts = self.__dict__[code]

        # date of the current wts in seconds
        np_seconds_site = at.decyear2seconds(wts.data[:, 0], rounding=rounding)

        # check for non duplicated dates
        if np.min(np.diff(np_seconds_site)) <= 0:
            ERROR("there is a date problem in gts: %s" % code)
            ERROR("most probably, your round option (%s) leads to two successive equal dates" % (rounding))
            ERROR(("%d" % np.min(np.diff(np_seconds_site))))
            print(np.diff(np_seconds_site))
            print(np_seconds_site)
            ERROR("", exit=True)

        # current gts date index in T

        lindex = np.intersect1d(np_seconds_sol, np_seconds_site, assume_unique=True, return_indices=True)[1]

        DEBUG("number of observation in ts: %s:%d" % (code, wts.data.shape[0]))
        DEBUG("number of observation in T : %d" % (T.shape[0]))
        DEBUG("number of observation to be filled in T from lindex : %d" % (lindex.shape[0]))
        DEBUG("size of T slice to be filled: %s" % (T[lindex, i, :].shape,))

        # fill T - ENU - SE SN SU
        T[lindex, i, 0] = wts.data[:, 2] * 1000.
        T[lindex, i, 1] = wts.data[:, 1] * 1000.
        T[lindex, i, 2] = wts.data[:, 3] * 1000.

        T[lindex, i, 3] = wts.data[:, 5] * 1000.
        T[lindex, i, 4] = wts.data[:, 4] * 1000.
        T[lindex, i, 5] = wts.data[:, 6] * 1000.

    # return
    return T, np_names, np_seconds_sol

