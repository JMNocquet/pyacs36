###############################################################################
def to_tspck(self, rounding='day', save=None, verbose=False):
###############################################################################
    """
    returns a pickle with the following numpy array

    T_OBS: observation array of dim 3. T_OBS[k,i,0] returns the East value at the k_th date for site i in mm
    np_names_t_obs: np_names_t_obs[i] is the code of site i
    np_obs_date_s: time series dates in pyacs seconds
    np_coo: for site i, np_coo[i,0] is longitude, np_coo[i,1] latitude and np_coo[i,1] ellipsoidal height

    :param sgts: Sgts instance
    :param rounding: 'day','hour','second'. all dates will be rounded to the nearest chosen day, hour or second. default is 'day'
    :param save: output file name for writing. Extension is .tspck
    :param verbose: boolean for verbose mode

    :return: object

    :note: tspck file to be read with read_tspck
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime as at

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.debug_message as DEBUG

    import inspect

    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)


    import pyeq.obs_tensor.sgts2obs_tensor

    T_OBS, np_names, np_obs_date_s = pyeq.obs_tensor.sgts2obs_tensor.sgts2tensor(self, rounding=rounding)

    # fills coordinates

    np_coo = np.zeros(( np_names.shape[0] , 3 ))
    for i in np.arange( np_names.shape[0] ):
        code = np_names[i]
        np_coo[i,0] = self.__dict__[code].lon
        np_coo[i,1] = self.__dict__[code].lat
        np_coo[i,2] = self.__dict__[code].h

    # create object

    class Tspck:
        pass

    wtspck = Tspck()

    wtspck.T_OBS = T_OBS
    wtspck.np_names = np_names
    wtspck.np_obs_date_s = np_obs_date_s
    wtspck.np_coo = np_coo

    # write file

    if save is not None:
        import pathlib
        if pathlib.Path(save).suffix != 'tspck':
            save = save + '.tspck'
        try:
            import pickle
            ofile = open( save, 'wb')
            pickle.dump(wtspck , ofile , pickle.DEFAULT_PROTOCOL)
            ofile.close()
        except:
            ERROR("Could not write %s" % save , exit=True)


    return wtspck


