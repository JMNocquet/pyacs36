###############################################################################
def to_tsnpz(self, ofile, rounding='day', verbose=False):
###############################################################################
    """Write Sgts to an npz file.

    Arrays saved: T_OBS (3D, mm), np_names_t_obs (site codes), np_obs_date_s (pyacs seconds),
    np_coo (lon, lat, height per site).

    Parameters
    ----------
    ofile : str
        Output file path (.npz).
    rounding : str, optional
        'day', 'hour', or 'second' for date rounding. Default is 'day'.
    verbose : bool, optional
        Verbose mode. Default is False.

    Returns
    -------
    None
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


    from pyacs.gts.lib.tensor_ts.sgts2obs_tensor import sgts2tensor

    T_OBS, np_names, np_obs_date_s = sgts2tensor(self, rounding=rounding)

    # fills coordinates

    np_coo = np.zeros(( np_names.shape[0] , 3 ))
    for i in np.arange( np_names.shape[0] ):
        code = np_names[i]
        np_coo[i,0] = self.__dict__[code].lon
        np_coo[i,1] = self.__dict__[code].lat
        np_coo[i,2] = self.__dict__[code].h


    # write file

    import pathlib
    if pathlib.Path( ofile ).suffix == 'npz':
        ofile = ofile.replace('.npz', '')
    MESSAGE("Writing % s " % ofile)
    try:
        np.savez( ofile, T_OBS=T_OBS, names=np_names, date_s=np_obs_date_s, coo=np_coo )
    except:
        ERROR("Could not write %s" % ofile , exit=True)

    try:
        npzfile = np.load( ofile )
    except:
        ERROR("Could not properly load %s" % ofile, exit=True)

    return npzfile


