def apply_coseismic(self , date, coseismic_en_gmt_file=None, coseismic_up_gmt_file=None, verbose=False, add=True ):
    """Apply coseismic offsets from GMT psvelo files.

    Parameters
    ----------
    date : float
        Coseismic date in decimal year.
    coseismic_en_gmt_file : str, optional
        GMT psvelo file with East-North displacement (mm). Default is None.
    coseismic_up_gmt_file : str, optional
        GMT psvelo file with Up displacement (mm). Default is None.
    verbose : bool, optional
        Verbose mode. Default is False.
    add : bool, optional
        True = add offset, False = subtract. Default is True.

    Returns
    -------
    Sgts
        New Sgts with coseismic offsets applied.
    """

    # import
    import numpy as np
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import inspect
    VERBOSE("Running Sgts.%s" % inspect.currentframe().f_code.co_name)


    # reads en coseismic file
    DEN = np.array([])
    DEN_NAME = np.array([])

    if coseismic_en_gmt_file is not None:
        DEN = np.atleast_2d(np.genfromtxt( coseismic_en_gmt_file, usecols=(2,3) ) * 1E-3)
        DEN_NAME = np.atleast_1d(np.genfromtxt( coseismic_en_gmt_file, usecols=7, dtype=str ))
    # reads up coseismic file
    DUP = np.array([])
    DUP_NAME = np.array([])
    if coseismic_up_gmt_file is not None:
        DUP = np.atleast_2d(np.genfromtxt( coseismic_up_gmt_file, usecols=(3) ) * 1E-3)
        DUP_NAME = np.atleast_1d(np.genfromtxt( coseismic_up_gmt_file, usecols=7, dtype=str ))

    # apply

    new_Sgts = self.copy()

    for i in np.arange( DEN.shape[0] ):
        code = DEN_NAME[i]
        if new_Sgts.has_ts( code ):
            np_offset = np.array([[date,DEN[i,1],DEN[i,0],0.0]])
            new_Sgts.__dict__[ code ] = new_Sgts.__dict__[ code ].apply_offsets( np_offset , opposite=add)
            VERBOSE("Updating %s with offset at %.4lf E %.1lf N %.1lf " % (code,date,DEN[i,1]*1E3,DEN[i,0]*1E3))

    for i in np.arange( DUP.shape[0] ):
        code = DUP_NAME[i]
        if new_Sgts.has_ts( code ):
            np_offset = np.array([[date,0.0,0.0,DUP[i]]])
            new_Sgts.__dict__[ code ] = new_Sgts.__dict__[ code ].apply_offsets( np_offset , opposite=add)
            VERBOSE("Updating %s with offset at %.4lf U %.1lf " % (code,date,DUP[i]*1E3))


    # return

    return new_Sgts
