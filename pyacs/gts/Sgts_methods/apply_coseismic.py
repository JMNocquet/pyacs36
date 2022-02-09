def apply_coseismic(self , date, coseismic_en_gmt_file=None, coseismic_up_gmt_file=None, verbose=False, add=True ):
    """

    :param self:
    :param date: date of the coseismic offset to be applied, in decimal year
    :param coseismic_en_gmt_file: coseismic eat-north displacement (in mm) gmt psvelo file
    :param coseismic_up_gmt_file: coseismic up displacement (in mm) gmt psvelo file
    :param verbose:
    :param add: boolean, add (True, default) or substract (False)
    :return: new Sgts with coseismic offsets applied
    """

    # import
    import numpy as np

    # reads en coseismic file
    DEN = np.array([])
    DEN_NAME = np.array([])

    if coseismic_en_gmt_file is not None:
        DEN = np.genfromtxt( coseismic_en_gmt_file, usecols=(2,3) ) * 1E-3
        DEN_NAME = np.genfromtxt( coseismic_en_gmt_file, usecols=7, dtype=str )
    # reads up coseismic file
    DUP = np.array([])
    DUP_NAME = np.array([])
    if coseismic_up_gmt_file is not None:
        DUP = np.genfromtxt( coseismic_up_gmt_file, usecols=(3) ) * 1E-3
        DUP_NAME = np.genfromtxt( coseismic_up_gmt_file, usecols=7, dtype=str )

    # apply

    new_Sgts = self.copy()

    for i in np.arange( DEN.shape[0] ):
        code = DEN_NAME[i]
        if new_Sgts.has_ts( code ):
            np_offset = np.array([[date,DEN[i,1],DEN[i,0],0.0]])
            new_Sgts.__dict__[ code ] = new_Sgts.__dict__[ code ].apply_offsets( np_offset , opposite=add)
            if verbose:
                print("Updating %s with offset at %.4lf E %.1lf N %.1lf " % (code,date,DEN[i,1]*1E3,DEN[i,0]*1E3))

    for i in np.arange( DUP.shape[0] ):
        code = DUP_NAME[i]
        if new_Sgts.has_ts( code ):
            np_offset = np.array([[date,0.0,0.0,DUP[i]]])
            new_Sgts.__dict__[ code ] = new_Sgts.__dict__[ code ].apply_offsets( np_offset , opposite=add)
            if verbose:
                print("Updating %s with offset at %.4lf U %.1lf " % (code,date,DUP[i]*1E3))


    # return

    return new_Sgts
