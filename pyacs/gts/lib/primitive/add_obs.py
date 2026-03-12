###################################################################
def add_obs(self,date,NEUSNSESUCNECNUCEU,in_place=False,check=True,verbose=False):
###################################################################
    """
    Add observation(s) as DN, DE, DU to a time series.

    Parameters
    ----------
    date : float or list or ndarray
        Date(s) in decimal year.
    NEUSNSESUCNECNUCEU : list or ndarray
        Values to add: at least NEU (North, East, Up). Optional: SN, SE, SU, CNE, CNU, CEU
        (standard deviations and correlations). If not provided, SN=SE=SU=0.001, CNE=CNU=CEU=0.
    in_place : bool, optional
        If True, add to current Gts; if False, return a new Gts.
    check : bool, optional
        Check time order and duplicate dates.
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    Gts
        New Gts or the modified Gts if in_place.

    Notes
    -----
    If .data_xyz exists, it will be set to None for consistency.
    """
    # import 
    import numpy as np
    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG
    import pyacs.lib.astrotime as at

    # check argument

    np_date = np.array(date).reshape(-1)
    np_data = np.array(NEUSNSESUCNECNUCEU)


    if np_data.ndim==1:
        if np_data.shape[0] not in [3,6,9]:
            ERROR("second argument must be of length 3, 6 or 9. Its shape is %s" % str(NEUSNSESUCNECNUCEU.shape),exit=True)
        else:
            np_data = np_data.reshape(1,-1)
    else:
        if np_data.shape[1] not in [3, 6, 9]:
            ERROR("second argument must have 3, 6 or 9 columns. Its shape is %s" % str(NEUSNSESUCNECNUCEU.shape), exit=True)

    # fill data as 2d
    np_data_2d = np.zeros(( np_data.shape[0],10 ))
    if np_data.shape[1] ==3:
        np_data_2d[:,1:4] = np_data
        np_data_2d[:,4:7] = 1.E-3
    if np_data.shape[1] ==6:
        np_data_2d[:,1:7] = np_data
    if np_data.shape[1] ==9:
        np_data_2d[:,1:] = np_data

    if np_date.shape[0] != np_data_2d.shape[0]:
        ERROR("date and observation have incompatible shape. date: %s observation %s" % (str(np_date.shape),str(np_data.shape)),exit=True)

    np_data_2d[:,0] = np_date

    # update .data
    
    if in_place:
        new_gts = self
    else:
        new_gts = self.copy()
    
    # set .data_xyz to None
    
    new_gts.data_xyz = None

    VERBOSE("updating Gts for code %s with %d new observations" % (new_gts.code,np_data_2d.shape[0]))

    if new_gts.data is None:
        new_gts.data = np_data_2d
    else:
        if new_gts.data.shape[1] == 7:
            new_gts.data = np.vstack ( ( new_gts.data , np_data_2d ) )
        else:
            new_gts.data = np.vstack ( ( new_gts.data , np_data_2d ) )
        
    # check order and duplicate
    if check:
        # reorder
        new_gts.reorder(verbose=verbose)
        # check for duplicates
        new_gts.correct_duplicated_dates(action='correct',tol= .00000001, in_place=True,verbose=verbose)
    
    return(new_gts)
