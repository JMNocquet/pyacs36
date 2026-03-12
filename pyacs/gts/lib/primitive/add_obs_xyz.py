###################################################################
def add_obs_xyz(self,date,XYZSXSYSZCXYCXZCYZ,in_place=False,check=True, neu=True , verbose=False):
###################################################################
    """
    Add observation(s) as XYZ to a time series.

    Parameters
    ----------
    date : float or list or ndarray
        Date(s) in decimal year.
    XYZSXSYSZCXYCXZCYZ : list or ndarray
        Values to add: at least X, Y, Z. Optional: SX, SY, SZ, CXY, CXZ, CYZ (standard
        deviations and correlations). If not provided, SX=SY=SZ=0.001, CXY=CXZ=CYZ=0.
    in_place : bool, optional
        If True, add to current Gts; if False, return a new Gts.
    check : bool, optional
        Check time order, duplicate dates, and regenerate NEU (.data).
    neu : bool, optional
        Regenerate .data from updated .data_xyz.
    verbose : bool, optional
        Verbose mode.

    Returns
    -------
    Gts
        New Gts or the modified Gts if in_place.

    Notes
    -----
    By default .data is updated from .data_xyz and X0, Y0, Z0 are updated.
    """
    
    # import 
    import numpy as np

    # check argument
    
    if isinstance(date,list) :
        if  (len(list)>1) :
            date=np.array(date)
        else:
            date=date[0]
            
        
    if isinstance(date,np.ndarray):
        
        if not ( isinstance(XYZSXSYSZCXYCXZCYZ,np.ndarray) and (np.ndim(XYZSXSYSZCXYCXZCYZ) != 2 ) ):
            raise TypeError('!!! ERROR: various dates require the second argument to be a 2D numpy array')
            return(self)
    
    if isinstance(XYZSXSYSZCXYCXZCYZ,list):
        XYZSXSYSZCXYCXZCYZ=np.array(XYZSXSYSZCXYCXZCYZ)
    
    if np.ndim(XYZSXSYSZCXYCXZCYZ) == 1:
        if XYZSXSYSZCXYCXZCYZ.shape[0] not in [3,6,9]:
                raise TypeError('!!! ERROR: second argument must be of length 3, 6 or 9')
    
    if np.ndim(XYZSXSYSZCXYCXZCYZ) == 2:
        if XYZSXSYSZCXYCXZCYZ.shape[1] not in [3,6,9]:
                raise TypeError('!!! ERROR: second argument must be a 2D array with 3, 6 or 9 columns')
    
    # creates array to be stacked to current .data array
    
    if isinstance(date,float):
        # single obs provided
        
        new_data=np.zeros((1,10))
        new_data[0,4] = .001
        new_data[0,5] = .001
        new_data[0,6] = .001
        
        new_data[0,0] = date


        new_data[0,1] = XYZSXSYSZCXYCXZCYZ[0] 
        new_data[0,2] = XYZSXSYSZCXYCXZCYZ[1] 
        new_data[0,3] = XYZSXSYSZCXYCXZCYZ[2]
        
        if XYZSXSYSZCXYCXZCYZ.shape[0] > 3:

            new_data[0,4] = XYZSXSYSZCXYCXZCYZ[3] 
            new_data[0,5] = XYZSXSYSZCXYCXZCYZ[4] 
            new_data[0,6] = XYZSXSYSZCXYCXZCYZ[5]
             
        if XYZSXSYSZCXYCXZCYZ.shape[0] > 6:

            new_data[0,7] = XYZSXSYSZCXYCXZCYZ[6] 
            new_data[0,8] = XYZSXSYSZCXYCXZCYZ[7] 
            new_data[0,9] = XYZSXSYSZCXYCXZCYZ[8]
         
    if isinstance(date,np.ndarray):
        # various obs provided
        
        new_data=np.zeros((date.shape[0],10))
        new_data[:,4] = new_data[:,4] + .001
        new_data[:,5] = new_data[:,5] + .001
        new_data[:,6] = new_data[:,6] + .001

        if XYZSXSYSZCXYCXZCYZ.shape[1] > 3:

            new_data[:,4] = XYZSXSYSZCXYCXZCYZ[:,3] 
            new_data[:,5] = XYZSXSYSZCXYCXZCYZ[:,4] 
            new_data[:,6] = XYZSXSYSZCXYCXZCYZ[:,5]
             
        if XYZSXSYSZCXYCXZCYZ.shape[1] > 6:

            new_data[:,7] = XYZSXSYSZCXYCXZCYZ[:,6] 
            new_data[:,8] = XYZSXSYSZCXYCXZCYZ[:,7] 
            new_data[:,9] = XYZSXSYSZCXYCXZCYZ[:,8]

    # update .data_xyz
    
    if in_place:
        new_gts = self
    else:
        new_gts = self.copy()
    
    if verbose:
        print('-- updating Gts for code ',new_gts.code,' with ',new_data.shape[0], ' new entries.')

    if new_gts.data_xyz is None:
        new_gts.data_xyz = new_data
    else:
        if new_gts.data_xyz.shape[1] == 7:
            new_gts.data_xyz = np.vstack ( ( new_gts.data_xyz , new_data[:,:7] ) )
        else:
            new_gts.data_xyz = np.vstack ( ( new_gts.data_xyz , new_data ) )

    # re-generate .data
    if neu:
        new_gts.xyz2neu(corr=True)

    if check:
        # reorder
        new_gts.reorder(verbose=False)
        # check for duplicates
        new_gts.correct_duplicated_dates(action='correct',tol= .000001, in_place=True,verbose=verbose)
    
    return  new_gts 
