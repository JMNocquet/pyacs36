###################################################################
def add_obs(self,date,NEUSNSESUCNECNUCEU,in_place=False,check=True,verbose=False):
###################################################################
    """
    Adds observation(s) as DN,DE,DU to a time series

    :param date: date in decimal year. float, a list or 1D numpy array
    :param NEUSNSESUCNECNUCEU: value to be added in the Gts, provided as a list, a 1D numpy array or a 2D numpy array.
                               requires at least NEU: North, East, UP values
                               optional: SN, SE, SU, CNE, CNU, CEU: standard deviations and correlation coefficient between North, East and Up components.
                               If not provided, SN=SE=SU=0.001 (1 mm) and CNE=CNU=CEU=0

    :param in_place: boolean, if True add_obs to the current Gts, if False, returns a new Gts
    :param check: check time order and duplicate dates
    :param verbose: verbose mode
    :return: new Gts or the modified Gts if in_place
    :note: if it exists, .data_xyz will be set to None for consistency.
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
        
        if not ( isinstance(NEUSNSESUCNECNUCEU,np.ndarray) and (np.ndim(NEUSNSESUCNECNUCEU) != 2 ) ):
            raise TypeError('!!! ERROR: various dates require the second argument to be a 2D numpy array')
            return(self)
    
    if isinstance(NEUSNSESUCNECNUCEU,list):
        NEUSNSESUCNECNUCEU=np.array(NEUSNSESUCNECNUCEU)
    
    if np.ndim(NEUSNSESUCNECNUCEU) == 1:
        if NEUSNSESUCNECNUCEU.shape[0] not in [3,6,9]:
                raise TypeError('!!! ERROR: second argument must be of length 3, 6 or 9')
    
    if np.ndim(NEUSNSESUCNECNUCEU) == 2:
        if NEUSNSESUCNECNUCEU.shape[1] not in [3,6,9]:
                raise TypeError('!!! ERROR: second argument must be a 2D array with 3, 6 or 9 columns')
    
    # creates array to be stacked to current .data array
    
    if isinstance(date,float):
        # single obs provided
        
        new_data=np.zeros((1,10))
        new_data[0,4] = .001
        new_data[0,5] = .001
        new_data[0,6] = .001
        
        new_data[0,0] = date


        new_data[0,1] = NEUSNSESUCNECNUCEU[0] 
        new_data[0,2] = NEUSNSESUCNECNUCEU[1] 
        new_data[0,3] = NEUSNSESUCNECNUCEU[2]
        
        if NEUSNSESUCNECNUCEU.shape[0] > 3:

            new_data[0,4] = NEUSNSESUCNECNUCEU[3] 
            new_data[0,5] = NEUSNSESUCNECNUCEU[4] 
            new_data[0,6] = NEUSNSESUCNECNUCEU[5]
             
        if NEUSNSESUCNECNUCEU.shape[0] > 6:

            new_data[0,7] = NEUSNSESUCNECNUCEU[6] 
            new_data[0,8] = NEUSNSESUCNECNUCEU[7] 
            new_data[0,9] = NEUSNSESUCNECNUCEU[8]
         
    if isinstance(date,np.ndarray):
        # various obs provided
        
        new_data=np.zeros((date.shape[0],10))
        new_data[:,4] = new_data[:,4] + .001
        new_data[:,5] = new_data[:,5] + .001
        new_data[:,6] = new_data[:,6] + .001

        if NEUSNSESUCNECNUCEU.shape[1] > 3:

            new_data[:,4] = NEUSNSESUCNECNUCEU[:,3] 
            new_data[:,5] = NEUSNSESUCNECNUCEU[:,4] 
            new_data[:,6] = NEUSNSESUCNECNUCEU[:,5]
             
        if NEUSNSESUCNECNUCEU.shape[1] > 6:

            new_data[:,7] = NEUSNSESUCNECNUCEU[:,6] 
            new_data[:,8] = NEUSNSESUCNECNUCEU[:,7] 
            new_data[:,9] = NEUSNSESUCNECNUCEU[:,8]

    # update .data
    
    if in_place:
        new_gts = self
    else:
        new_gts = self.copy()
    
    ### set .data_xyz to None
    
    new_gts.data_xyz = None

    if verbose:
        print('-- updating Gts for code ',new_gts.code,' with ',new_data.shape[0], ' new entries.')

    if new_gts.data is None:
        new_gts.data = new_data
    else:
        if new_gts.data.shape[1] == 7:
            new_gts.data = np.vstack ( ( new_gts.data , new_data[:,:7] ) )
        else:
            new_gts.data = np.vstack ( ( new_gts.data , new_data ) )
        
    # check order and duplicate
    if check:
        # reorder
        new_gts.reorder(verbose=verbose)
        # check for duplicates
        new_gts.correct_duplicated_dates(action='correct',tol= .00000001, in_place=True,verbose=verbose)
    
    return(new_gts)
