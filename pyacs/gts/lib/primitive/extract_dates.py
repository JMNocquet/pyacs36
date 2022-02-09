###################################################################
def extract_dates(self,dates, tol= 0.05 , in_place=False, verbose=True):
###################################################################
    """
    Returns a time series extracted for a given list of dates
    
    :param dates: dates either as a list or  1D numpy array of decimal dates
    :param tol: date tolerance in days to assert that two dates are equal (default 0.05 day)
    :param in_place: if True, will make change in place, if False, return s a new time series
    :param verbose: boolean, verbose mode

    """

    # import 
    import inspect
    import numpy as np
    import pyacs.gts

    # check data is not None
    from pyacs.gts.lib.errors import GtsInputDataNone
    
    try:
        if self.data is None:
            # raise exception
            raise GtsInputDataNone(inspect.stack()[0][3],__name__,self)
    except GtsInputDataNone as error:
        # print PYACS WARNING
        print( error )
        return( self )

    # working gts
    new_gts = self.copy()
  
    # case .data_xyz is None
    
    if new_gts.data_xyz is None:
        new_gts.neu2xyz(corr=True)

    else:
        # check data/data_xyz consistency
        try:
            if not new_gts.cdata(data=True):
                # raise exception
                from pyacs.gts.lib.errors import GtsCDataError
                raise GtsCDataError( inspect.stack()[0][3],__name__,self )
        except GtsCDataError as error:
            print( error )
            return( self )
    
  
    new_data=None
    
    # extract dates
        
    index = np.array( pyacs.gts.Gts.get_index_from_dates(dates, self.data, tol=tol) )
    
    if verbose:
        print('-- Extracting ',index.shape[0],' entries from Gts or code: ',self.code)
    
    if index.shape[0] > 0:
            new_data_xyz= self.data_xyz[index,:]
            new_sigma   = self.data[index,4:]
    else:
        new_data=None
        if verbose:
            print("-- time series ",self.code," does not have dates at the requested dates ")

  
    # handles outliers

    if new_data is not None:    
        ldate_outliers=self.data[:,0][self.outliers]
        lupdated_outliers=pyacs.gts.Gts.get_index_from_dates(ldate_outliers, new_data, tol=tol)
    else:
        lupdated_outliers = []

    if verbose:
        print('-- Transmitting ',len(lupdated_outliers),' outliers to the extracted Gts ')
    
    # case observations
    
    new_gts.data_xyz = new_data_xyz
      
    # handle outliers

    ldate_outliers=self.data[:,0][self.outliers]
    lupdated_outliers=pyacs.gts.Gts.get_index_from_dates(ldate_outliers, new_data_xyz, tol=0.05)
    
    # handles offsets_date
    
    upd_offsets=[]
    for offset_date in self.offsets_dates:
        if offset_date>=new_data_xyz[0,0] and offset_date<=new_data_xyz[-1,0]:
            upd_offsets.append(offset_date)
    
    # handles X0,Y0,Z0
    
    new_gts.X0 = new_data_xyz[0,1]
    new_gts.Y0 = new_data_xyz[0,2]
    new_gts.Z0 = new_data_xyz[0,3]
    
    # re-generate NEU time series
    new_gts.xyz2neu(corr=False)

    # re-populate the uncertainties columns
    new_gts.data[:,4:] = new_sigma
    
    # offsets & outliers
        
    new_gts.offsets_dates=upd_offsets
    new_gts.outliers=lupdated_outliers
    
    if in_place:
        self = new_gts
        return(self)
    else:
        return(new_gts)
