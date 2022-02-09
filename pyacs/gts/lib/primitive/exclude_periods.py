
###################################################################
def exclude_periods(self,lperiod,in_place=False,verbose=False):
###################################################################
    """
    exclude periods of a Gts
    
    :param lperiod: a list [start_date,end_date] or a list of periods e.g. periods=[[2000.1,2003.5],[2009.3,2010.8]]
    :param in_place: if True, will make change in place, if False, return s a new time series
    
    :note 1: X0,Y0,Z0 attributes will be changed if necessary
    :note 2: handles both .data and .data_xyz
    """

    ### import 
    import inspect
    import numpy as np
    import pyacs.gts.Gts
    import pyacs.lib.utils


    ### working gts
    new_gts = self.copy()

    ### case .data_xyz is None
    
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

    # ensure lperiod is a list of lists
    lperiod = pyacs.lib.utils.__ensure_list_of_list(lperiod)

    
    # check lperiod is not [[]]
    if lperiod==[[]]:
        if in_place:
            return(self)
        else:
            return( self.copy() )
    
    # initialize
    
    lperiod.sort()
    wlindex=np.array([],dtype=int)
    
    # loop on periods

    for period in lperiod:
        
        # make actual extraction - case data_xyz now handled
        start_date_period = period[0]
        end_date_period   = period[1]

        lindex = np.where( (new_gts.data_xyz[:,0] >= start_date_period) & \
                           (new_gts.data_xyz[:,0] <= end_date_period ) )
        
        wlindex=np.append(lindex, wlindex)
    # end loop periods
        
        
    new_data_xyz   = np.delete( new_gts.data_xyz , wlindex , axis=0 )
    new_sigma = np.delete( new_gts.data , wlindex , axis=0)[:,4:]
    
    # case no observation in periods
    if new_data_xyz.shape[0] == 0:
        if verbose: 
            print("!!! ",self.code," has no data for the selected list of periods ",lperiod)
        if in_place:
            self.data     = None
            self.data_xyz = None
            return( self )
        else:
            new_gts.data = None
            new_gts.data_xyz = None
            return(new_gts)

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

