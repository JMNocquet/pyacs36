###################################################################
def extract_periods(self,lperiod,in_place=False,verbose=False, no_reset=False, ignore_data_xyz=False):
###################################################################
    """
    Extract periods of a Gts.

    Parameters
    ----------
    lperiod : list
        A list [start_date,end_date] or a list of periods e.g. [[2000.1,2003.5],[2009.3,2010.8]].
    in_place : bool, optional
        If True, will make change in place; if False, returns a new time series.
    verbose : bool, optional
        Verbose mode.
    no_reset : bool, optional
        If True, do not reset X0, Y0, Z0 to first epoch of extracted data.
    ignore_data_xyz : bool, optional
        If True, work on .data only and ignore .data_xyz.

    Notes
    -----
    1. X0, Y0, Z0 attributes will be changed if necessary.
    2. Handles both .data and .data_xyz.
    """

    # import 
    import inspect
    import numpy as np
    import pyacs.gts.Gts
    import pyacs.lib.utils

    import logging
    import pyacs.message.message as MESSAGE
    import pyacs.message.verbose_message as VERBOSE
    import pyacs.message.error as ERROR
    import pyacs.message.warning as WARNING
    import pyacs.message.debug_message as DEBUG

    import pyacs.lib.astrotime as at


    # ensure lperiod is a list of lists
    lperiod = pyacs.lib.utils.__ensure_list_of_list(lperiod)

    # check lperiod is not [[]]
    if lperiod==[[]]:
        if in_place:
            return(self)
        else:
            return( self.copy() )

    # CASE 1: ignore_data_xyz is True
    if ignore_data_xyz:
        # check if self has attribute data_xyz
        if hasattr(self, 'data_xyz'):
            # remove data_xyz
            delattr(self, 'data_xyz')
        # converts lperiod into lperiods_seconds
        lperiod_seconds = []
        for period in lperiod:
            start_sec = at.decyear2seconds(period[0])
            end_sec = at.decyear2seconds(period[1])
            lperiod_seconds.append([start_sec,end_sec])

        # convert time series dates into seconds
        t_seconds = at.decyear2seconds(self.data[:,0])

        # loop on lperiod_seconds and extract the data
        all_indices = []
        for period in lperiod_seconds:
            # get the indices of the data in the period
            indices = np.where( (t_seconds >= period[0]) & (t_seconds <= period[1]) )
            # append indices successively
            all_indices.extend(indices[0])
        
        # raise error if the number of data is 0
        if len(all_indices) == 0:
            ERROR("No data found for %s in the selected list of periods %s " % (self.code, str(lperiod)))
            data = None
        else:
            # extract the data
            data = self.data[all_indices]
        
        # create a new gts
        new_gts = self.copy()
        # add the data to the new gts
        new_gts.data = data
        # return
        return(new_gts)

    # CASE 2: ignore_data_xyz is False

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
            ERROR( error )
            return( self )

    
    # initialize
    
    lperiod.sort()
    new_data_xyz=None
    new_sigma = None
    
    # loop on periods

    for period in lperiod:
        
        # make actual extraction - case data_xyz now handled
        start_date_period = period[0]
        end_date_period   = period[1]

        lindex = np.where( (new_gts.data_xyz[:,0] >= start_date_period) & \
                           (new_gts.data_xyz[:,0] <= end_date_period ) )

        sel_data_xyz   = np.copy(new_gts.data_xyz[ lindex ])
        sel_data_sigma = new_gts.data[ lindex  ][:,4:]

        # populates new_data_xyz
        if not isinstance(new_data_xyz, np.ndarray):
            new_data_xyz = sel_data_xyz
            new_sigma = sel_data_sigma
        else:
            new_data_xyz = np.vstack((new_data_xyz,sel_data_xyz))
            new_sigma = np.vstack((new_sigma,sel_data_sigma))

    # end loop periods
    
    # case no observation in periods
    if new_data_xyz.shape[0] == 0:
        VERBOSE(" %s has no data for the selected list of periods %s " % (self.code , str(lperiod)))

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

    #ldate_outliers=self.data[:,0][self.outliers]
    #lupdated_outliers=pyacs.gts.Gts.get_index_from_dates(ldate_outliers, new_data_xyz, tol=0.05)
    
    # handles offsets_date
    
    upd_offsets=[]
    for offset_date in self.offsets_dates:
        if offset_date>=new_data_xyz[0,0] and offset_date<=new_data_xyz[-1,0]:
            upd_offsets.append(offset_date)
    
    # handles X0,Y0,Z0
    if not no_reset:
        new_gts.X0 = sel_data_xyz[0,1]
        new_gts.Y0 = sel_data_xyz[0,2]
        new_gts.Z0 = sel_data_xyz[0,3]
    else:
        new_gts.X0 = self.X0
        new_gts.Y0 = self.Y0
        new_gts.Z0 = self.Z0

# re-generate NEU time series
    new_gts.xyz2neu(corr=True, ref_xyz = [new_gts.X0,new_gts.Y0,new_gts.Z0])

    # re-populate the uncertainties columns
    new_gts.data[:,4:] = new_sigma
    
    # offsets & outliers
        
    new_gts.offsets_dates=upd_offsets
    #new_gts.outliers=lupdated_outliers
    
    if in_place:
        self = new_gts
        return(self)
    else:
        return(new_gts)
