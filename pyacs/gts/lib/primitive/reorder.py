###################################################################
def reorder(self,verbose=False):
###################################################################
    """
    reorder data and/or data_xyz by increasing dates
    always in place
    
    :param verbose: verbose mode
    """
    
    # import
    import inspect
    import numpy as np
    import pyacs.gts.Gts
    from pyacs.gts.lib.errors import GtsCDataError
    
    # check data
    try:
        if not self.cdata(verbose=verbose):
            raise GtsCDataError(inspect.stack()[0][3],__name__,self)
    except GtsCDataError as error:
        print( error )
        return( self )

        
    # test whether reorder is required
    
    if self.data is not None:
        if self.data.shape[0]>1:
            diff = np.diff(self.data[:,0])
            
            if np.min(diff) < 0:
                if verbose:
                    print('-- data needs to be re-ordered')

            else:
                if verbose:
                    print('-- data does not need to be re-ordered')
                return(self)

        else:
            if verbose:
                print('-- data does not need to be re-ordered')
            
            return(self)
            

    
    else:
        if self.data_xyz.shape[0]>1:
            diff = np.diff(self.data_xyz[:,0])
            
            if np.min(diff) < 0:
                if verbose:
                    print('-- data needs to be re-ordered')

            else:
                if verbose:
                    print('-- data does not need to be re-ordered')
                return(self)

        else:
            if verbose:
                print('-- data does not need to be re-ordered')
            return(self)
        
        
    # save outliers dates before re-ordering - outliers are assumed for .data
    if ( self.outliers != [] ) and ( self.data is not None ):
        ldate_outliers=self.data[:,0][self.outliers]
    else:
        ldate_outliers = []
    
    # reorder
    if isinstance(self.data,np.ndarray):
        lindex=np.argsort(self.data[:,0])
        self.data=self.data[lindex,:]
        
    if isinstance(self.data_xyz,np.ndarray):
        lindex=np.argsort(self.data_xyz[:,0])
        self.data_xyz=self.data_xyz[lindex,:]
    
    # check data
    if not self.cdata(verbose=verbose):
        print('!!! ERROR: .data found to be incorrect by .cdata after reordering. Check time series and results.')
        return(self)

    if verbose:
        print('-- data order OK')
    
    # deal with outliers

    if verbose:
        print('-- updating outliers list for ',self.code)
    
    lupdated_outliers=pyacs.gts.Gts.get_index_from_dates(ldate_outliers, self.data, tol=0.1)
    self.outliers=lupdated_outliers
    
    return(self)
